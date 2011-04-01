#!/usr/local/bin/perl

##################################################################
# Description: HTAB.pm
#
# Parses hmmer 3 htab files and fetches HMM annotation information
# from a SQLite database (hmm3.db)
# --------------
# Date:   Dec 23, 2010  
##################################################################

package CAMERA::Parser::HTAB;

use strict;
use warnings;

use DBI;
use Carp;
use CAMERA::Polypeptide;
use CAMERA::PolypeptideSet;
use CAMERA::AnnotationData::Polypeptide;
use DBM::Deep;
use File::Copy;

BEGIN { $SIG{__DIE__}= sub{ Carp::confess @_ } }

my $to_text;

## allows text to be directly written out to a filehandle in text mode
my $text_fh = undef;

sub new {
    my ($class, %args) = @_;

    my $self = {};

    ## path to sqlite database
    my $sqliteDB = "$args{'work_dir'}/hmm3.db";

	## open sqlite database
	$self->{'db'} = DBI->connect( "dbi:SQLite:$sqliteDB", "", "", {PrintError=>1,RaiseError=>1,AutoCommit=>0} );
	if ( !defined $self->{'db'} ) {
		die "could not connect to sqlite database: $sqliteDB" . $DBI::errstr;
	}
    
    ## pass in a reference to a CAMERA::PolypeptideSet object
    ## can be null if no other annotation has been parsed yet
    if (ref($args{'polypeptides'}) eq 'CAMERA::PolypeptideSet') {
        $self->{'polypeptides'} = $args{'polypeptides'};
    } else {
        $self->{'polypeptides'} = new CAMERA::PolypeptideSet();
    }

    ## initialize the hmm data hash
    ##_initialize_hmm_data();
    
    return bless $self, $class;
}

sub parse {	
    my ($self, $file, $text_mode, $txt_fh) = @_;

    if ($text_mode) {
        $to_text = 1;
    }
    if ($txt_fh) {
        $text_fh = $txt_fh;
    }

    my $rank = {};

    my %has_trusted_cutoffs = (
                        'TIGRFAM::FullLength::Equivalog'                => 1,
                        'TIGRFAM::FullLength::HypotheticalEquivalog'    => 1,
                        'TIGRFAM::FullLength::Exception'                => 1,
                        'TIGRFAM::FullLength::Subfamily'                => 1,
                        'TIGRFAM::FullLength::Superfamily'              => 1,
                        'TIGRFAM::FullLength::EquivalogDomain'          => 1,
                        'TIGRFAM::FullLength::HypotheticalEquivalogDomain'  => 1,
                        'TIGRFAM::FullLength::SubfamilyDomain'          => 1,
                        'TIGRFAM::FullLength::Domain'                   => 1,
                        'PFAM::FullLength::Equivalog'                   => 1,
                        'PFAM::FullLength::HypotheticalEquivalog'       => 1,
                        'PFAM::FullLength::Subfamily'                   => 1,
                        'PFAM::FullLength::Superfamily'                 => 1,
                        'PFAM::FullLength::EquivalogDomain'             => 1,
                        'PFAM::FullLength::HypotheticalEquivalogDomain' => 1,
                        'PFAM::FullLength::SubfamilyDomain'             => 1,
                        'PFAM::FullLength::Domain'                      => 1,
                        'PFAM::FullLength::Uncategorized'               => 1,
			'TIGRFAM::FRAG::HypotheticalEquivalog'         => 1,
			'TIGRFAM::FRAG::Exception'                     => 1,
			'TIGRFAM::FRAG::Subfamily'                     => 1 ,
                        'TIGRFAM::FRAG::Superfamily'                    => 1,
                        'TIGRFAM::FRAG::EquivalogDomain'               => 1,
                        'TIGRFAM::FRAG::HypotheticalEquivalogDomain'   => 1,
                        'TIGRFAM::FRAG::SubfamilyDomain'                => 1,
                        'TIGRFAM::FRAG::Domain'                        => 1,
			'PFAM::FRAG::Equivalog'                 => 1,
                        'PFAM::FRAG::HypotheticalEquivalog'             => 1,
                        'PFAM::FRAG::Subfamily'                         => 1,
                        'PFAM::FRAG::Superfamily'                       => 1,            
                        'PFAM::FRAG::EquivalogDomain'                   => 1,
                        'PFAM::FRAG::HypotheticalEquivalogDomain'       => 1,
                        'PFAM::FRAG::SubfamilyDomain'                  => 1,
                        'PFAM::FRAG::Domain'                          => 1,
                        'PFAM::FRAG::Uncategorized'                     => 1,
                              );
    
    my $infh = _open_file_read($file);
  
    my $counter;
    my $result = '';
    while (<$infh>) {
        chomp;
		print "Reading htab line=$_\n";
        
        ## fix for concatenated output not from htab.pl -m
        if (/^No hits/ || /^No domain/) {
            next;
        }
        
        my @t = split("\t", $_);
        my $hmm_id = $t[0];
        
        ## database results
        my $hmm3 = $self->get($hmm_id);
        
        my $pep_id = $t[5];
        
        unless ($pep_id) { 
            confess "failed parsing htab file";
        }

        $rank->{$pep_id}++;
        
        my $pep = undef;
        if ($self->{'polypeptides'}->exists($pep_id)) {
            $pep = $self->{'polypeptides'}->get($pep_id);
        } else {
            $pep = new CAMERA::Polypeptide('id' => $pep_id);
            $self->{'polypeptides'}->add($pep);
        }
        
        my $hmm_type = _get_hmm_type($hmm3);
    		
        ## if we don't recognize the hmm_type then skip
        unless ($hmm_type) {
            next;
        }

		print "Adding annotation for hmm_id=$hmm_id pep_id=$pep_id\n";
        
        $result .= _add_annotation($pep, $hmm3, $hmm_type, $rank->{$pep_id});
   }
   return $result;
}

sub _open_file_read {
    my ($file) = @_;

    open (my $infh, $file) || die "Failed to open file '$file' for reading: $!";

    return $infh;
}

sub _add_annotation {
    my ($pep, $hmm3, $type, $rank) = @_;
    
	my $hmm_id = $hmm3->{hmm_acc};

    my $annotation = new CAMERA::AnnotationData::Polypeptide(
                'id'        => $pep->{'id'},
                'type'      => $type,
                'source'    => $hmm_id,
                'rank'      => $rank,
                                           );
                                               
    my $success = 0;
                               
    $success += _set_common_name($annotation,$hmm3); 
    $success += _set_gene_symbol($annotation,$hmm3);
    $success += _set_ec_numbers($annotation,$hmm3);
    $success += _set_go_ids($annotation,$hmm3);
     
    if ($success) {
        if ($to_text) {
            if ($text_fh) {
                print $text_fh $annotation->to_string();
                return '';
            } else {
                return $annotation->to_string();
            }
        } else {
            $pep->add_annotation($annotation);
            return '';
        }
    }
    
}

sub _get_hmm_type {
    my $hmm3 = shift;

	my $hmm_id = $hmm3->{hmm_acc};

    my $hmm_lib = '';

    if ($hmm_id =~ /^TIGR/) {
        $hmm_lib = 'TIGRFAM';
    } elsif ($hmm_id =~ /^PF/) {
        $hmm_lib = 'PFAM'; 
    } else {
        confess "Unrecognized HMM type for '$hmm_id'";
    }

    my $hmm_length = 'FullLength';

    my $iso_type = $hmm3->{iso_type};
    my $cutoff 	 = $hmm3->{trusted_cutoff};
   
    if (! $iso_type) {
        confess "No iso_type for '$hmm_id'";
    }
    
    my $type = $hmm_lib.'::'.$hmm_length.'::';
	
    if ($iso_type =~ /^(equivalog)$|^(PFAM_equivalog)$/) {
        $type .= 'Equivalog';
    } elsif ($iso_type =~ /^(hypoth_equivalog)$/) {
        $type .= 'HypotheticalEquivalog';
    } elsif ($iso_type =~ /^(exception)$/) {
        $type .= 'Exception';
    } elsif ($iso_type =~ /^(subfamily)$/) {
        $type .= 'Subfamily';
    } elsif ($iso_type =~ /^(superfamily)$/) {
        $type .= 'Superfamily';
    } elsif ($iso_type =~ /^(equivalog_domain)$|^(PFAM_equivalog_domain)$/) {
        $type .= 'EquivalogDomain';
    } elsif ($iso_type =~ /^(hypoth_equivalog_domain)$/) {
        $type .= 'HypotheticalEquivalogDomain';
    } elsif ($iso_type =~ /^(subfamily_domain)$/) {
        $type .= 'SubfamilyDomain';
    } elsif ($iso_type =~ /^(domain)$/) {
        $type .= 'Domain';
    } elsif ($iso_type =~ /^(PFAM)$/) {
        $type .= 'Uncategorized';
    } else {
        $type = '';
    }
    
    return $type;
}

sub _set_common_name {
    my ($annotation,$hmm3) = @_;
    
    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($hmm3->{com_name}) {
		return 0;
    }

    my $common_name = $hmm3->{com_name};

    $annotation->add_attribute('common_name', $common_name);
    
    return 1;
}

sub _set_gene_symbol {
    my ($annotation,$hmm3) = @_;
    
    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($hmm3->{gene_sym}) {
		return 0;
    }
   
   	my $geneSymbol = $hmm3->{gene_sym};
    
  	
	$annotation->add_attribute('gene_symbol', $geneSymbol);

	return 1;  
}

sub _set_ec_numbers {
    my ($annotation,$hmm3) = @_;
    
    my $hmm_id = $hmm3->{hmm_acc};
 
    unless ($hmm3->{ec_num}) {
		return 0;
    }
       
    my $ec_num = $hmm3->{ec_num};

	my @ec_numbers = split(' ',$ec_num);
	
	foreach my $ec(@ec_numbers) {
	    $annotation->add_attribute('EC', $ec);
	}
	
	return 1; 
}

sub _set_go_ids {
    my ($annotation,$hmm3) = @_;

    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless($hmm3->{go}) {
		return 0;
    }
    
	my @go_ids = @{$hmm3->{go}};
    
	foreach my $go(@go_ids) {   
	    $annotation->add_attribute('GO', @{$go}[0]);
	}

	return 1;
    
}

sub get() {
	my ($self,$accession) = @_;
	
	my $db_proc = $self->{'db'};
	
	my $hmm3Result = $self->{'db'}->selectrow_hashref(
			"select hmm_acc, hmm_len, trusted_cutoff, noise_cutoff, gathering_cutoff, hmm_com_name as com_name, "
				. "gene_sym, ec_num, iso_type, trusted_cutoff2 "
				. "from hmm3 "
				. "where hmm_acc = ? and is_current = 1",
			undef, $accession );
		
	## custom PFAM iso-type mapping provided by Dan Haft	
	if($hmm3Result->{hmm_acc} =~ /^PF/) {
		if($hmm3Result->{iso_type} =~ /^Domain$/) {
			$hmm3Result->{iso_type} = 'domain';
		}
		elsif($hmm3Result->{iso_type} =~ /^Motif$/) {
			$hmm3Result->{iso_type} = 'domain';
		}
		elsif($hmm3Result->{iso_type} =~ /^Family$/) {
			$hmm3Result->{iso_type} = 'PFAM';
		}		
	}	
			
			
	my $goResults = $self->{'db'}->selectall_arrayref(
			"select go_term from hmm_go_link where hmm_acc = ? ",
			undef, $accession );	
	

	if($goResults) {
		$hmm3Result->{go} =$goResults;
	}		
	
	printHmm($hmm3Result);
	
	return $hmm3Result;
}

sub printHmm() {
	my ($hmm3) = @_;
	
	print "accession:\t".$hmm3->{hmm_acc}."\n";
	print "common name:\t".$hmm3->{com_name}."\n";
	
	if(defined $hmm3->{gene_sym} ) {
		print "gene symbol:\t".$hmm3->{gene_sym}."\n";
	}

	if(defined $hmm3->{ec}) {		
		my $ec_num = $hmm3->{ec_num};

		my @ec_numbers = split(' ',$ec_num);
	 		
		foreach my $ec(@ec_numbers) {
		   print "ec:\t".$ec."\n";
		}
	}
	if(defined $hmm3->{go}) {
		my @go_ids = @{$hmm3->{go}};	    
		foreach my $go(@go_ids) {	   
		    print "go:\t".@{$go}[0]."\n";
		}	
	}	
}
1;