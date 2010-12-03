#!/usr/local/bin/perl

package CAMERA::Parser::HTAB;

use strict;
use warnings;

use Carp;
use CAMERA::Polypeptide;
use CAMERA::PolypeptideSet;
use CAMERA::AnnotationData::Polypeptide;
use DBM::Deep;
use File::Copy;

BEGIN { $SIG{__DIE__}= sub{ Carp::confess @_ } }

#my $hmm_db_source = '../../../data';
my $hmm_db_source = '/usr/local/projects/jorvis/metagenomic_annotation/SOM/data/hmm-index';
my $hmm_db_name = 'ALL_LIB.HMM.dump';




my $hmm_db_source_file = $hmm_db_source."/".$hmm_db_name;

print "Using dbm souce file=" . $hmm_db_source_file . "\n";

my $hmm_db_file = '';

my $hmm_data = {};

my $hmm_frag = 0;

my $to_text;

## allows text to be directly written out to a filehandle in text mode
my $text_fh = undef;

sub new {
    my ($class, %args) = @_;

    my $self = {};

    my $hmm_db_dir = ($args{'work_dir'}) ? $args{'work_dir'} : '';

    $hmm_db_file = "$hmm_db_dir/$hmm_db_name";
    print "$hmm_db_source_file\n";
    if (! -e $hmm_db_file) {
        copy($hmm_db_source_file, $hmm_db_file);
    }

    unless (-e $hmm_db_file) {
        confess "HMM data file '$hmm_db_file' does not exist";
    }

    if ($args{'frag'}) {
        $hmm_frag = 1;
    }
    ## pass in a reference to a CAMERA::PolypeptideSet object
    ## can be null if no other annotation has been parsed yet
    if (ref($args{'polypeptides'}) eq 'CAMERA::PolypeptideSet') {
        $self->{'polypeptides'} = $args{'polypeptides'};
    } else {
        $self->{'polypeptides'} = new CAMERA::PolypeptideSet();
    }

    ## initialize the hmm data hash
    _initialize_hmm_data();
    
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
        if (/^No hits/) {
            next;
        }
        my @t = split("\t", $_);
        my $hmm_id = $t[0];
        my $pep_id = $t[5];
        my $score = $t[11];
        
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
        
        my $hmm_type = _get_hmm_type($hmm_id);
   
        ## if we don't recognize the hmm_type then skip
        unless ($hmm_type) {
            next;
        }
        
        ## filter on trusted cutoffs
        if ($has_trusted_cutoffs{$hmm_type}) {
            my $cutoff = $hmm_data->{$hmm_id}->[2];
	    if($hmm_frag){
		if($hmm_type =~ /PFAM/){
#		    $cutoff = $hmm_data->{$hmm_id}->[9];
		}
	    }
            unless ($score >= $cutoff) {
		print "Skipping because score=$score is less than cutoff=$cutoff\n";
                next;
            }
        }

	print "Adding annotation for hmm_id=$hmm_id pep_id=$pep_id\n";
        
        $result .= _add_annotation($pep, $hmm_id, $hmm_type, $rank->{$pep_id});
   }
   return $result;
}

sub _open_file_read {
    my ($file) = @_;

    open (my $infh, $file) || die "Failed to open file '$file' for reading: $!";

    return $infh;
}

sub _add_annotation {
    my ($pep, $hmm_id, $type, $rank) = @_;

    my $annotation = new CAMERA::AnnotationData::Polypeptide(
                'id'        => $pep->{'id'},
                'type'      => $type,
                'source'    => $hmm_id,
                'rank'      => $rank,
                                           );
                                               
    my $success = 0;
                                           
    $success += _set_common_name($annotation); 
    $success += _set_gene_symbol($annotation);
    $success += _set_ec_numbers($annotation);
    $success += _set_go_ids($annotation);
    $success += _set_TIGR_roles($annotation);
    
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
    my ($hmm_id) = @_;

    my $hmm_lib = '';

    if ($hmm_id =~ /^TIGR/) {
        $hmm_lib = 'TIGRFAM';
    } elsif ($hmm_id =~ /^PF/) {
        $hmm_lib = 'PFAM';
    } else {
        confess "Unrecognized HMM type for '$hmm_id'";
    }

    my $hmm_length = ($hmm_frag) ? 'FRAG' : 'FullLength';

    my $hmm_type = $hmm_data->{$hmm_id}->[0];
    my $iso_type = $hmm_data->{$hmm_id}->[1];
    my $cutoff = $hmm_data->{$hmm_id}->[2];
   
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
    my ($annotation) = @_;
    
    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($hmm_data->{$hmm_id}->[4]) {
	return 0;
    }

    my $common_name = $hmm_data->{$hmm_id}->[4];

    $annotation->add_attribute('common_name', $common_name);
    
    return 1;
}

sub _set_gene_symbol {
    my ($annotation) = @_;
    
    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($hmm_data->{$hmm_id}->[5]) {
	return 0;
    }
    my $data_type = ref($hmm_data->{$hmm_id}->[5]);
    
    if(!$data_type){
	print "there is no data type\n";
	#return 1;
    }else{

	my $gc = $hmm_data->{$hmm_id}->[5];
	$annotation->add_attribute('gene_symbol', $gc);

#	my @gene_symbol = @{$hmm_data->{$hmm_id}->[5]};
#
#	foreach my $gc(@gene_symbol) {
#	    $annotation->add_attribute('gene_symbol', $gc);
#	}

	return 1;
    }
}

sub _set_ec_numbers {
    my ($annotation) = @_;
    
    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($hmm_data->{$hmm_id}->[6]) {
	return 0;
    }
    my $data_type = ref($hmm_data->{$hmm_id}->[6]);
    
    if(!$data_type){
	print "there is no data type\n";
	#return 1;
    }else{

	my @ec_numbers = @{$hmm_data->{$hmm_id}->[6]};
	
	foreach my $ec(@ec_numbers) {
	    $annotation->add_attribute('EC', $ec);
	}
	
	return 1;
    }  
}

sub _set_go_ids {
    my ($annotation) = @_;

    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($hmm_data->{$hmm_id}->[7]) {
	return 0;
    }
    my $data_type = ref($hmm_data->{$hmm_id}->[7]);
    if(!$data_type){
	print "htere is no data type\n";

    }else{
	my @go_ids = @{$hmm_data->{$hmm_id}->[7]};
    
	foreach my $go(@go_ids) {
	    print "$hmm_id,$go\n";
	    $annotation->add_attribute('GO', $go);
	}

	return 1;
    }
    
}

sub _set_TIGR_roles {
    my ($annotation) = @_;
    
    my $hmm_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($hmm_data->{$hmm_id}->[8]) {
	return 0;
    }
    my $data_type = ref($hmm_data->{$hmm_id}->[8]);
    if(!$data_type){
	print "there is no data type\n";
	
    }else{
	my @roles = @{$hmm_data->{$hmm_id}->[8]};
    
	foreach my $roles_ref(@roles) {
	    $annotation->add_attribute('TIGR_role', $roles_ref->[0]);
	}

	return 1;
    }
   
}

sub _initialize_hmm_data {
   
    my %hmm_db; 
    if (-e $hmm_db_file) {
	
        # ## the database already exists
        #tie %hmm_db, 'MLDBM', $hmm_db_file, O_RDONLY, 0444
        #  or die "Couldn't tie database file '$hmm_db_file': $!";
        #$hmm_data = \%hmm_db;
        print "Using $hmm_db_file for hmm data\n";
        $hmm_data = new DBM::Deep(file => $hmm_db_file, locking => 0, autoflush => 0);

    } else {

	die "Could not find file $hmm_db_file\n";

    }
}

1;
