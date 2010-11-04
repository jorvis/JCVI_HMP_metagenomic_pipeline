#!/usr/local/bin/perl

package CAMERA::Parser::BTAB;

use strict;
use warnings;

use CAMERA::Polypeptide;
use CAMERA::PolypeptideSet;
use CAMERA::AnnotationData::Polypeptide;
use Split_DB_File;
use Carp;
use File::Copy;
use File::Basename;

## source directory containing parsed panda headers
my $pandaHeaderDir = "../../../data/panda-header-index";

## Panda header data serialized with Split_DB_File
my $panda_header_file;
my $panda_header_offset_db;

my %pandaHeaderOffsets;


my $to_text;
my $text_fh = undef;

sub new {
    my ($class, %args) = @_;

    my $self = {};
    
    my $panda_header_dir = ($args{'work_dir'}) ? $args{'work_dir'} : '';

    $panda_header_file = "$panda_header_dir/parsedPanda_header";
    $panda_header_offset_db = "$panda_header_dir/parsedPanda";

    if (! -e $panda_header_file) {
        my @files = glob($pandaHeaderDir."/parsedPanda*");
        foreach my $file(@files) {
            copy($file, $panda_header_dir."/".basename($file));
        }
    }

    unless (-e $panda_header_file) {
        confess "parsedPanda_header file does not exist in dir '$panda_header_dir'";
    }

    ## pass in a reference to a CAMERA::PolypeptideSet object
    ## can be null if no other annotation has been parsed yet
    if (ref($args{'polypeptides'}) eq 'CAMERA::PolypeptideSet') {
        $self->{'polypeptides'} = $args{'polypeptides'};
    } else {
        $self->{'polypeptides'} = new CAMERA::PolypeptideSet();
    }
    
    tie(%pandaHeaderOffsets, 'Split_DB_File', $panda_header_offset_db)
        or die("Unable to tie Panda header offset hash");

    return bless $self, $class;
}

sub parse {
    my ($self, $file, $text_mode, $txt_fh) = @_;

    if ($text_mode) {
        $to_text = 1;
    }
    $text_fh = $txt_fh;

    my $infh = _open_file_read($file);
 
    my $rank = {};

    my $annotation_counter = {};

    my $result = ''; 
    while (<$infh>) {
        chomp;
        my (
            $pep_id, 
            $analysis_date, 
            $query_length, 
            $search_method,
            $database_name, 
            $subject_id, 
            $query_start,
            $query_end,
            $subject_start,
            $subject_end,
            $percent_identity,
            $percent_similarity,
            $score,
            $file_offset1,
            $file_offset2,
            $description,
            $frame,
            $query_strand,
            $subject_length,
            $expect,
            $pvalue,
           ) = split("\t", $_);
        
        unless ($pep_id) { 
            confess "failed parsing btab file";
        }

        $rank->{$pep_id}->{$database_name}++;
        
        ## swap start and end positions if necessary
        if ($query_start > $query_end) {
            ($query_start, $query_end) = ($query_end, $query_start);
        }
        ## swap start and end positions if necessary
        if ($subject_start > $subject_end) {
            ($subject_start, $subject_end) = ($subject_end, $subject_start);
        }
        
        my $query_hit_length = $query_end - $query_start;
        my $percent_coverage = sprintf("%.2f", $query_hit_length / $query_length * 100);
        
        ## need create a sub is_characterized that uses a dump from the CHAR database
        #my $characterized = _is_characterized($pep_id);
        my $characterized = 0;            
       
        ## check the hit type 
        my $hit_type = _get_hit_type($characterized, $percent_identity, $percent_coverage);
      
        ## skip if we don't recognize the type
        unless ($hit_type) {
            next;
        }
           
        ## skip if we've already stored annotation of this type
        if ($annotation_counter->{$pep_id}->{$hit_type}) {
            next;
        }
        
        ## fetch a panda header array ref for the subject
        my $panda_header_arr_ref = _get_panda_header($subject_id);
      
        ## to make the header values more clear
        my $panda = {
                      'accession'   =>  $panda_header_arr_ref->[0],
                      'common_name' =>  $panda_header_arr_ref->[1],
                      'ec_numbers'  =>  $panda_header_arr_ref->[2],
                      'gene_symbol' =>  $panda_header_arr_ref->[3],
                      'species'     =>  $panda_header_arr_ref->[4],
                      'go_ids'      =>  $panda_header_arr_ref->[5],
                      'TIGR_roles'  =>  $panda_header_arr_ref->[6],
                      'exp'         =>  $panda_header_arr_ref->[7],
                      'wgp'         =>  $panda_header_arr_ref->[8],
                      'cg'          =>  $panda_header_arr_ref->[9],
                    };
        
        my $pep = undef;
        if ($self->{'polypeptides'}->exists($pep_id)) {
            $pep = $self->{'polypeptides'}->get($pep_id);
        } else {
            $pep = new CAMERA::Polypeptide('id' => $pep_id);
            $self->{'polypeptides'}->add($pep);
        }

        $result .= _add_annotation($pep, $subject_id, $hit_type, $rank->{$pep_id}->{$database_name}, $panda);
        $annotation_counter->{$pep_id}->{$hit_type}++;

    }
    return $result;
}

sub _get_hit_type {
    my ($characterized, $pct_id, $pct_cov) = @_;

    if ($characterized && $pct_id >= 35 && $pct_cov >= 80) {
        "PandaBLASTP::Characterized";
    } elsif ($pct_id >= 35 && $pct_cov >= 80) {
        "PandaBLASTP::HighConfidence";
    } elsif ($pct_id < 35 && $pct_cov >= 80) {
        "PandaBLASTP::Putative";
    } elsif ($pct_id >= 35 && $pct_cov < 80) {
        "PandaBLASTP::ConservedDomain";
    }
}

sub _add_annotation {
    my ($pep, $subject_id, $type, $rank, $panda) = @_;    
    
    my $annotation = new CAMERA::AnnotationData::Polypeptide(
                                'id'        => $pep->{'id'},
                                'type'      => $type,
                                'source'    => $subject_id,
                                'rank'      => $rank,
                                                            );
                                                           
    my $success = 0;

    $success += _set_common_name($annotation, $panda);
    $success += _set_gene_symbol($annotation, $panda);
    $success += _set_ec_numbers($annotation, $panda);
    $success += _set_go_ids($annotation, $panda);
    $success += _set_TIGR_roles($annotation, $panda);

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

sub _set_common_name {
    my ($annotation, $panda) = @_;

    my $source_id = $annotation->{'source'};
    my $type = $annotation->{'type'};

    unless ($panda->{'common_name'}) {
        return 0;
    }

    my $common_name = $panda->{'common_name'};

    ## remove taxon:TAXON_ID from end of common_name
    $common_name =~ s/ taxon:\d+//;
    
    $annotation->add_attribute('common_name', $common_name);

    return 1;
}

sub _set_gene_symbol {
    my ($annotation, $panda) = @_;

    my $source_id = $annotation->{'source'};
    my $type = $annotation->{'type'};
    
    unless ($panda->{'gene_symbol'}) {
        return 0;
    }

    my $gene_symbol = $panda->{'gene_symbol'};

    $annotation->add_attribute('gene_symbol', $gene_symbol);
    
    return 1;
    
}

sub _set_ec_numbers {
    my ($annotation, $panda) = @_;

    my $source_id = $annotation->{'source'};
    my $type = $annotation->{'type'};
    
    unless ($panda->{'ec_numbers'}) {
        return 0;
    }

    my $ec_numbers = $panda->{'ec_numbers'};

    my @ec = split(" ", $ec_numbers);

    foreach my $ec_number(@ec) {
        if ($ec_number ne '') {
            $annotation->add_attribute('EC', $ec_number);
        }
    }
    
    return 1;
}

sub _set_go_ids {
    my ($annotation, $panda) = @_;

    my $source_id = $annotation->{'source'};
    my $type = $annotation->{'type'};
    
    unless ($panda->{'go_ids'}) {
        return 0;
    }

    my $go_ids = $panda->{'go_ids'};

    my @go = split(" ", $go_ids);
    
    foreach my $go_id(@go) {
        if ($go_id ne '') {
            $annotation->add_attribute('GO', $go_id);
        }
    }
    
    return 1;
}

sub _set_TIGR_roles {
    my ($annotation, $panda) = @_;

    my $source_id = $annotation->{'source'};
    my $type = $annotation->{'type'};
    
    unless ($panda->{'TIGR_roles'}) {
        return 0;
    }

    my $TIGR_roles = $panda->{'TIGR_roles'};

    my @roles = split(" ", $TIGR_roles);

    foreach my $role(@roles) {
        if ($role ne '') {
            $annotation->add_attribute('TIGR_role', $role);
        }
    }
    
    return 1;
}

sub _open_file_read {
    my ($file) = @_;

    open (my $infh, $file) || die "Failed to open file '$file' for reading: $!";

    return $infh;
}

sub _get_panda_header {
    my ($accession) = @_;

    my $result = $pandaHeaderOffsets{$accession};

    unless($result) {
        confess("Could not find Panda header for accession '$accession'");
    }

    open(my $infh, $panda_header_file) || die("Couldn't open '$panda_header_file': $!");

    my ($pos, $len) = split(/\s+/, $result);
    
    seek($infh, $pos, 0);
    read($infh, my $values, $len);

    $" = "\n";
    my @vals = split(/\t/, $values);
   
    for (my $i = 0; $i < scalar(@vals); $i++) {
        $vals[$i] =~ s/^\s+|\s+$//g;
    }
    
    return \@vals;
}

1;
