#!/usr/local/bin/perl

package CAMERA::Parser::BTAB;

use strict;
use warnings;

#use lib "/share/apps/ergatis/current/lib";
use CAMERA::Polypeptide;
use CAMERA::PolypeptideSet;
use CAMERA::AnnotationData::Polypeptide;
use Split_DB_File;
use Carp;
use File::Copy;
use File::Basename;
use UniRef::UniRefDao;

my $to_text;
my $text_fh = undef;

sub new {
    my ($class, %args) = @_;

    my $self = {};
    
    ## pass in a reference to a CAMERA::PolypeptideSet object
    ## can be null if no other annotation has been parsed yet
    if (ref($args{'polypeptides'}) eq 'CAMERA::PolypeptideSet') {
        $self->{'polypeptides'} = $args{'polypeptides'};
    } else {
        $self->{'polypeptides'} = new CAMERA::PolypeptideSet();
    }
    
    ## path to sqlite database
    my $sqliteDB = "$args{'work_dir'}/uniref100.db";
   
   	$self->{'unirefDao'} = new UniRef::UniRefDao($sqliteDB );

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
        
        ## parse Uniref100 header
        $description =~  /^(.*)\sn=\d*\sTax=(.*)\sRepID=.*/; 
      	my $common_name= $1;
      	my $taxon_name = $2;
      	       	    
      	## stripping off unire100 prefix for datbase lookup  
      	$subject_id =~ s/UniRef100_//;
  
        ## fetch uniref lookup information
        my $unirefEntry = $self->{'unirefDao'}->get($subject_id);   
        
        ## fetch information if entry id reviewed or not
        my $isReviewed = $unirefEntry->isReviewed();      
              
        ## check the hit type 
        my $hit_type = _get_hit_type($isReviewed, $percent_identity, $percent_coverage);
      
        ## skip if we don't recognize the type
        unless ($hit_type) {
            next;
        }
           
        ## skip if we've already stored annotation of this type
        if ($annotation_counter->{$pep_id}->{$hit_type}) {
            next;
        }
            
        my $pep = undef;
        
        ## add peptide to polypeptide set if it does not exist
        if ($self->{'polypeptides'}->exists($pep_id)) {
            $pep = $self->{'polypeptides'}->get($pep_id);
        } else {
            $pep = new CAMERA::Polypeptide('id' => $pep_id);
            $self->{'polypeptides'}->add($pep);
        }

        $result .= _add_annotation($pep, "UniRef100_$subject_id", $hit_type, $rank->{$pep_id}->{$database_name}, $unirefEntry,$common_name,$taxon_name);
        $annotation_counter->{$pep_id}->{$hit_type}++;

    }
    
    ## close store
    $self->{'unirefDao'}->close();
    
    return $result;
}

sub _get_hit_type {
    my ($isReviewed, $pct_id, $pct_cov) = @_;

	if ($pct_id >= 35 && $pct_cov >= 80) {
		if($isReviewed) {
        	return "UnirefBLASTP::Reviewed";
		}
		else {
			return "UnirefBLASTP::HighConfidence";
		}
    } elsif ($pct_id < 35 && $pct_cov >= 80) {
        return "UnirefBLASTP::Putative";
    } elsif ($pct_id >= 35 && $pct_cov < 80) {
        return "UnirefBLASTP::ConservedDomain";
    } 
    else {
    	return "UnirefBLASTP::LowConfidence";
    }   
}

sub _add_annotation {
    my ($pep, $subject_id, $type, $rank, $unirefEntry,$common_name,$taxon_name) = @_;    
    
    my $annotation = new CAMERA::AnnotationData::Polypeptide(
                                'id'        => $pep->{'id'},
                                'type'      => $type,
                                'source'    => $subject_id,
                                'rank'      => $rank,);
                                
    ## add UniRef100 annotations fetched from the defline                            
    $annotation->add_attribute('common_name', $common_name);  
    $annotation->add_attribute('TAXON', $taxon_name);
           
    ## add UniRef100 annotations fetched from the SQLite database                                                      
	if(defined $unirefEntry) {   
		
		## fetch fields from UniRef object
	    my $geneSymbol  = $unirefEntry->getGeneSymbol();
	     
	    my @goIds	= @{$unirefEntry->getGoIds()};
	    my @ecIds	= @{$unirefEntry->getEcIds()};
	    my @cazyIds	= @{$unirefEntry->getCazyIds()};
    		
    	## add gene symbol annotation	
    	$annotation->add_attribute('gene_symbol', $geneSymbol);
    
	    ## add GO annotations
	    foreach my $goId( @goIds) {
	    	$annotation->add_attribute('GO', $goId);
	    }
	    ## add EC annotations
	    foreach my $ecId( @ecIds) {
	    	$annotation->add_attribute('EC', $ecId);
	    } 
	    ## add CAZY annotations
	    foreach my $cazyId( @cazyIds) {
	    	$annotation->add_attribute('CAZY', $cazyId);
	    } 
	}
         
    if ($unirefEntry->hasAnnotation() || $common_name || $taxon_name) {
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

sub _open_file_read {
    my ($file) = @_;

    open (my $infh, $file) || die "Failed to open file '$file' for reading: $!";

    return $infh;
}
1;