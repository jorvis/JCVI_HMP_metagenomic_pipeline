#!/usr/local/bin/perl

##################################################################
# Description: PredictedProtein.pm
# This module is responsible for implemementing the annotation 
# rules for CAMERA predicted proteins
#
# Input to the module is a PolypeptideSet containing Polypeptide
# objects with associated AnnotationData::Polypeptide objects
#
# The attributes contained in the AnnotationData objects are expected
# to have been already processed during the parsing of the source data,
# so the attributes should not need cleanup of any kind.
#
# For example, all leading and trailing whitespace should have already 
# been removed from all attributes; all GO ids should be in the form 
# 'GO:#######'; all EC numbers should be in the form '#.#.#.#'; all 
# multi-value attributes (eg: GO, EC) for a given annotation data type
# should already have been split and should be present as single array
# elements in the attribute arrays; all attempts at consistency of
# capitalization (or other formatting) should have been (and should be) 
# made in the parsing scripts.
# --------------
# Date:   Dec 21, 2010  
##################################################################

package CAMERA::AnnotationRules::PredictedProtein;

use strict;
use warnings;
use Carp;
use CAMERA::AnnotationData::Polypeptide;
use CAMERA::AnnotationRules::Util;
use CAMERA::Parser::SynonymsDao;
use CAMERA::AnnotationData::PolypeptideDataTypes;
use Data::Dumper;

{
	
    ## constructor             
    sub new { 
        my ($class, $synonyms) = @_;
            
        my $self = {
       		_synonymsDao => new CAMERA::Parser::SynonymsDao($synonyms),
   		};
   		
        return bless $self, $class;
    }	
	
    ## ordered array containing all or a subset of the annotation data types
    ## specified in CAMERA::AnnotationData::PolypeptideDataTypes
    ## order should be by order of annotation attribute assignment from 
    ## highest priority to lowest
    ## format is:
    ## ANNOTATION_TYPE|X|attribute1 attribute2 ...
    ## where X is either '=' for assign, '+' for append, '-' for replace attribute
    my @annotation_order = (

                ## equiv. level tigrfam and pfam - full length get top priority
                ## equiv. pfam almost as good as equiv. tigrfam hits
                ## exceptions are the absolute best hits

                'TIGRFAM::FullLength::Exception|=|common_name gene_symbol GO EC',
                'TIGRFAM::FullLength::Equivalog|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::Equivalog|=|common_name gene_symbol GO EC',
                'TIGRFAM::FullLength::HypotheticalEquivalog|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::HypotheticalEquivalog|=|common_name gene_symbol GO EC',

                'TIGRFAM::FRAG::Exception|=|GO',
                'TIGRFAM::FRAG::Equivalog|=|GO',
                'TIGRFAM::FRAG::HypotheticalEquivalog|=|GO',
                
                'UnirefBLASTP::Reviewed|=|GO EC CAZY TAXON',
                
                'TIGRFAM::FullLength::Domain|=|GO',
                
                ## equivalog level hits vs tigrfam frag, only use equivalog level for TIGRFAM
                'TIGRFAM::FRAG::Exception|=|common_name gene_symbol GO EC',
                'TIGRFAM::FRAG::Equivalog|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::Equivalog|=|common_name gene_symbol GO EC',
                'TIGRFAM::FRAG::HypotheticalEquivalog|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::HypotheticalEquivalog|=|common_name gene_symbol GO EC',
               
                ## characterized high confidence blast hit
                'UnirefBLASTP::Reviewed|=|common_name gene_symbol',
                
                ## pfam and non-equivalog tigrfams - full length
                'TIGRFAM::FullLength::Subfamily|=|common_name gene_symbol GO EC',
                'TIGRFAM::FullLength::Superfamily|=|common_name gene_symbol GO EC',
                'TIGRFAM::FullLength::EquivalogDomain|=|common_name gene_symbol GO EC',
                'TIGRFAM::FullLength::HypotheticalEquivalogDomain|=|common_name gene_symbol GO EC',
                'TIGRFAM::FullLength::SubfamilyDomain|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::Subfamily|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::Superfamily|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::EquivalogDomain|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::HypotheticalEquivalogDomain|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::SubfamilyDomain|=|common_name gene_symbol GO EC',
                
                ## pfam - fragment models, use all
                'PFAM::FRAG::Subfamily|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::Superfamily|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::EquivalogDomain|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::HypotheticalEquivalogDomain|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::SubfamilyDomain|=|common_name gene_symbol GO EC',
				
				## UniRef High Confidence
                'UnirefBLASTP::HighConfidence|=|common_name gene_symbol EC CAZY TAXON',

                'TIGRFAM::FullLength::Domain|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::Domain|=|common_name gene_symbol GO EC',
                'PFAM::FullLength::Uncharacterized|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::Domain|=|common_name gene_symbol GO EC',
                'PFAM::FRAG::Uncharacterized|=|common_name gene_symbol GO EC',

				## UniRef Medium Confidence
                'UnirefBLASTP::Putative|=|common_name gene_symbol TAXON',
                'UnirefBLASTP::ConservedDomain|=|common_name gene_symbol TAXON',

                'TMHMM|=|common_name',
                'LipoproteinMotif|=|common_name',
                'Hypothetical|=|common_name',
                
                ## UniRef Low Confidence
                'UnirefBLASTP::LowConfidence|=|TAXON',
                 );
   
    ## the handler hash defines which methods are responsible for doing the annotation operations
    ## current operations are '=' assign, '+' append, '-' overwrite, but other operations can
    ## be defined here and used in the annotation order array above
    my %handler = (
                    '=' =>  \&_assign,
                    '+' =>  \&_append,
                    '-' =>  \&_overwrite,
                  );


    
    ## annotate method annotates a polypeptide with associated AnnotationData objects
    ## passed as an argument
    sub annotate {
        my ($self, $pep,$outfh) = @_;
    
        my $pep_id = $pep->{'id'};
        
        ## will store our merged final annotation
        my $final_annotation = new CAMERA::AnnotationData::Polypeptide(type => 'final', source => {});
        $final_annotation->{'id'} = $pep_id;
            
        ## go through all of the annotation types in order
        foreach my $annot_rule(@annotation_order) {
                
            ## parse the annotation rule into type, action, attributes involved
            my ($annotation_type, $action, $attribs) = split(/\|/, $annot_rule);
            my @attributes = split(" ", $attribs);
   
            ## get annotation of the specified type
            my $annot_arr_ref = $pep->get_annotation($annotation_type);

            ## if we have any annotation of this type
            if ($annot_arr_ref && scalar @{$annot_arr_ref}) {
               
                ## sort annotation data of this type by rank
                my @sorted_annotation = sort _annotation_ref_sort @{$annot_arr_ref};
                my $annotation = $sorted_annotation[0];
                
                ## for each attribute involved, do the action
                foreach my $attribute(@attributes) {
                    $handler{$action}($final_annotation, $annotation, $attribute);
                } 

            }
        }
        
        ###### attribute-specific post-merging steps #############

        ## collapse child and parent ec numbers
        _do_ec_collapse($final_annotation);

        ## assign go terms based on final annotations
        
        #_do_tigr_role2go($final_annotation);
        #_do_ec2go($final_annotation);

        ## remove duplicate go ids
        _do_go_collapse($final_annotation);

        ## kgalens
        ## Added these function to implement some hacks
        ## These should be added to the parser functionality, but since I don't want
        ## rerun the parses, I put them here for now.  If this makes it into SVN
        ## feel free to delete it.
        _hack_this($final_annotation);
               
        #_hack_ec($final_annotation);

		my $com_name_arr_ref = $final_annotation->get_attribute('common_name');
		my @names = grep(/conserved hypothetical/, @{$com_name_arr_ref});
		if( @names ) {
			print Dumper( $final_annotation );
			exit(11);
		}

        ## set common_name to 'protein of unassigned function' if empty
        unless ($final_annotation->has_attribute('common_name')) {
            $final_annotation->add_attribute('common_name', 'protein of unassigned function');
        }

		## run the Protein Naming Utility on the final names
        _clean_common_names($self->{_synonymsDao},$final_annotation);

        ##########################################################
           
        ## output final annotation object to text table
        _to_table($final_annotation,$outfh);       
    }

    sub _hack_this {
        my ($final_annotation) = @_;
		my $names_arr_ref = $final_annotation->get_attribute('common_name');
       	
       	foreach my $common_name( @{$names_arr_ref} ) {
        
	        unless( defined( $common_name ) ) {
	            _hack_make_hypo( $final_annotation );
	            return;
	        }
	
	        #If the name is less than three characters
	        if( length($common_name) < 3 ) {
	            _hack_make_hypo( $final_annotation );
				last;
	        }
	
	        #If the name starts with a dash
	        if( $common_name =~ /^[^\w\d]/ ) {
	            _hack_make_hypo( $final_annotation );
				last;
	        }
	
	        #If the common name is more than 75% numbers
	        my $no_space = $common_name;
	        $no_space =~ s/\s+//;
	        
	        my $full_length = length($no_space);
	        $no_space =~ s/\D//g;
	        
	        my $num_length = length($no_space);
	        my $num_percent = int(($num_length/$full_length)*100);
	
	        if($num_percent > 60 ) {
	            _hack_make_hypo( $final_annotation );
				last;
        	}

	        #Reverse hypo.  It happens...
	        if( $common_name =~ /lacitehtopyh/ ) {
	            _hack_make_hypo( $final_annotation );
			last;
	        } 

			#Conserved Hypothetical
			if( $common_name =~ /conserved hypothetical/ ) {
				_hack_make_hypo( $final_annotation );
				last;
			}  
		
			#Common name is taxon:(\d+)
			if( $common_name =~ /taxon\:\d+/ ) {
				_hack_make_hypo( $final_annotation );
				last;
			}
		}
    }

    sub _hack_ec {
        my ($final_annotation) = @_;
        return unless( $final_annotation->has_attribute('EC') );
        
        my $valid_ec_file = "/usr/local/db/calit_db/EC/enzyme.tab";
        open( my $ec_fh, "< $valid_ec_file") or die("Could not open $valid_ec_file ($!)");

        my $ec_num = $final_annotation->get_attribute('EC');
        $ec_num =~ s/./\./g;
        $ec_num =~ s/-/./g;

        my @lines = grep(/$ec_num/,<$ec_fh>);
        close( $ec_fh );
        
        unless( @lines ) {
            _hack_make_hypo($final_annotation);
        }
    }
    
    sub _clean_common_names() {
	    my ($pnu,$final_annotation) = @_;
	    my @cleanNames = ();
	    
       	my $names_arr_ref = $final_annotation->get_attribute('common_name');
        
       	foreach my $common_name( @{$names_arr_ref} ) {
       		## lower case any words that start with a capital letter followed
       		## by a lower case letter
       		$common_name =~ s/\b([A-Z][a-z])(\w+)\b/\l$1$2/g;
       		
			push(@cleanNames,$pnu->rename($common_name));        
        }
        
        $final_annotation->replace_attribute('common_name', \@cleanNames);
    }
    
    ## Another hack subroutine
    sub _hack_make_hypo {
        my ($final_annotation) = @_;

        my @att_list = @{CAMERA::AnnotationData::ProteinDataTypes::attribute_types};
        foreach my $attribute ( @att_list ) {
            if( $final_annotation->has_attribute( $attribute ) ) {
                $final_annotation->clear_attribute( $attribute );
            }
        }

        $final_annotation->replace_attribute( 'common_name', ['hypothetical protein'] );
        
    }

    ## sort annotation object array ref by rank
    sub _annotation_ref_sort {
        $$a{'rank'} <=> $$b{'rank'};
    }

    ## assign method does the assign action of the annotation rules
    ## will assign attributes of a given type of they are currently unassigned
    sub _assign {
        my ($final_annotation, $annotation, $attribute) = @_;

        ## if hasn't been assigned yet in final_annotation then assign
        unless ($final_annotation->has_attribute($attribute)) {
            my $attrib_arr_ref = $annotation->get_attribute($attribute);
            if ($attrib_arr_ref) {
                foreach my $att(@{$attrib_arr_ref}) {
                    $final_annotation->add_attribute($attribute, $att);
                    push (@{$final_annotation->{'source'}->{$attribute}}, $annotation->{'source'});
                }
            }
        }
    }
    
    ## append method does the append action of the annotation rules
    ## adds annotation attribute to the final annotation even if it currently has an assigned value
    sub _append {
        my ($final_annotation, $annotation, $attribute) = @_;
        
        ## if has been assigned in final_annotation then append
        my $attrib_arr_ref = $annotation->get_attribute($attribute);
        if ($attrib_arr_ref) {
            foreach my $att(@{$attrib_arr_ref}) {
                $final_annotation->add_attribute($attribute, $att);
                push (@{$final_annotation->{'source'}->{$attribute}}, $annotation->{'source'});
            }
        }
    }
    
    ## overwrite method does the overwrite action of the annotation rules
    ## overwrites any currently assigned attribute value
    sub _overwrite { 
        my ($final_annotation, $annotation, $attribute) = @_;

        ## if has been assigned in final_annotation then overwrite
        my $attrib_arr_ref = $annotation->get_attribute($attribute);
        if ($attrib_arr_ref) {
            $final_annotation->clear_attribute($attribute);
            foreach my $att(@{$attrib_arr_ref}) {
                $final_annotation->add_attribute($attribute, $att);
                push (@{$final_annotation->{'source'}->{$attribute}}, $annotation->{'source'});
            }
        }
 

    }

    ## collapse down ec numbers to the set of most specific levels available
    sub _do_ec_collapse {
        my ($annotation) = @_;

        ## WARNING: this needs to be fixed for maintaining source associations/counts
        ##
        ## source should correspond to:
        ##
        ##  1) the most specific EC number
        ##  2) the higher ranked annotation (first in the array)

        my $ec_arr_ref = $annotation->get_attribute('EC');
        my $ec_source_arr_ref = $annotation->{'source'}->{'EC'};
        
        if ($ec_arr_ref && scalar(@{$ec_arr_ref}) > 1) {
        
            my $collapsed_ec_arr_ref = _ec_collapse($ec_arr_ref, $ec_source_arr_ref);

            $annotation->replace_attribute('EC', $collapsed_ec_arr_ref);
        }
    }
    
    sub _ec_collapse {
        my ($ec_arr_ref, $ec_source_arr_ref) = @_;

        return CAMERA::AnnotationRules::Util::ec_collapse($ec_arr_ref, $ec_source_arr_ref);
    }
    
    ## collapse down set of assigned go ids to remove duplicates
    sub _do_go_collapse {
        my ($annotation) = @_;
        
        ## WARNING: this needs to be fixed for maintaining source associations/counts
        ##
        ## source should correspond to:
        ##
        ##  1) the most specific EC number
        ##  2) the higher ranked annotation (first in the array)

        my $go_arr_ref = $annotation->get_attribute('GO');
        my $go_source_arr_ref = $annotation->{'source'}->{'GO'};

        if ($go_arr_ref && scalar(@{$go_arr_ref}) > 1) {
        
            #my $collapsed_go_arr_ref = _go_collapse($go_arr_ref, $go_source_arr_ref);
	    my $collapsed_go_arr_ref = [];
	    my %tmp;
	    foreach my $go ( @{$go_arr_ref} ) {
		$tmp{$go} = 1;	
	    }
            @{$collapsed_go_arr_ref} = keys( %tmp );
            $annotation->replace_attribute('GO', $collapsed_go_arr_ref);
        }
    }
    
    sub _go_collapse {
        my ($go_arr_ref, $go_source_arr_ref) = @_;

        return CAMERA::AnnotationRules::Util::go_collapse($go_arr_ref, $go_source_arr_ref);
    }

    
    ## converts final annotation objects into table rows
    #
    ## this should be moved to a different class maybe? though source is specific to the usage in this class
    sub _to_table {
        my ($annotation,$out) = @_;

        my $row = $annotation->{'id'};
        
        foreach my $key(@CAMERA::AnnotationData::PolypeptideDataTypes::attribute_list) {
		
            my $att = ($annotation->{'atts'}->{$key}) ? join(" || ", @{$annotation->{'atts'}->{$key}}) : '';
            my $source = ($annotation->{'source'}->{$key}) ? join(" || ", @{$annotation->{'source'}->{$key}}) : '';
            
            $row .= "\t".$key."\t".$att."\t".$source;
        }
	

	## took out prining to the output file because somehow
	## the filehandle was getting closed and throwing errors.
	$row .= "\n";	
        print $row;
    }
}
1;
