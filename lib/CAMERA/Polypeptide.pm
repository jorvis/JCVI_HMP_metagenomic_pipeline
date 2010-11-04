#!/usr/local/bin/perl

package CAMERA::Polypeptide;

use strict;
use warnings;
use Carp;
use CAMERA::AnnotationData::PolypeptideDataTypes;

{
    
    my @pep_data_types = @CAMERA::AnnotationData::PolypeptideDataTypes::type_list;
    unless(@pep_data_types) {
        confess "Polypeptide annotation data types array is empty";
    }
   
    my %valid_type = ();
    foreach my $type(@pep_data_types) {
        $valid_type{$type} = 1;
    }
    
    sub new {
        my ($class, %args) = @_;
        my $self = {};

        $self->{'id'} = $args{'id'};        
        $self->{'annotation'} = {};
        
        return bless $self, $class; 
    }
    
    ## return true if the protein accession indicates its a public database protein
    sub is_public_database_protein {
        my ($self) = @_;

        if ($self->{'id'} =~ /^NCBI_|^SPROT_/) {
            return 1;
        } else {
            return 0;
        }
    }

    ## add an annotation object to the protein
    sub add_annotation {
        my ($self, $annotation_obj_ref) = @_;

        #print join("\t", caller())."\t".$annotation_obj_ref->{'source'}."\t".$annotation_obj_ref->{'type'}."\n";
        
        unless ($valid_type{$annotation_obj_ref->{'type'}}) {
            confess "Attempted to add invalid annotation object type '$annotation_obj_ref->{type}'";
        }
        
        push(@{$self->{'annotation'}->{$annotation_obj_ref->{'type'}}}, $annotation_obj_ref);
        
        return 1;
    }

    ## check if the protein has any annotation of a specific type
    sub has_annotation {
        my ($self, $type) = @_;

        if ($self->{'annotation'}->{$type}) {
            return 1;
        } else {
            return 0;
        }
    }
    
    ## return the array reference for annotation of specific type
    sub get_annotation {
        my ($self, $type) = @_;
        
        return $self->{'annotation'}->{$type};
    }
    
    ## return the array reference for annotation of specific type
    sub get_annotation_keys {
        my ($self, $type) = @_;
        
        return keys(%{$self->{'annotation'}});
    }
}

1;
