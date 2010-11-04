#!/usr/local/bin/perl

package CAMERA::AnnotationData;

use strict;
use warnings;
use Carp;

my $attribute_types = {};

{

    sub new {
        my ($class, %args) = @_;
        
        my $id = ($args{'id'}) ? $args{'id'} : undef;
        my $type = ($args{'type'}) ? $args{'type'} : undef;
        my $source = ($args{'source'}) ? $args{'source'} : undef;
        my $rank = ($args{'rank'}) ? $args{'rank'} : 1;
       
        my $annotation = {  'id'           => $id,
                            'type'         => $type,
                            'source'       => $source,
                            'rank'         => $rank,
                         };
       
        return bless $annotation, ref($class) || $class;
    }

    sub has_attribute {
        my ($self, $type) = @_;

        if ($self->{'atts'}->{$type}) {
            return 1;
        } else {
            return 0;
        }
    }
    
    sub get_attribute {
        my ($self, $type) = @_;

        if ($self->{'atts'}->{$type}) {
            return $self->{'atts'}->{$type};
        } else {
            return undef;
        }
    }

    sub add_attribute {
        my ($self, $key, $value) = @_;

        push (@{$self->{'atts'}->{$key}}, $value);

        return 1;
    }
    
    sub replace_attribute {
        my ($self, $key, $attrib_arr_ref) = @_;
    
        unless ($attrib_arr_ref) {
            confess "replace_attribute called with invalid arguments";
        }
        
        $self->{'atts'}->{$key} = $attrib_arr_ref;

        return 1;
    }
    
    sub clear_attribute {
        my ($self, $key) = @_;

        $self->{'atts'}->{$key} = [];

        return 1;
    }

    sub unset_attribute {
        my ($self, $key) = @_;

        $self->clear_attribute($key);
        
        return 1;
    }

    sub get_attribute_keys {
        my ($self) = @_;

        return keys(%{$self->{'atts'}});
    }

}

1;
