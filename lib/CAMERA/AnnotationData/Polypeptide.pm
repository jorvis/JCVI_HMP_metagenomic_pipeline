#!/usr/local/bin/perl

package CAMERA::AnnotationData::Polypeptide;
@ISA = qw( CAMERA::AnnotationData );

use strict;
use warnings;
use CAMERA::AnnotationData;
use CAMERA::AnnotationData::PolypeptideDataTypes;
use Carp;

my $attribute_types = {};
foreach my $type(@CAMERA::AnnotationData::PolypeptideDataTypes::attribute_list) {
    $attribute_types->{$type} = 1;
}

sub add_attribute {
    my ($self, $key, $value) = @_;

    unless ($attribute_types->{$key}) {
        confess "Attempted to add an unrecognized attribute type '$key'";
    }
    
    if ($value =~ / \|\| /) {
        confess "Can't attribute '$key' value contains record seperator string ' || ':\n$value";
    }

    $self->SUPER::add_attribute($key, $value);
}

sub to_string {
    my ($self) = @_;

    my $id = $self->{'id'};
    my $source = $self->{'source'};
    my $type = $self->{'type'};
    my $rank = $self->{'rank'};

    my $row = "$id\t$type\t$source\t$rank";
    
    foreach my $key(@CAMERA::AnnotationData::PolypeptideDataTypes::attribute_list) {
        my $att = ($self->{'atts'}->{$key}) ? join(" || ", @{$self->{'atts'}->{$key}}) : '';
        $row .= "\t".$key."\t".$att;
    }
    $row .= "\n";
    
    return $row;

}

sub from_string {
    my ($self, $string) = @_;

    chomp($string);

    my ($id, $type, $source, $rank, $attributes) = split(/\t/, $string, 5);

    my @attribs = split(/\t/, $attributes);
    if (@attribs % 2) {
        push(@attribs, '');
    }
    my %attribs = @attribs;
 
    $self->{'id'} = $id;
    $self->{'type'} = $type;
    $self->{'source'} = $source;
    $self->{'rank'} = $rank;
    foreach my $key(keys(%attribs)) {
        my @vals = split(/ \|\| /, $attribs{$key});
        foreach my $val(@vals) {
            $self->add_attribute($key, $val);
        }
    }
}

1;
