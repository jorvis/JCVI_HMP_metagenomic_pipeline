package BSML::BsmlCrossReference;
@ISA = qw( BSML::BsmlElement );

=head1 NAME

BsmlCrossReference.pm - Bsml Cross Reference Class

=head1 VERSION

This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS



=head1 DESCRIPTION

=head2 Overview

  The BsmlCrossReference class provides method and storage for Bsml Cross-reference elements

=head2 Constructor and initialization





=over 4

=cut

BEGIN {
use BSML::BsmlElement;
}
use XML::Writer;
use strict;
use warnings;

sub new
  {
    my $class = shift;
    my $self = {};
    bless $self, $class;
    
    $self->init();
    return $self;
  }


sub init
  {
    my $self = shift;

    $self->{ 'attr' } = {};
    $self->{ 'BsmlAttr' } = {};
    $self->{ 'BsmlLink' } = [];

  }

sub write {
    my $self = shift;
    my $writer = shift;
    
    $writer->startTag( 'Cross-reference', %{$self->{'attr'}} );
    
    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $link (@{$self->{'BsmlLink'}}){
	$writer->startTag( "Link", %{$link} );
	$writer->endTag( "Link" );
    }

    $writer->endTag( "Cross-reference" );
    
}

1

