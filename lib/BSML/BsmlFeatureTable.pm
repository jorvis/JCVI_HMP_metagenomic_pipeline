package BSML::BsmlFeatureTable;
@ISA = qw( BSML::BsmlElement );

=head1 NAME

BsmlFeatureTable.pm - Bsml Feature Table Class

=head1 VERSION

This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS

  Adding a feature and reference to a feature table

  my $feature = $table->returnBsmlFeatureR( $table->addBsmlFeature() );
  my $ref = $table->returnBsmlReferenceR( $table->addBsmlReference() );

=head1 DESCRIPTION

=head2 Overview

  The BsmlFeatureTable class provides methods and storage for Bsml Feature and Resource elements.

=head2 Constructor and initialization

  Typically a BsmlFeatureTable is created by the BsmlSequence object it is contained within and
  manipulated as a reference.

  my $FTable = $seq->returnFeatureTableR( addBsmlFeatureTable() );

=head2 Class and object methods

=over 4

=cut

BEGIN {
use BSML::BsmlReference;
use BSML::BsmlFeature;
use BSML::BsmlElement;
}
use XML::Writer;

use strict;
use warnings;

# A feature table contains a list of references and a list of features

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
    $self->{ 'BsmlFeatures' } = [];
    $self->{ 'BsmlReferences' } = [];
    $self->{ 'BsmlLink' } = [];
  }

=item $FTable->addBsmlFeature()

B<DESCRIPTION:> adds a feature to the feature table

B<Parameters:> None

B<Returns:> index of the created feature (position of the feature in the feature list)

=cut

sub addBsmlFeature
  {
    my $self = shift;
    push( @{$self->{'BsmlFeatures'}}, new BSML::BsmlFeature );

    my $index = @{$self->{'BsmlFeatures'}} - 1;
    return $index;
  }

=item $FTable->dropBsmlFeature( $index )

B<DESCRIPTION:> removes a feature in the feature table with given index

B<Parameters:> $index, position of the feature in the feature list

B<Returns:> None

=cut

sub dropBsmlFeature
  {
    my $self = shift;
    my ($index) = @_;

    my $newlist = [];

    for(  my $i=0;  $i< @{$self->{'BsmlFeatures'}}; $i++ ) 
      {
	if( $i != $index )
	  {
	    push( @{$newlist}, $self->{'BsmlFeatures'}[$i] );
	  }
      }

    $self->{'BsmlFeatures'} = $newlist;
  }

=item $FTable->returnBsmlFeatureListR()

B<DESCRIPTION:> returns a reference to the list of feature objects contained in the feature table

B<Parameters:> None

B<Returns:> list reference

=cut

sub returnBsmlFeatureListR
  {
    my $self = shift;

    return $self->{'BsmlFeatures'};
  }

=item $FTable->returnBsmlFeatureR( $index )

B<DESCRIPTION:> returns a reference to the BsmlFeature at position $index in the feature list

B<Parameters:> $index - the features position in the feature list

B<Returns:> reference to a BsmlFeature

=cut

sub returnBsmlFeatureR
  {
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlFeatures'}[$index];
  }

=item $FTable->addBsmlReference()

B<DESCRIPTION:> adds a reference to the feature table

B<Parameters:> None

B<Returns:> index of the created reference (position of the reference in the reference list)

=cut

sub addBsmlReference
  {
    my $self = shift;
    push( @{$self->{'BsmlReferences'}}, new BSML::BsmlReference );

    my $index = @{$self->{'BsmlReferences'}} - 1;
    return $index;
  }

=item $FTable->dropBsmlReference( $index )

B<DESCRIPTION:> removes a reference in the reference table with given index

B<Parameters:> $index, position of the reference in the reference list

B<Returns:> None

=cut

sub dropBsmlReference
  {
    my $self = shift;
    my ($index) = @_;

    my $newlist = [];

    for(  my $i=0;  $i < @{$self->{'BsmlReferences'}}; $i++ ) 
      {
	if( $i != $index )
	  {
	    push( @{$newlist}, $self->{'BsmlReferences'}[$i] );
	  }
      }

    $self->{'BsmlReferences'} = $newlist;
  }

=item $FTable->returnBsmlReferenceListR()

B<DESCRIPTION:> returns a reference to the list of BsmlReference objects contained in the feature table

B<Parameters:> None

B<Returns:> reference to a list of BsmlReferences

=cut

sub returnBsmlReferenceListR
  {
    my $self = shift;
    
    return $self->{'BsmlReferences'};
  }

=item $FTable->returnBsmlReferenceR( $index )

B<DESCRIPTION:> returns a reference to the BsmlReference at position $index in the reference list

B<Parameters:> $index - the references position in the reference list

B<Returns:> reference to a BsmlReference

=cut

sub returnBsmlReferenceR
  {
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlReferences'}[$index];
  }

=item $seq->write()

  B<Description:> writes the BSML elements encoded by the class to a file using XML::Writer. This method should only be called through the BsmlDoc->write() process.

  B<Parameters:> None

  B<Returns:> None

=cut 

sub write
  {
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Feature-table", %{$self->{'attr'}} );

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $bsmlfeature (@{$self->{'BsmlFeatures'}})
      {
	$bsmlfeature->write( $writer );
      }

    foreach my $bsmlref (@{$self->{'BsmlReferences'}})
      {
	$bsmlref->write( $writer );
      }

    foreach my $link (@{$self->{'BsmlLink'}})
      {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
      }

    $writer->endTag( "Feature-table" );
  }

1
