package BSML::BsmlFeatureGroup;
@ISA = qw( BSML::BsmlElement );

=head1 NAME

BsmlFeatureGroup.pm - Bsml Feature Group Class

=head1 VERSION

This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS

  Adding a feature group member to a feature group

  addBsmlFeatureGroupMember( 'Cald1', 'GENE', '', '' );

=head1 DESCRIPTION

=head2 Overview

  The BsmlFeatureGroup class provides method and storage for Bsml Feature-group and Feature-group-member elements

=head2 Constructor and initialization

  Typically  BsmlFeatureGroup is created by the BsmlSequence object it is contained within and manipulated as a 
  reference.

  my $FGroup = $seq->returnBsmlFeatureGroupR( $seq->addBsmlFeatureGroup() );

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

# A BsmlFeatureGroup stores a list of feature-group-members. Feature-group-members
# are represented as Perl hashes containing the keys for: feature-id, feature-type, 
# group-type, and the CDATA text for the element 

sub init
  {
    my $self = shift;

    $self->{ 'attr' } = {};
    $self->{ 'BsmlAttr' } = {};
    $self->{ 'BsmlFeatureGroupMembers' } = [];
    $self->{ 'text' } = '';
    $self->{ 'BsmlLink' } = [];

    #ParentSequenceId is set in BsmlSequence::addBsmlFeatureGroup()...
    #In order to retain the relationship of genes to the assembly on which they
    #are contained in a memory efficient manner consistent with the document level
    #lookups, the sequence id is embedded in each feature group. 

    $self->{ 'ParentSequenceId' } = ''; 
  }

=item $FGroup->addBsmlFeatureGroupMember( $feature_id, $feature_type, $group_type, $cdata )

B<Description:> adds a feature group memeber to the feature group

B<Parameters:> ( $feature_id, $feature_type, $group_type, $cdata ) - 
  $feature_id - the id of the feature to be added to the feature group (unique at the document level)
  $feature_type - the type of feature
  $group_type - the group type the feature belongs to
  $cdata - the character data associated with the feature

B<Returns:> The index of the added feature group member

=cut

sub addBsmlFeatureGroupMember
  {
    my $self = shift;
    my ( $feat, $feature_type, $group_type, $text ) = @_;

    my $href = { 'feature' => $feat, 'feature-type' => $feature_type, 'group-type' => $group_type, 'text' => $text };

    push( @{$self->{'BsmlFeatureGroupMembers'}}, $href );

    my $index = @{$self->{'BsmlFeatureGroupMembers'}} - 1;
    return $index;
  }

=item $FGroup->dropBsmlFeatureGroupMember( $index )

B<Description:> deletes a feature group member from the feature group

B<Parameters:> $index - the index of the feature group member to be deleted

B<Returns:> none

=cut

sub dropBsmlFeatureGroupMember
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlFeatureGroupMembers'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlFeatureGroupMembers'}[$i] );
	  }
      }

    $self->{'BsmlFeatureGroupMembers'} = \@newlist;    
  }

=item $FGroup->setText( $text )

B<Description:> sets the CDATA associated with the feature group

B<Parameters:> The text to be set

B<Returns:> none

=cut

sub setText
  {
    my $self = shift;
    my ($text) = @_;

    $self->{'text'} = $text;
  }

=item $FGroup->returnText()

B<Description:> returns the CDATA associated with the feature group

B<Parameters:> None

B<Returns:> the CDATA associated with the feature group

=cut

sub returnText
  {
    my $self = shift;
    return $self->{'text'};
  }

=item $FGroup->returnFeatureGroupMemberListR()

B<Description:> returns a list of feature-group-members (hash references)

B<Parameters:> none

B<Returns:> returns a list of feature-group-members (hash references)

=cut

sub returnFeatureGroupMemberListR
  {
    my $self = shift;
    return $self->{'BsmlFeatureGroupMembers'};
  }

=item $FGroup->returnFeatureGroupMemberR( $index )

B<Description:> returns a feature-group-member (hash reference)

B<Parameters:> $index - the index of the feature-group-member

B<Returns:> returns a feature-group-member (hash references)

=cut

sub returnFeatureGroupMemberR
  {
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlFeatureGroupMembers'}[$index];
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

    $writer->startTag( 'Feature-group', %{$self->{'attr'}} );

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $bsmlfeaturemember ( @{$self->{'BsmlFeatureGroupMembers'}} )
      {
	my $h = {};

	if( $bsmlfeaturemember->{'feature'} ){ $h->{'featref'} = $bsmlfeaturemember->{'feature'};}
	if( $bsmlfeaturemember->{'feature-type'} ){$h->{'feature-type'} = $bsmlfeaturemember->{'feature-type'};}
	if( $bsmlfeaturemember->{'group-type'} ){$h->{'group-type'} = $bsmlfeaturemember->{'group-type'};}

	$writer->startTag( "Feature-group-member", %{ $h } );
	if( $bsmlfeaturemember->{'text'} ){ $writer->characters( $bsmlfeaturemember->{'text'} ); }
	$writer->endTag( "Feature-group-member" );
      }

    foreach my $link (@{$self->{'BsmlLink'}})
      {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
      }
    
    $writer->endTag( "Feature-group" );

  }

sub returnParentSequenceId
  {
    my $self = shift;
    return $self->{'ParentSequenceId'};
  }

1

