package BSML::BsmlReference;
@ISA = qw( BSML::BsmlElement );

=head1 NAME

  BsmlReference.pm - Bsml API Object representing the Bsml Reference Element

=head1 VERSION

  This document refers to version 1.0 of the BSML Object Layer

=head1 DESCRIPTION

=head2 Overview

  The BsmlReference class allows for storage and manipulation of Bsml Reference elements 

=head2 Constructor and initialization

=head2 Class and object methods

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
    $self->{ 'BsmlRefAuthors' } = '';
    $self->{ 'BsmlRefTitle' } = '';
    $self->{ 'BsmlRefJournal' } = '';
    $self->{ 'BsmlLink' } = [];
  }

sub addBsmlRefAuthors
  {
    my $self = shift;
    my ( $authors ) = @_;

    $self->{ 'BsmlRefAuthors' } = $authors;
  }

sub setBsmlRefAuthors
  {
    my $self = shift;
    my ( $authors ) = @_;
    
    $self->addBsmlRefAuthors( $authors );
  }

sub dropBsmlRefAuthors
  {
    my $self = shift;
    
    $self->{ 'BsmlRefAuthors' } = '';
  }

sub returnBsmlRefAuthors
  {
    my $self = shift;
    return $self->{ 'BsmlRefAuthors' };
  }

sub addBsmlRefTitle
  {
    my $self = shift;
    my ( $title ) = @_;

    $self->{ 'BsmlRefTitle' } = $title;
  }

sub setBsmlRefTitle
  {
    my $self = shift;
    my ( $title ) = @_;

    $self->addBsmlRefTitle( $title );
  }

sub dropBsmlRefTitle
  {
    my $self = shift;
    $self->{ 'BsmlRefTitle' } = '';
  }

sub returnBsmlRefTitle
  {
    my $self = shift;
    return $self->{'BsmlRefTitle'};
  }


sub addBsmlRefJournal
  {
    my $self = shift;
    my ( $journal ) = @_;

    $self->{ 'BsmlRefJournal' } = $journal;
  }

sub setBsmlRefJournal
  {
    my $self = shift;
    my ( $journal ) = @_;

    $self->addBsmlRefJournal( $journal );
  }

sub dropBsmlRefJournal
  {
    my $self = shift;
    $self->{'BsmlRefJournal'} = '';
  }

sub returnBsmlRefJournal
  {
    my $self = shift;
    return $self->{'BsmlRefJournal'};
  }

sub write
  {
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Reference", %{$self->{'attr'}} );

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    if( $self->{'BsmlRefAuthors'} )
      {
	$writer->startTag( "RefAuthors" );
	$writer->characters( $self->{'BsmlRefAuthors'} );
	$writer->endTag( "RefAuthors" );
      }

    if( $self->{'BsmlRefTitle'} )
      {
	$writer->startTag( "RefTitle" );
	$writer->characters( $self->{'BsmlRefTitle'} );
	$writer->endTag( "RefTitle" );
      }

     if( $self->{'BsmlRefJournal'} )
      {
	$writer->startTag( "RefJournal" );
	$writer->characters( $self->{'BsmlRefJournal'} );
	$writer->endTag( "RefJournal" );
      }

    foreach my $link (@{$self->{'BsmlLink'}})
      {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
      }
    
    $writer->endTag( "Reference" );
  }

1
