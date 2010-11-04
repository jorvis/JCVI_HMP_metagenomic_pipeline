package BSML::BsmlAlignmentSummary;
@ISA = qw( BSML::BsmlElement );

BEGIN {
use BSML::BsmlAlignedSequence;
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
    
    $self->{'attr'} = {};
    $self->{'BsmlAttr'} = {};
    $self->{'BsmlLink'} = [];
    $self->{'BsmlAlignedSequences'} = [];
}

sub addBsmlAlignedSequence
{
    my $self = shift;
    push( @{$self->{'BsmlAlignedSequences'}}, new BSML::BsmlAlignedSequence );

    my $index = @{$self->{'BsmlAlignedSequences'}} - 1;
    return $index;
}

sub dropBsmlAlignedSequence
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlAlignedSequences'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlAlignedSequences'}[$i] );
	  }
      }

    $self->{'BsmlAlignedSequences'} = \@newlist;    
  }

sub returnBsmlAlignedSequenceListR
{
    my $self = shift;
    return $self->{'BsmlAlignedSequences'};
}

sub returnBsmlAlignedSequenceR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlAlignedSequences'}[$index];
}

sub write
{
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Alignment-summary", %{$self->{'attr'}} );
     
    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $alignedSequence ( @{$self->{'BsmlAlignedSequences'}} )
    {
	$alignedSequence->write( $writer );
    }

    foreach my $link (@{$self->{'BsmlLink'}})
    {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
    }

    $writer->endTag( "Alignment-summary" );
}
    
1
