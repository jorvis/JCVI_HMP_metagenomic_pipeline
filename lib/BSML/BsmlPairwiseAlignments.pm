package BSML::BsmlPairwiseAlignments;
@ISA = qw( BSML::BsmlElement );

BEGIN {
use BSML::BsmlElement;
use BSML::BsmlAlignedPair;
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
    $self->{'BsmlAlignedPairs'} = [];
}

sub addBsmlAlignedPair
{
    my $self = shift;
    push( @{$self->{'BsmlAlignedPairs'}}, new BSML::BsmlAlignedPair );

    my $index = @{$self->{'BsmlAlignedPairs'}} - 1;
    return $index;
}

sub dropBsmlAlignedPair
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlAlignedPairs'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlAlignedPairs'}[$i] );
	  }
      }

    $self->{'BsmlAlignedPairs'} = \@newlist;    
  }

sub returnBsmlAlignedPairListR
{
    my $self = shift;
    return $self->{'BsmlAlignedPairs'};
}

sub returnBsmlAlignedPairR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlAlignedPairs'}[$index];
}

sub write
{
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Pairwise-alignments", %{$self->{'attr'}} );
     
    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $alignedPair ( @{$self->{'BsmlAlignedPairs'}} )
    {
	$alignedPair->write( $writer );
    }

    foreach my $link (@{$self->{'BsmlLink'}})
    {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
    }

    $writer->endTag( "Pairwise-alignments" );
}
    
1
