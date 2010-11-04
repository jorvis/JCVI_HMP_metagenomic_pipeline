package BSML::BsmlSegmentSet;
@ISA = qw( BSML::BsmlElement );

use XML::Writer;
use strict;
use warnings;

BEGIN {
use BSML::BsmlSegment;
use BSML::BsmlElement;
}

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
    $self->{'BsmlSegments'} = [];
}

sub addBsmlSegment
{
    my $self = shift;
    
    push( @{$self->{'BsmlSegments'}}, new BSML::BsmlSegment );
    
    my $index = @{$self->{'BsmlSegment'}} - 1;

    return $index;    
}

sub dropBsmlStrain
{
    my $self = shift;
    
    my ($index) = @_;
    
    my $newlist;
    
    for(  my $i=0;  $i< @{$self->{'BsmlSegments'}}; $i++ ) 
    {
	if( $i != $index )
	{
	    push( @{$newlist}, $self->{'BsmlSegments'}[$i] );
	}
    }
    
    $self->{'BsmlSegments'} = $newlist;
}

sub returnBsmlStrainListR
{
    my $self = shift;
    return $self->{'BsmlSegments'};
}

sub returnBsmlStrainR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlSegments'}[$index];
}

sub write
{
    my $self = shift;
    my $writer = shift;
    
    $writer->startTag( "Segment-set", %{$self->{'attr'}} );

    foreach my $segment (@{$self->{'BsmlSegments'}})
    {
	$segment->write( $writer );
    }

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $link (@{$self->{'BsmlLink'}})
    {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
    }

    $writer->endTag( "Segment-set" );
}

1
