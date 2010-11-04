package BSML::BsmlStrain;
@ISA = qw( BSML::BsmlElement );

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
    
    $self->{'attr'} = {};
    $self->{'BsmlAttr'} = {};
    $self->{'BsmlLink'} = [];
}

sub write
{
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Strain", %{$self->{'attr'}} );

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

    $writer->endTag( "Strain" );
}

1
