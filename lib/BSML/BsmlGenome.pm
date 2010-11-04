package BSML::BsmlGenome;
@ISA = qw( BSML::BsmlElement );

BEGIN {
use BSML::BsmlOrganism;
use BSML::BsmlCrossReference;
use BSML::BsmlElement;
}
use XML::Writer;
use Data::Dumper;


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
    $self->{'BsmlOrganism'} = undef;
    $self->{'BsmlChromosomes'} = [];
    $self->{'BsmlCrossReference'} = [];
}

sub addBsmlOrganism
{
    my $self = shift;

    $self->{'BsmlOrganism'} = new BSML::BsmlOrganism;
    return $self->{'BsmlOrganism'};
}

sub returnBsmlOrganismR
{
    my $self = shift;

    return $self->{'BsmlOrganism'};
}



sub addBsmlChromosome
{
    my $self = shift;
}

sub dropBsmlChromosome
{
    my $self = shift;

}

sub returnBsmlChromosomeListR
{
    my $self = shift;
    return $self->{'BsmlChromosomes'};
}

sub returnBsmlChromosomeR
{
    my $self = shift;
    my $index = shift;

    return $self->{'BsmlChromosomes'}->[$index];
}

sub write
{
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Genome", %{$self->{'attr'}} );

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    if( my $org = $self->{'BsmlOrganism'} )
    {
	$org->write( $writer );
    }

#    print Dumper $self;

    foreach my $xref (@{$self->{'BsmlCrossReference'}}){


	$xref->write( $writer );
    }


    foreach my $link (@{$self->{'BsmlLink'}})
    {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
    }

    $writer->endTag( "Genome" );
}
