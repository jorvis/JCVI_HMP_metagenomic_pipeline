package BSML::BsmlSequenceData;
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
    $self->{'seqAlignmentDat'} = '';
    $self->{'BsmlLink'} = [];
}

# Set the alignment sequence data
sub addSequenceAlignmentData
{
    my $self = shift;

    my ($seqDat) = @_;
    $self->{'seqAlignmentDat'} = $seqDat;
}

# return the sequence alignment

sub returnSequenceAlignmentData
{
    my $self = shift;
    return $self->{'seqAlignmentDat'};
}

sub write
{
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Sequence-data", %{$self->{'attr'}} );
     
    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }


    $writer->characters( $self->{'seqAlignmentDat'} );

    foreach my $link (@{$self->{'BsmlLink'}})
    {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
    }

    $writer->endTag( "Sequence-data" );
}
    
1
