package BSML::BsmlSequenceAlignment;
@ISA = qw( BSML::BsmlElement );

BEGIN {
use BSML::BsmlSequenceData;
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
    $self->{'BsmlSequenceData'} = [];
    $self->{'Alignment-consensus'} = '';
}

sub addBsmlAlignmentConsensus
{
    my $self = shift;
    my $consensus = shift;

    $self->{'Alignment-consensus'} = $consensus;
}

sub returnBsmlAlignmentConsensus
{
    my $self = shift;
    return $self->{'Alignment-consensus'};
}

sub addBsmlSequenceData
{
    my $self = shift;
    push( @{$self->{'BsmlSequenceData'}}, new BSML::BsmlSequenceData );

    my $index = @{$self->{'BsmlSequenceData'}} - 1;
    return $index;
}

sub dropBsmlSequenceData
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlSequenceData'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlSequenceData'}[$i] );
	  }
      }

    $self->{'BsmlSequenceData'} = \@newlist;    
  }

sub returnBsmlSequenceDataListR
{
    my $self = shift;
    return $self->{'BsmlSequenceData'};
}

sub returnBsmlSequenceDataR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlSequenceData'}[$index];
}

sub write
{
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Sequence-alignment", %{$self->{'attr'}} );
     
    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $seqDat ( @{$self->{'BsmlSequenceData'}} )
    {
	$seqDat->write( $writer );
    }

    $writer->startTag( "Alignment-consensus" );
    $writer->characters( $self->{'Alignment-consensus'} );
    $writer->endTag( "Alignment-consensus" );

    

    foreach my $link (@{$self->{'BsmlLink'}})
    {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
    }

    $writer->endTag( "Sequence-alignment" );
}
    
1
