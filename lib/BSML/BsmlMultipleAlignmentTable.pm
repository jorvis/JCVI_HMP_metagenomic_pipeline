package BSML::BsmlMultipleAlignmentTable;
@ISA = qw( BSML::BsmlElement );

BEGIN {
use BSML::BsmlSequenceAlignment;
use BSML::BsmlPairwiseAlignments;
use BSML::BsmlAlignmentSummary;
use BSML::BsmlCrossReference;
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
    $self->{'BsmlAlignmentSummaries'} = [];
    $self->{'BsmlPairwiseAlignments'} = [];
    $self->{'BsmlSequenceAlignments'} = [];
    $self->{'BsmlCrossReference'} = undef;
}

sub addBsmlAlignmentSummary
{
    my $self = shift;
    push( @{$self->{'BsmlAlignmentSummaries'}}, new BSML::BsmlAlignmentSummary );

    my $index = @{$self->{'BsmlAlignmentSummaries'}} - 1;
    return $index;
}

sub dropBsmlAlignmentSummary
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlAlignmentSummaries'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlAlignmentSummaries'}[$i] );
	  }
      }

    $self->{'BsmlAlignmentSummaries'} = \@newlist;    
  }

sub returnBsmlAlignmentSummaryListR
{
    my $self = shift;
    return $self->{'BsmlAlignmentSummaries'};
}

sub returnBsmlAlignmentSummaryR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlAlignmentSummaries'}[$index];
}

sub addBsmlPairwiseAlignments
{
    my $self = shift;
    push( @{$self->{'BsmlPairwiseAlignments'}}, new BSML::BsmlPairwiseAlignments );

    my $index = @{$self->{'BsmlPairwiseAlignments'}} - 1;
    return $index;
}

sub dropBsmlPairwiseAlignments
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlPairwiseAlignments'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlPairwiseAlignments'}[$i] );
	  }
      }

    $self->{'BsmlPairwiseAlignments'} = \@newlist;    
  }

sub returnBsmlPairwiseAlignmentsListR
{
    my $self = shift;
    return $self->{'BsmlPairwiseAlignments'};
}

sub returnBsmlPairwiseAlignmentsR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlPairwiseAlignments'}[$index];
}

sub addBsmlSequenceAlignment
{
    my $self = shift;
    push( @{$self->{'BsmlSequenceAlignments'}}, new BSML::BsmlSequenceAlignment );

    my $index = @{$self->{'BsmlSequenceAlignments'}} - 1;
    return $index;
}

sub dropBsmlSequenceAlignment
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlSequenceAlignments'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlSequenceAlignments'}[$i] );
	  }
      }

    $self->{'BsmlSequenceAlignments'} = \@newlist;    
  }

sub returnBsmlSequenceAlignmentListR
{
    my $self = shift;
    return $self->{'BsmlSequenceAlignments'};
}

sub returnBsmlSequenceAlignmentR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlSequenceAlignments'}[$index];
}

sub write
{
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Multiple-alignment-table", %{$self->{'attr'}} );

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $summary ( @{$self->{'BsmlAlignmentSummaries'}} )
    {
	$summary->write( $writer );
    }

    foreach my $pairwiseAlignments ( @{$self->{'BsmlPairwiseAlignments'}} )
    {
	$pairwiseAlignments->write( $writer );
    }

    foreach my $seqAln ( @{$self->{'BsmlSequenceAlignments'}} )
    {
	$seqAln->write( $writer );
    }

    foreach my $link (@{$self->{'BsmlLink'}})
    {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
    }

    if ( my $xref = $self->{'BsmlCrossReference'})
    {
	$xref->write( $writer );
    }

    $writer->endTag( "Multiple-alignment-table" );
}
    
1
