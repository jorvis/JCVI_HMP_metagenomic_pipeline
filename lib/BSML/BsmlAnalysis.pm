package BSML::BsmlAnalysis;
@ISA = qw( BSML::BsmlElement );

BEGIN {
use BSML::BsmlCrossReference;
use BSML::BsmlElement;
}
use XML::Writer;
use strict;
use warnings;

sub new 
  {
    my $class = shift;
    my ($logger_conf) = @_;
    my $self = {};
    bless $self, $class;
    
    $self->init( $logger_conf );
    return $self;
  }

sub init
  {
    my $self = shift;

    $self->{ 'attr' } = {};
    $self->{ 'BsmlAttr' } = {};
    $self->{ 'BsmlLink' } = [];
    $self->{ 'BsmlCrossReference' } = undef;

  }

sub write
  {
    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Analysis", %{$self->{'attr'}} );
       
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

     if ( my $xref = $self->{'BsmlCrossReference'})
    {
	$xref->write( $writer );
    } 


    $writer->endTag( "Analysis" );
  }

1
