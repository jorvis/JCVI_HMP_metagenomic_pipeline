package BSML::BsmlFeature;
@ISA = qw( BSML::BsmlElement );

use XML::Writer;
BEGIN {
use BSML::BsmlCrossReference;
use BSML::BsmlElement;
use BSML::Logger;
}
use strict;
use warnings;
use Data::Dumper;

my $logger = BSML::Logger::get_logger("Logger::BSML");


=head1 NAME

BsmlFeature.pm - The Bsml Feature class encodes feature qualifier and feature location elements

=head1 VERSION

This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS

  Adding a feature qualifier and site location to a feature created with BsmlFeatureTable

  my $FTable = new BsmlFeatureTable;
  my $feature = $Ftable->returnBsmlFeatureR( $Ftable->addBsmlFeature() );
  $feature->addattr( 'id', '_1396a' );
  $feature->addBsmlQualifier( 'gene', 'pf14' );
  $feature->addBsmlSiteLoc( '1002' );
  
=head1 DESCRIPTION

=head2 Overview

  This file implements a class providing support for BSML Feature elements and their children
  including Bsml Qualifier, Site-Loc, and Interval-Loc elements.

=head2 Constructor and initialization

  Typically a BsmlFeature is created by the BsmlFeatureTable it is contained within and manipulated
  as a reference.

=head2 Class and object methods

=over 4

=cut

sub new
  {

      $logger->debug("Instantiating BSML::BsmlFeature") if $logger->is_debug;

    my $class = shift;
    my $self = {};
    bless $self, $class;
    
    $self->init();
    return $self;
  }

sub init
  {

      $logger->debug("Init'ing BSML::BsmlFeature") if $logger->is_debug;

    my $self = shift;

    $self->{ 'attr' } = {};
    $self->{ 'BsmlAttr' } = {};
    $self->{ 'BsmlAttributeList' } = undef;
    $self->{ 'BsmlSite-Loc' } = [];
    $self->{ 'BsmlInterval-Loc' } = [];
    $self->{ 'BsmlQual' } = {};
    $self->{ 'BsmlLink' } = [];
    $self->{ 'BsmlCrossReference' } = [];
  }
    
=item $feature->addBsmlQualifier( $valuetype, $value )

B<Description:> adds a Qualifier element to the feature

B<Parameters:> ( $valuetype, $value ) - the key value pair describing the Qualifer

B<Returns:> None

=cut 

sub addBsmlQualifier
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $valuetype, $value ) = @_;

    $self->{'BsmlQual'}->{ $valuetype } = $value;
  }

=item $feature->setBsmlQualifier( $valuetype, $value )

B<Description:> sets a Qualifier element on a feature

B<Parameters:> ( $valuetype, $value ) - the key value pair describing the Qualifer

B<Returns:> None

=cut

sub setBsmlQualifier
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $valuetype, $value ) = @_;

    $self->addBsmlQualifier( $valuetype, $value );
  }

=item $feature->setBsmlQualifiersHash( $hashref )

B<Description:> sets Qualifiers according to the key/value pairs (valuetype, value) defined in the hash pointed to by $hashref

B<Parameters:> ($hashref) a reference pointing to a hash containing the key value pairs to be added as qualifiers

B<Returns:> None

sub setBsmlQualifiersHash
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    my ($href) = @_;

    foreach my $key (keys(%{$href}))
      {
	$self->addBsmlQualifier( $key, $href->{$key} );
      }
  }

=item $feature->dropBsmlQualifier( $valuetype )

B<Description:> deletes a Qualifier element from a feature

B<Parameters:> ( $valuetype ) - the key describing the Qualifer

B<Returns:> None

=cut

sub dropBsmlQualifier
  {


      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $valuetype ) = @_;

    delete( $self->{'BsmlQual'}->{ $valuetype } );
  }

=item $feature->returnBsmlQualifierHash( )

B<Description:> returns a hash reference to the key value pairs representing a features qualifiers

B<Parameters:> None

B<Returns:> hash reference

=cut

sub returnBsmlQualifierHashR
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    return $self->{'BsmlQual'};
  }

=item $feature->returnBsmlQualifier( $key )

B<Description:> return a the value associated with a Qualifiers key

B<Parameters:> ($key) - Qualifier key

B<Returns:> scaler

=cut

sub returnBsmlQualifier
{

      $logger->debug("") if $logger->is_debug;

  my $self = shift;
  my ($key) = @_;

  return $self->{'BsmlQual'}->{$key};
}

=item $feature->addBsmlSiteLoc( $position, $complement )

B<Description:> add a Site-Loc tag to the feature

B<Parameters:> ($position, $complement) - sequence position of the feature, optional 0 or 1 reflecting the strand containing the feature  

B<Returns:> None

=cut

sub addBsmlSiteLoc
  {


      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $sitepos, $comp, $class ) = @_;

    if( !( $comp ) ){ $comp = 0; }

    my $href = {'sitepos' => $sitepos, 'complement' => $comp};
    if( $class ){ $href->{'class'} = $class;}


    push( @{$self->{'BsmlSite-Loc'}}, $href );
  }

=item $feature->setBsmlSiteLoc( $position, $complement )

B<Description:> add a Site-Loc tag to the feature

B<Parameters:> ($position, $complement) - sequence position of the feature, optional 0 or 1 reflecting the strand containing the feature  

B<Returns:> None

=cut

sub setBsmlSiteLoc
  {

      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    my ( $sitepos, $comp, $class ) = @_;

    my @newlist;

    foreach my $site ( @{$self->{'BsmlSite-Loc'}} )
      {
	if( !($site->{'sitepos'} == $sitepos ) )
	  {
	    push( @newlist, $site );
	  }
      }

    $self->{'BsmlSite-Loc'} = \@newlist;

    $self->addBsmlSiteLoc( $sitepos, $comp, $class );
  }

=item $feature->dropBsmlSiteLoc( $position  )

B<Description:> drop a Site-Loc tag at a given location

B<Parameters:> ($position) - the sequence position of the Site-Loc to be dropped 

B<Returns:> None

=cut

sub dropBsmlSiteLoc
  {

      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    my ( $sitepos, $comp ) = @_;

    my $newlist = [];
    
    foreach my $site ( @{$self->{'BsmlSite-Loc'}} )
      {
	if( !($site->{'sitepos'} == $sitepos ) )
	  {
	    push( @{$newlist}, $site );
	  }
      }

    $self->{'BsmlSite-Loc'} = $newlist;
  }

=item $feature->returnBsmlSiteLocList( )

B<Description:> returns a reference to a list containing the Site-Loc elements of the feature

B<Parameters:> None

B<Returns:> list reference

=cut

sub returnBsmlSiteLocListR
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    return $self->{'BsmlSite-Loc'};
  }

# see bug 2307
# Added $startopen and $endopen

=item $feature->addBsmlIntervalLoc( $startpos, $endpos, $complement, $startopen, $endopen )

B<Description:> add a Interval-Loc tag to the feature

B<Parameters:> ($startpos, $endpos, $complement, $startopen, $endopen) - sequence start position, sequence end position, optional 0 or 1 reflecting the strand containing the feature, 0 or 1 reflecting if the actual start position is unbounded lower than $startpos, 0 or 1 reflecting if the actual end position is unbounded higher than $endpos

B<Returns:> None

=cut

sub addBsmlIntervalLoc
  {

      $logger->debug("") if $logger->is_debug;
      
    my $self = shift;
    my ( $startpos, $endpos, $comp, $startopen, $endopen ) = @_;

    if( !( $comp ) ){ $comp = 0; }

    push( @{$self->{'BsmlInterval-Loc'}}, { 'startpos' => $startpos, 'endpos' => $endpos, 'complement' => $comp } );
    # only add startopen or endopen if they're true to cut down on the size of the BSML
    # they are added to the last hash added to @{$self->{'BsmlInterval-Loc'}}
    if ($startopen) {
	my $idx = $#{$self->{'BsmlInterval-Loc'}};
	${$self->{'BsmlInterval-Loc'}}[$idx]->{'startopen'} = $startopen;
    }
    if ($endopen) {
	my $idx = $#{$self->{'BsmlInterval-Loc'}};
	${$self->{'BsmlInterval-Loc'}}[$idx]->{'endopen'} = $endopen;
    }
  }


# see bug 2307
# Added $startopen and $endopen

=item $feature->setBsmlIntervalLoc( $startpos, $endpos, $complement, $startopen, $endopent )

B<Description:> add a Interval-Loc tag to the feature

B<Parameters:> ($startpos, $endpos, $complement, $startopen, $endopen) - sequence start position, sequence end position, optional 0 or 1 reflecting the strand containing the feature, 0 or 1 reflecting if the actual start position is unbounded lower than $startpos, 0 or 1 reflecting if the actual end position is unbounded higher than $endpos

B<Returns:> None

=cut

sub setBsmlIntervalLoc
  {

      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    my ( $startpos, $endpos, $comp, $startopen, $endopen ) = @_;

    my @newlist;
    
    foreach my $site ( @{$self->{'BsmlInterval-Loc'}} )
      {
	if( !($site->{'startpos'} == $startpos ) && !( $site->{'endpos'} == $endpos) )
	  {
	    push( @newlist, $site );
	  }
      }

    $self->{'BsmlInterval-Loc'} = \@newlist;

    $self->addBsmlIntervalLoc( $startpos, $endpos, $comp, $startopen, $endopen );
  }




=item $feature->dropBsmlIntervalLoc( $startpos, $endpos  )

B<Description:> drop an Interval-Loc tag from the feature

B<Parameters:> ($startpos, $endpos ) - sequence start position, sequence end position  

B<Returns:> None

=cut

sub dropBsmlIntervalLoc
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $startpos, $endpos, $comp ) = @_;

    my $newlist = [];
    
    foreach my $site ( @{$self->{'BsmlInterval-Loc'}} )
      {
	if( !($site->{'startpos'} == $startpos ) && !( $site->{'endpos'} == $endpos) )
	  {
	    push( @{$newlist}, $site );
	  }
      }

    $self->{'BsmlInterval-Loc'} = $newlist;
  }

=item $feature->returnBsmlIntervalLocList( )

B<Description:> return a reference to a list containing the features interval locations

B<Parameters:> None

B<Returns:> reference to a list of hash references

=cut

sub returnBsmlIntervalLocListR
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    return $self->{'BsmlInterval-Loc'};
  }

=item $seq->write()

  B<Description:> writes the BSML elements encoded by the class to a file using XML::Writer. This method should only be called through the BsmlDoc->write() process.

  B<Parameters:> None

  B<Returns:> None

=cut 

# #
# # Added 2004-11-02 Bugzilla case 1808
# #
# sub addBsmlAttributeList {

#     my $self = shift;
#     push ( @{$self->{'BsmlAttributeList'}}, new BSML::BsmlAttributeList );
    
#     my $index = @{$self->{'BsmlAttributeList'}} - 1;
#     return $index;
# }

# sub dropBsmlAttributeList {

#     my $self = shift;
#     my ($index) = @_;

#     my $newlist = [];

#     for(  my $i=0;  $i< @{$self->{'BsmlAttributeList'}}; $i++ ) {
# 	if( $i != $index ){
# 	    push( @{$newlist}, $self->{'BsmlAttributeList'}[$i] );
# 	}
#     }

#     $self->{'BsmlAttributeList'} = $newlist;
# }

# sub returnBsmlAttributeListoListR {

#     my $self = shift;
#     return $self->{'BsmlAttributeList'};

# }

# sub returnBsmlAttributeListR {

#     my $self = shift;
#     my ($index) = @_;

#     return $self->{'BsmlAttributeList'}[$index];
# }





sub write
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my $writer = shift;

    $writer->startTag( "Feature", %{$self->{'attr'}} );

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    foreach my $bsmlqual (keys( %{$self->{'BsmlQual'}} ))
      {
	my %tmp = ( 'value-type' => $bsmlqual, 'value' => $self->{'BsmlQual'}->{$bsmlqual} );
	
	$writer->startTag( "Qualifier", %tmp );
	$writer->endTag( "Qualifier" );
      }

    foreach my $bsmlsite (@{$self->{'BsmlSite-Loc'}})
      {
	$writer->startTag( "Site-loc", %{$bsmlsite} );
	$writer->endTag( "Site-loc" );
      }

    foreach my $bsmlsite (@{$self->{'BsmlInterval-Loc'}})
      {
	$writer->startTag( "Interval-loc", %{$bsmlsite} );
	$writer->endTag( "Interval-loc" );
      }


    foreach my $xref (@{$self->{'BsmlCrossReference'}}){
	$xref->write( $writer );
    }
    
    foreach my $link (@{$self->{'BsmlLink'}})
      {
        $writer->startTag( "Link", %{$link} );
        $writer->endTag( "Link" );
      }

#    use Data::Dumper;

#    print Dumper $self->{'BsmlAttributeList'};
   
  foreach my $listref (@{$self->{'BsmlAttributeList'}}){
	  $writer->startTag( "Attribute-list");
	    foreach my $hash ( @{$listref} ){ 
		    $writer->startTag( "Attribute",  'name' => $hash->{'name'}, 'content' => $hash->{'content'} );
		    $writer->endTag( "Attribute" );
		}
	   $writer->endTag( "Attribute-list" );
      }
        
    
    $writer->endTag( "Feature" );
}

1
