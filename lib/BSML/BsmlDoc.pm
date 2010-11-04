package BSML::BsmlDoc;
@ISA = qw( BSML::BsmlElement );

=head1 NAME

BsmlDoc.pm - Bsml Document Class

=head1 VERSION

This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS

  Reading a BSML document and writing it back out from the object layer.

  my $doc = new BsmlDoc;
  my $parser = new BsmlParserTwig;
  $parser->parse( \$doc, $inputfile );
  $doc->write( $outputfile );

  Creating BSML objects and writing them out to a BSML file.

  my $doc = new BsmlDoc;
  my $seq = $doc->returnBsmlSequenceR( $doc->addBsmlSequence() );
  $seq->addBsmlSeqData( "ggtaccttctgaggcggaaagaaccagccggatccctcgaggg" );
  $doc->write();

=head1 DESCRIPTION

=head2 Overview

  This file provides a document level class for storing, creating, 
  and writing BSML elements

=head2 Constructor and initialization

When a BsmlDoc is created an empty document skeleton is created which currently supports 
document level attributes, BSML Attribute elements, and sequence data. Future support for 
Genome BSML elements as well as Research elements is expected at here.

If a path to a valid Log4Perl configuration file is passed to the constructor, the default
logging will be over-ridden.

=head2 Class and object methods

=over 4

=cut


use XML::Writer;
use warnings;

BEGIN {
use BSML::BsmlAnalysis;
use BSML::BsmlSequence;
use BSML::BsmlMultipleAlignmentTable;
use BSML::BsmlSegmentSet;
use BSML::BsmlGenome;
use BSML::BsmlSeqPairAlignment;
use BSML::BsmlElement;
}
use IO::File;

my $bsml_logger = BSML::Logger::get_logger( "Bsml" );


# The default links to the BSML dtd maintained by Labbook

my $default_dtd_pID = '-//EBI//Labbook, Inc. BSML DTD//EN';
my $default_dtd_sID = 'http://www.labbook.com/dtd/bsml3_1.dtd';

# The document level id lookup tables

$BsmlIdLookups = [];
$BsmlSeqAlignmentLookups = [];
$BsmlFeatureGroupLookups = [];
$BsmlTableIdCount = 0;
$BsmlCurrentTableId = 0;

# a bsml document stores a list of annotated sequences, 
# document level attributes, and Bsml Attribute Elements

sub new
  {
    my $class = shift;
    my ($logger_conf) = @_;
    my $self = {};
    $self->{'GZIP'}=0; #toggle for default gzipped output
    $self->{'EXT_DTD'}=0; #toggle for default use of external DTD at labbook
    bless $self, $class;    
    $self->init( $logger_conf );
    return $self;
  }

sub init
  {
    my $self = shift;
    my ($logger_conf) = @_;

    $self->{ 'attr' } = {};
    $self->{ 'BsmlAttr' } = {};
    $self->{ 'BsmlSequences' } = [];
    $self->{ 'BsmlSegmentSets' } = [];
    $self->{ 'BsmlSeqPairAlignments' } = [];
    $self->{ 'BsmlMultipleAlignmentTables' } = [];
    $self->{ 'BsmlAnalyses' } = [];
    $self->{ 'BsmlGenomes' } = [];

    $self->{'UID_COUNTER'} = 0;

    # initialize a namespace table
    $self->{ 'BsmlTableId' } = $BsmlTableIdCount;
    @{$BsmlIdLookups}[$BsmlTableIdCount] = {};
    @{$BsmlSeqAlignmentLookups}[$BsmlTableIdCount] = {};
    @{$BsmlFeatureGroupLookups}[$BsmlTableIdCount] = {};
    $BsmlTableIdCount++;

    $bsml_logger->info( "Created new BsmlDoc - BsmlDocID: $self->{'BsmlTableId'}" );
    $bsml_logger->level($DEBUG);
  }

sub DESTROY
  {
    my $self = shift;

    #free the memory associated with the global lookup table

    @{$BsmlIdLookups}[$self->{'BsmlTableId'}] = undef;
    @{$BsmlSeqAlignmentLookups}[$self->{'BsmlTableId'}] = undef;
    @{$BsmlFeatureGroupLookups}[$self->{'BsmlTableId'}] = undef;

    $bsml_logger->info( "Destructor Called - BsmlDocID: $self->{'BsmlTableId'}" );
    $bsml_logger->level($DEBUG);
  }

# return a unique, document level id. Counter is maintained as a 
# member variable of BsmlDoc

sub getUID
{
    return "Bsml".$self->{'UID_COUNTER'}++;
}

sub BsmlSetDocumentLookup
  {
    my ($BsmlID, $BsmlRef) = @_;
    $BsmlIdLookups->[$BsmlCurrentTableId]->{$BsmlID} = $BsmlRef;
  }

sub BsmlSetAlignmentLookup
  {
    my ($Seq1, $Seq2, $aln) = @_;

    if( !( $BsmlSeqAlignmentLookups->[$BsmlCurrentTableId]->{$Seq1}->{$Seq2} ) )
    {
	$BsmlSeqAlignmentLookups->[$BsmlCurrentTableId]->{$Seq1}->{$Seq2} = [];
    }

    push( @{$BsmlSeqAlignmentLookups->[$BsmlCurrentTableId]->{$Seq1}->{$Seq2}}, $aln );

  }

sub BsmlSetFeatureGroupLookup
  {
    my ($id, $fgroupid) = @_;
    
    if( (my $listref = $BsmlFeatureGroupLookups->[$BsmlCurrentTableId]->{$id}) )
      {
	push( @{$listref}, $fgroupid );
      }
    else
      {
	$BsmlFeatureGroupLookups->[$BsmlCurrentTableId]->{$id} = [$fgroupid];
      }
  }

sub BsmlReturnDocumentLookup
  {
    my ($BsmlID) = @_;
    return $BsmlIdLookups->[$BsmlCurrentTableId]->{$BsmlID};
  }

sub BsmlReturnAlignmentLookup
  {
    my ($Seq1, $Seq2) = @_;

    if( $Seq2 )
      {
	#returns a reference to a list of alignment objects
	return $BsmlSeqAlignmentLookups->[$BsmlCurrentTableId]->{$Seq1}->{$Seq2};
      }
    else
      {
	#returns a hash reference where each key specifies a reference to a list of alignment objects
	return $BsmlSeqAlignmentLookups->[$BsmlCurrentTableId]->{$Seq1};
      }
  }

sub BsmlReturnFeatureGroupLookup
  {
    my ($BsmlId) = @_;
    return $BsmlFeatureGroupLookups->[$BsmlCurrentTableId]->{$BsmlId};
  }

sub BsmlReturnFeatureGroupLookupIds
  {
    return keys( %{$BsmlFeatureGroupLookups->[$BsmlCurrentTableId]} );
  }

sub BsmlSetCurrentDocumentLookupTable
  {
    my ($index) = @_;
    $BsmlCurrentTableId = $index;
  }

sub makeCurrentDocument
  {
    my $self = shift;
    BsmlSetCurrentDocumentLookupTable( $self->{'BsmlTableId'} );
  }

=item $doc->addBsmlSequence()

B<Description:> This is a method to add a Bsml Sequence Object to the document.

B<Parameters:> None - A $ref parameter will probably be added to allow a sequence reference to be returned directly by the function.

B<Returns:> The reference index of the added sequence

=cut

sub addBsmlSequence
  {
    my $self = shift;
     
    push( @{$self->{'BsmlSequences'}}, new BSML::BsmlSequence );

    my $index = @{$self->{'BsmlSequences'}} - 1;

    $bsml_logger->info( "Added BsmlSequence: $index" );

    return $index;    
  }

sub addBsmlSegmentSet
{
    my $self = shift;
    
    push( @{$self->{'BsmlSegmentSets'}}, new BSML::BsmlSegmentSet );

    my $index = @{$self->{'BsmlSegmentSets'}} - 1;

    $bsml_logger->info( "Added BsmlSegmentSet: $index" );

    return $index;
}

=item $doc->dropBsmlSequence()

B<Description:> Delete a Bsml Sequence from the document.

B<Parameters:> ($index) - the sequence index returned from addBsmlSequence (position of the sequence in the reference list)

B<Returns:> None

=cut

sub dropBsmlSequence
  {
    my $self = shift;

    my ($index) = @_;

    my $newlist;

    for(  my $i=0;  $i< @{$self->{'BsmlSequences'}}; $i++ ) 
      {
	if( $i != $index )
	  {
	    push( @{$newlist}, $self->{'BsmlSequences'}[$i] );
	  }
      }

    $self->{'BsmlSequences'} = $newlist;

    $bsml_logger->info( "Dropped BsmlSequence: $index" );
  }

sub dropBsmlSegmentSet
  {
    my $self = shift;

    my ($index) = @_;

    my $newlist;

    for(  my $i=0;  $i< @{$self->{'BsmlSegmentSets'}}; $i++ ) 
      {
	if( $i != $index )
	  {
	    push( @{$newlist}, $self->{'BsmlSegmentSets'}[$i] );
	  }
      }

    $self->{'BsmlSegmentSets'} = $newlist;

    $bsml_logger->info( "Dropped BsmlSegmentSet: $index" );
  }

=item $doc->returnBsmlSequenceListR()

B<Description:> Return a list of references to all the sequence objects contained in the document.

B<Parameters:> None

B<Returns:> a list of BsmlSequence object references

=cut 

sub returnBsmlSequenceListR
  {
    my $self = shift;

    return $self->{'BsmlSequences'};
  }

sub returnBsmlSegmentSetListR
{
    my $self = shift;
    return $self->{'BsmlSegmentSets'};
}

=item $doc->returnBsmlSequenceR()

B<Description:> Return a reference to a sequence object given its index

B<Parameters:> ($index) - the sequence index returned from addBsmlSequence (position of the sequence in the reference list)

B<Returns:> a BsmlSequence object reference

=cut

sub returnBsmlSequenceR
  {
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlSequences'}[$index];  
  }

sub returnBsmlSegmentSetR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlSegmentSets'}[$index];
}

sub returnBsmlSequenceByIDR
  {
    my $self = shift;
    my ($id) = @_;

    return $BsmlIdLookups->[$BsmlCurrentTableId]->{$id};
  }

sub addBsmlSeqPairAlignment
  {
    my $self = shift;
     
    push( @{$self->{'BsmlSeqPairAlignments'}}, new BSML::BsmlSeqPairAlignment );

    my $index = @{$self->{'BsmlSeqPairAlignments'}} - 1;

    $bsml_logger->info( "Added BsmlSeqPairAlignment: $index" );

    return $index;    
  }

sub dropBsmlSeqPairAlignment
{
   my $self = shift;
   my ($index) = @_;

   my $newlist;
   for(  my $i=0;  $i< @{$self->{'BsmlSeqPairAlignments'}}; $i++ ) 
      {
	if( $i != $index )
	  {
	    push( @{$newlist}, $self->{'BsmlSeqPairAlignments'}[$i] );
	  }
      }

    $self->{'BsmlSeqPairAlignements'} = $newlist;

    $bsml_logger->info( "Dropped BsmlSeqPairAlignment: $index" );
}

sub returnBsmlSeqPairAlignmentListR
{
  my $self = shift;
  return $self->{'BsmlSeqPairAlignments'};
}

sub returnBsmlSeqPairAlignmentR
{
  my $self = shift;
  my ($index) = @_;

  return $self->{'BsmlSeqPairAlignments'}[$index];
}

sub addBsmlMultipleAlignmentTable
{
    my $self = shift;
    push( @{$self->{'BsmlMultipleAlignmentTables'}}, new BSML::BsmlMultipleAlignmentTable );

    my $index = @{$self->{'BsmlMultipleAlignmentTables'}} - 1;
    return $index;
}

sub dropBsmlMultipleAlignmentTable
  {
    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlMultipleAlignmentTables'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlMultipleAlignmentTables'}[$i] );
	  }
      }

    $self->{'BsmlMultipleAlignmentTables'} = \@newlist;    
  }

sub returnBsmlMultipleAlignmentTableListR
{
    my $self = shift;
    return $self->{'BsmlMultipleAlignmentTables'};
}

sub returnBsmlMultipleAlignmentTableR
{
    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlMultipleAlignmentTables'}[$index];
}


sub addBsmlAnalysis
  {
    my $self = shift;
     
    push( @{$self->{'BsmlAnalyses'}}, new BSML::BsmlAnalysis );

    my $index = @{$self->{'BsmlAnalyses'}} - 1;

    $bsml_logger->info( "Added BsmlAnalysis: $index" );

    return $index;    
  }

sub dropBsmlAnalysis
{
   my $self = shift;
   my ($index) = @_;

   my $newlist;
   for(  my $i=0;  $i< @{$self->{'BsmlAnalyses'}}; $i++ ) 
      {
	if( $i != $index )
	  {
	    push( @{$newlist}, $self->{'BsmlAnalyses'}[$i] );
	  }
      }

    $self->{'BsmlAnalyses'} = $newlist;

    $bsml_logger->info( "Dropped BsmlAnalyses: $index" );
}

sub returnBsmlAnalysisListR
{
  my $self = shift;

  return $self->{'BsmlAnalyses'};
}

sub returnBsmlAnalysisR
{
  my $self = shift;
  my ($index) = @_;

  return $self->{'BsmlAnalyses'}[$index];
}

sub addBsmlGenome
{
    my $self = shift;

    push( @{$self->{'BsmlGenomes'}}, new BSML::BsmlGenome );

    my $index = @{$self->{'BsmlGenomes'}} - 1;

    return $index;    
}

sub returnBsmlGenomeListR
{
  my $self = shift;
  return $self->{'BsmlGenomes'};
}

sub returnBsmlGenomeR
{
  my $self = shift;
  my ($index) = @_;

  return $self->{'BsmlGenomes'}[$index];
}

sub dropBsmlGenome
{
   my $self = shift;
   my ($index) = @_;

   my $newlist;
   for(  my $i=0;  $i< @{$self->{'BsmlGenomes'}}; $i++ ) 
      {
	if( $i != $index )
	  {
	    push( @{$newlist}, $self->{'BsmlGenomes'}[$i] );
	  }
      }

    $self->{'BsmlGenomes'} = $newlist;
}


=item $doc->write()

B<Description:> Writes the document to a file or stream

B<Parameters:> ($fileOrHandle, $dtd) - output file name or the string 'STDOUT' or a file handle, optional user specified dtd which will override the librarys default

B<Returns:> None

=cut

sub write
  {
    my $self = shift;
    my ($fileOrHandle, $dtd, $gzip) = @_;

    $bsml_logger->debug( "Attempting to write BsmlDoc" );
    my $output;

    # Maintain reverse compatibility; allow use of string 'STDOUT' as a filename
    if ($fileOrHandle eq 'STDOUT') 
    {
	$output = \*STDOUT;
    }
    # $fileOrHandle is a handle or glob
    elsif (ref($fileOrHandle) && ($fileOrHandle->isa("IO::Handle") || $fileOrHandle->isa("GLOB"))) 
    {
	$output = $fileOrHandle;
    }
    # otherwise assume $fileOrHandle is a filename
    else 
    {
	if($gzip || $self->{'GZIP'}){
	    $fileOrHandle .= ".gz" unless ($fileOrHandle =~ /\.gz$/);
	    open $output,">:gzip","$fileOrHandle" or die "could not open gzip output stream for - $fileOrHandle $!\n";;
	}
	else{
	    $output = new IO::File( ">$fileOrHandle" ) or die "could not open output file - $fileOrHandle $!\n";;
	}

	if( !( $output ) )
	  {
	    $bsml_logger->fatal( "Could not open output file - $fileOrHandle $!" );
	    die "could not open output file - $fileOrHandle $!\n";
	  }
      }

    # Setting DATA_MODE to 1 enables XML::Writer to insert newlines around XML elements for 
    # easier readability. DATA_INDENT specifies an indent of two spaces for child elements

    my $writer = new XML::Writer(OUTPUT => $output, DATA_MODE => 1, DATA_INDENT => 2);

    $writer->xmlDecl();

    # If a local file dtd is not specified the default will be used.
    # The default points to the publicly available dtd at Labbook

    if( $dtd ){$writer->doctype( "Bsml", "", "file:$dtd" );}
    elsif($self->{'EXT_DTD'}){
      $writer->doctype( "Bsml", $default_dtd_pID, $default_dtd_sID );
      $bsml_logger->debug( "DTD not specified - using $default_dtd_sID" );
    }
   
    # write the root node

    $writer->startTag( "Bsml", %{$self->{'attr'}} );

    # write any Bsml Attribute children

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    # write the Definitions section; current API only supports Sequences

    $writer->startTag( "Definitions" );

    if( @{$self->{'BsmlGenomes'}} )
    {
	$writer->startTag( "Genomes" );
	
	foreach my $genome ( @{$self->{'BsmlGenomes'}} )
	{
	    $genome->write( $writer );
	}
	
	$writer->endTag( "Genomes" );
    }

    # write the sequence and segment-set elements

    if( @{$self->{'BsmlSequences'}} || @{$self->{'BsmlSegmentSets'}} )
    {
	$writer->startTag( "Sequences" );

	foreach my $seq ( @{$self->{'BsmlSequences'}} )
	{
	    $seq->write( $writer );
	}

	foreach my $segmentSet ( @{$self->{'BsmlSegmentSets'}} )
	{
	    $segmentSet->write( $writer );
	}

	$writer->endTag( "Sequences" );
    }

    if( @{$self->{'BsmlSeqPairAlignments'}} || @{$self->{'BsmlMultipleAlignmentTables'}} )
    {
      $writer->startTag( "Tables", 'id' => 'BsmlTables' );

      foreach my $seqAlignment ( @{$self->{'BsmlSeqPairAlignments'}} )
      {
        $seqAlignment->write( $writer );
      }

      foreach my $malnTable ( @{$self->{'BsmlMultipleAlignmentTables'}} )
      {
	  $malnTable->write( $writer );
      }

      $writer->endTag( "Tables" );
    }

    $writer->endTag( "Definitions" );

    if( @{$self->{'BsmlAnalyses'}} )
      {
	$writer->startTag( 'Research' );
	$writer->startTag( 'Analyses' );
	
	foreach my $analysis ( @{$self->{'BsmlAnalyses'}} )
	  {
	    $analysis->write( $writer );
	  }

	$writer->endTag( 'Analyses' );
	$writer->endTag( 'Research' );
      }
	
    
    $writer->endTag( "Bsml" );

    #clean up open fhs
    $writer->end();
    
    # NOTE: one might also want to omit the call to close() if $output->isa("IO::String")
    if(($output != \*STDOUT) && ($output != \*STDERR))
      {
	  $output->close();
	  $bsml_logger->info("Output handle $output closed for file/handle $fileOrHandle $!");
      }

    $bsml_logger->debug( "BsmlDoc successfully written" );
  }

1
