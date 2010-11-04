package BSML::BsmlParserTwig;
use Data::Dumper;

=head1 NAME

  BsmlParserTwig.pm - Bsml Parsing Class built on XML::Twig

=head1 VERSION

This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS

  Parsing a BSML Document into the BSML Object Layer

  my $doc = new BsmlDoc;
  my $parser = new BsmlParserTwig;

  $parser->parse( $doc, $filenameOrHandle );

=head1 DESCRIPTION

=head2 Overview

  This file provides a parsing class written on top of XML::Twig to populate 
  a BsmlDoc with data provided in a BSML file. XML::Twig is used to process 
  subtrees at the Sequence level. XML::Twig does not apply XML validation.

=head2 Constructor and initialization

  my $parser = new BsmlParserTwig;

=head2 Class and object methods

=over 4

=cut

use strict;
use warnings;
use XML::Twig;
BEGIN {
use BSML::BsmlDoc;
use BSML::Logger;
}
use Data::Dumper;

my $bsmlDoc;

sub new
  {
    my $class = shift;
    my $self = {};
    bless $self, $class;
    
    return $self;
  }

=item $parser->parse( $bsml_doc, $filename )

B<Description:> parse the contents of $fileOrHandle into the BSML Document Object $bsml_doc

B<Parameters:> ($bsml_doc, $fileOrHandle) - a BsmlDoc object, a filename or file handle pointing to a BSML document

B<Returns:> None

=cut

sub parse
  {
    my $self = shift;
    my ( $bsml_doc, $fileOrHandle ) = @_;
    my $bsml_logger = BSML::Logger::get_logger( "Bsml" );

    $bsmlDoc = ${$bsml_doc}; 

    if( !( $fileOrHandle ) ){   
      $bsml_logger->fatal( "Filename or file handle not provided for Bsml Parsing." );
    }

    if( !( $bsmlDoc ) ){
      $bsml_logger->fatal( "No BsmlDoc object to populate" );
    }

    #insure the lookup tables are using the right namespace
    $bsmlDoc->makeCurrentDocument();

    # Set a Twig Handler on the BSML Sequence Object. The handler is called
    # each time an XML subtree rooted at a Sequence element is completely
    # parsed

    my $twig = new XML::Twig( TwigHandlers => 
			  { 'Sequence' => \&sequenceHandler, 
			    'Seq-pair-alignment' => \&seqPairAlignmentHandler, 
			    'Analysis' => \&analysisHandler, 
			    'Multiple-alignment-table' => \&multiAlignmentTableHandler, 
			    'Genome' => \&genomeHandler }
			  );
    
    # parsefile will die if an xml syntax error is encountered or if
    # there is an io problem

    $bsml_logger->debug( "Attempting to Parse Bsml Document: $fileOrHandle" );
    if (ref($fileOrHandle) && ($fileOrHandle->isa("IO::Handle") || $fileOrHandle->isa("GLOB"))) {
	$twig->parse( $fileOrHandle );
    } else {
	$twig->parsefile( $fileOrHandle );
    }
    $bsml_logger->info( "Successfully Parsed Bsml Document: $fileOrHandle" );

    # Twig documentation claims circular references in the twig class prevent garbage collection. Could
    # this be the source of the file handle problem? It seems like the 'twig' will go out of scope reguardlessly.

    $twig->dispose();

    # Don't want to keep the extra reference around in a class variable.

    $bsmlDoc = undef;
  }

# This is a private method implemented as an XML::Twig handler object. It is 
# called each time XML::Twig successfully and completely parses the subtree
# rooted at a Sequence element. The primary parse tree should be purged each time this
# method completes to handle memory efficiently. This method will need to be extended 
# to add support for additional BSML elements as they are chosen to be incorporated into
# the api.

sub sequenceHandler
  {
    my $bsml_logger = BSML::Logger::get_logger( "Bsml" );
    $bsml_logger->info( "Parsing Sequence Twig" );

    my ($twig, $seq) = @_;

    # add a new Sequence object to the bsmlDoc

    my $bsmlseq = $bsmlDoc->{'BsmlSequences'}[$bsmlDoc->addBsmlSequence()];
    
    # add the sequence element's attributes

    my $attr = $seq->atts();

    $bsml_logger->debug("Processing all //Sequence/\@attributes") if $bsml_logger->is_debug;

    foreach my $key ( keys( %{$attr} ) )
      {
	$bsmlseq->addattr( $key, $attr->{$key} );
      }

    #put the sequence into the general object lookup table if an id was specified

    if( my $id = $bsmlseq->returnattr('id')){
      BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $bsmlseq );}
      
    # add Sequence level Bsml Attribute elements 

    $bsml_logger->debug("Processing all //Sequence/Attribute elements") if $bsml_logger->is_debug;

    foreach my $BsmlAttr ( $seq->children( 'Attribute' ) )
      {
	my $attr = $BsmlAttr->atts();
	$bsmlseq->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
      }

    # add any Bsml Links

    $bsml_logger->debug("Processing all //Sequence/Link elements") if $bsml_logger->is_debug;

    foreach my $BsmlLink ( $seq->children( 'Link' ) )
      {
	my $attr = $BsmlLink->atts();
	$bsmlseq->addBsmlLink( $attr->{'title'}, $attr->{'href'} );
      }

    # add raw sequence data if found 

    $bsml_logger->debug("Processing all //Sequence/Seq-data elements") if $bsml_logger->is_debug;

    my $seqDat = $seq->first_child( 'Seq-data' );
    
    if( $seqDat ){
      $bsmlseq->addBsmlSeqData( $seqDat->text() );
    }

    # add appended sequence data if found

    $bsml_logger->debug("Processing all //Sequence/Seq-data-import elements") if $bsml_logger->is_debug;

    my $seqDatImport = $seq->first_child( 'Seq-data-import' );
    
    if( $seqDatImport )
      {
	my $attr = $seqDatImport->atts();
	$bsmlseq->addBsmlSeqDataImport( $attr->{'format'}, $attr->{'source'}, $attr->{'id'},$attr->{'identifier'});
      }

    # add numbering information

    $bsml_logger->debug("Processing all //Sequence/Numbering/ elements") if $bsml_logger->is_debug;

    my $numbering = $seq->first_child( 'Numbering' );

    if( $numbering )
    {
	my $attr = $numbering->atts();
	my $bsmlNumbering = $bsmlseq->returnBsmlNumberingR( $bsmlseq->addBsmlNumbering() );

	$bsmlNumbering->addattr( 'seqref', $attr->{'seqref'} );
	$bsmlNumbering->addattr( 'use-numbering', $attr->{'use_numbering'} );
	$bsmlNumbering->addattr( 'type', $attr->{'type'} );
	$bsmlNumbering->addattr( 'units', $attr->{'units'} );
	$bsmlNumbering->addattr( 'a', $attr->{'a'} );
	$bsmlNumbering->addattr( 'b', $attr->{'b'} );
	$bsmlNumbering->addattr( 'dec-places', $attr->{'dec_places'} );
	$bsmlNumbering->addattr( 'refnum', $attr->{'refnum'} );
	$bsmlNumbering->addattr( 'has-zero', $attr->{'has_zero'} );
	$bsmlNumbering->addattr( 'ascending', $attr->{'ascending'} );
	$bsmlNumbering->addattr( 'names', $attr->{'names'} );
	$bsmlNumbering->addattr(' from-aligns', $attr->{'from_aligns'} );
	$bsmlNumbering->addattr( 'aligns', $attr->{'aligns'} );


	$bsml_logger->debug("Processing //Sequence/Numbering/Attribute elements") if $bsml_logger->is_debug;

	foreach my $BsmlAttr ( $numbering->children( 'Attribute' ) )
	{
	    my $attr = $BsmlAttr->atts();
	    $bsmlNumbering->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
	}
    }

    # add Feature Tables with Feature, Reference, and Feature-group Elements

    my $BsmlFTables = $seq->first_child( 'Feature-tables' );

    if( $BsmlFTables )
      {
	  $bsml_logger->debug("Processing first_child //Feature-table/ element") if $bsml_logger->is_debug;

	  foreach my $BsmlFTable ( $BsmlFTables->children( 'Feature-table' ) )
	  {
	      my $table = $bsmlseq->{'BsmlFeatureTables'}[$bsmlseq->addBsmlFeatureTable()];
	      
	      my $attr = $BsmlFTable->atts();
	      
	      $bsml_logger->debug("Processing //Feature-table/\@attibutes") if $bsml_logger->is_debug;

	      foreach my $key ( keys( %{$attr} ) )
	      {
		  $table->addattr( $key, $attr->{$key} );
	      }
	      
	      #if an id has been specified, add the table to the general object lookups
	      if( my $id = $table->returnattr('id')){
		  BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $table );}
	      
	      $bsml_logger->debug("Processing //Feature-table/Link elements") if $bsml_logger->is_debug;

	      foreach my $BsmlLink ( $BsmlFTable->children( 'Link' ) )
	      {
		  my $attr = $BsmlLink->atts();
		  $table->addBsmlLink( $attr->{'title'}, $attr->{'href'} );
	      }
	      
	      $bsml_logger->debug("Processing //Feature-table/Reference/ elements") if $bsml_logger->is_debug;

	      foreach my $BsmlRef ($BsmlFTable->children( 'Reference' ) )
	      {
		my $ref = $table->{'BsmlReferences'}[$table->addBsmlReference()];
		my $attr = $BsmlRef->atts();

		foreach my $key ( keys( %{$attr} ) )
		  {
		    $ref->addattr( $key, $attr->{$key} );
		  }

                #if an id has been specified, add the reference to the general object lookups
	        if( my $id = $ref->returnattr('id')){
	            BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $ref );}

		foreach my $BsmlAuthor ($BsmlRef->children( 'RefAuthors' ))
		  {
		    $ref->addBsmlRefAuthors( $BsmlAuthor->text() );
		  }

		foreach my $BsmlRefTitle ($BsmlRef->children( 'RefTitle' ))
		  {
		    $ref->addBsmlRefTitle( $BsmlRefTitle->text() );
		  }

		foreach my $BsmlRefJournal ($BsmlRef->children( 'RefJournal' ))
		  {
		    $ref->addBsmlRefJournal( $BsmlRefJournal->text() );
		  }

		foreach my $BsmlLink ( $BsmlRef->children( 'Link' ) )
		  {
		    my $attr = $BsmlLink->atts();
		    $ref->addBsmlLink( $attr->{'title'}, $attr->{'href'} );
		  }
	      }
	  

	      $bsml_logger->debug("Attempting to process //Feature-table/Feature/ elements") if $bsml_logger->is_debug;

	      foreach my $BsmlFeature ($BsmlFTable->children( 'Feature' ))
	      {
		  my $feat = $table->{'BsmlFeatures'}[$table->addBsmlFeature()];
		  my $attr = $BsmlFeature->atts();

		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/\@attributes") if $bsml_logger->is_debug;
		  
		  foreach my $key ( keys( %{$attr} ) )
		  {
		      $feat->addattr( $key, $attr->{$key} );
		  }
		  
		 #if an id has been specified, add the feature to the general object lookups
	         if( my $id = $feat->returnattr('id')){
	           BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $feat );}


		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/Cross-reference elements") if $bsml_logger->is_debug;

		  foreach my $crossRef( $BsmlFeature->children( 'Cross-reference' )) {
		    
		      my $attr = $crossRef->atts();
		      
		      my $xref = $feat->returnBsmlCrossReferenceR( $feat->addBsmlCrossReference );
		    
		      foreach my $key ( keys( %{$attr} ) ) {
			  
			  $xref->addattr( $key, $attr->{$key});
		      }
		  }
		
		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/Qualifier elements") if $bsml_logger->is_debug;

		  foreach my $BsmlQualifier ($BsmlFeature->children( 'Qualifier' ))
		  {
		      my $attr = $BsmlQualifier->atts();
		      $feat->addBsmlQualifier( $attr->{'value-type'} , $attr->{'value'} ); 
		  }
		  
		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/Interval-loc elements") if $bsml_logger->is_debug;

		  foreach my $BsmlIntervalLoc ($BsmlFeature->children( 'Interval-loc' ))
		  {
		      my $attr = $BsmlIntervalLoc->atts();
		      if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };

		      $feat->addBsmlIntervalLoc( $attr->{'startpos'} , $attr->{'endpos'}, $attr->{'complement'} ); 
		  }

		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/Site-loc elements") if $bsml_logger->is_debug;

		  foreach my $BsmlSiteLoc ($BsmlFeature->children( 'Site-loc' ))
		  {
		      my $attr = $BsmlSiteLoc->atts();
		      
		      if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };
		      $feat->addBsmlSiteLoc( $attr->{'sitepos'} , $attr->{'complement'}, $attr->{'class'} ); 
		  }

		  # add Feature level Bsml Attribute elements 

		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/Attribute elements") if $bsml_logger->is_debug;

		  foreach my $BsmlAttr ( $BsmlFeature->children( 'Attribute' ) )
		  {
		      my $attr = $BsmlAttr->atts();
		      $feat->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
		  }

		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/Link elements") if $bsml_logger->is_debug;

		  foreach my $BsmlLink ( $BsmlFeature->children( 'Link' ) )
		  {
		      my $attr = $BsmlLink->atts();
		      $feat->addBsmlLink( $attr->{'rel'}, $attr->{'href'} );
		  }

		  $bsml_logger->debug("Attempting to process //Feature-table/Feature/Attribute-list elements") if $bsml_logger->is_debug;

		  foreach my $attlist( $BsmlFeature->children( 'Attribute-list' )){
		     
		      my $listref = [];

		      foreach my $BsmlAttr ( $attlist->children( 'Attribute' ) ) {
			  
			  my $attr = $BsmlAttr->atts();
	    
			  push (@{$listref}, $attr);
		      }

		      $feat->addBsmlAttributeList( $listref );
		  }   

		      $bsml_logger->debug("called feat->addBsmlAttributeList") if $bsml_logger->is_debug;
		      
	      }
	  }

	  $bsml_logger->debug("Attempting to process all //Feature-table/Feature-group elements") if $bsml_logger->is_debug;

	  foreach my $BsmlFGroup  ($BsmlFTables->children( 'Feature-group' )) 
	  {
	    my $group = $bsmlseq->{'BsmlFeatureGroups'}[$bsmlseq->addBsmlFeatureGroup()];
	    
	    my $attr = $BsmlFGroup->atts();
	    
	    foreach my $key ( keys( %{$attr} ) )
	      {
		$group->addattr( $key, $attr->{$key} );
	      }

	    #if an id has been specified, add the table to the general object lookups
	    if( my $id = $group->returnattr('id')){
	      BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $group );}

	    if( $BsmlFGroup->text() ){
	      $group->setText( $BsmlFGroup->text() ); 
	    }

	    foreach my $BsmlLink ( $BsmlFGroup->children( 'Link' ) )
	      {
		my $attr = $BsmlLink->atts();
		$group->addBsmlLink( $attr->{'rel'}, $attr->{'href'} );
	      }

	    foreach my $BsmlFGroupMember ( $BsmlFGroup->children( 'Feature-group-member' ))
	      {
		my $attr = $BsmlFGroupMember->atts();
		my $text = $BsmlFGroupMember->text();

		$group->addBsmlFeatureGroupMember( $attr->{'featref'}, $attr->{'feature-type'}, $attr->{'group-type'}, $text );
	      }

	    #if the feature group is part of a group-set, put it into the lookup tables
	    #this is the basis for returning all the transcripts (feature-groups) associated 
	    #with a gene

	    if( my $groupset = $group->returnattr('group-set'))
		{
		  BSML::BsmlDoc::BsmlSetFeatureGroupLookup( $groupset, $group );
		}
	    
	  }
	
      }

    $bsml_logger->debug("Attempting to process all //Sequence/Cross-reference/ elements") if $bsml_logger->is_debug;

    foreach my $crossRef( $seq->children( 'Cross-reference' )) {
	
	my $attr = $crossRef->atts();
	
	my $xref = $bsmlseq->returnBsmlCrossReferenceR( $bsmlseq->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) ) {
	    $xref->addattr( $key, $attr->{$key});
        }
    }

    $bsml_logger->debug("Attempting to process //Sequence/Attribute-list elements") if $bsml_logger->is_debug;
    
    foreach my $attlist( $seq->children( 'Attribute-list' )){
	
	my $listref = [];

	foreach my $BsmlAttr ( $attlist->children( 'Attribute' ) ) {

	    my $attr = $BsmlAttr->atts();
	    
	    push (@{$listref}, $attr);
	}
	
	$bsmlseq->addBsmlAttributeList( $listref );
    }

    $bsml_logger->debug("Attempting to process //Sequence/Link elements") if $bsml_logger->is_debug;

    foreach my $BsmlLink ( $seq->children( 'Link' ) )
    {
	my $attr = $BsmlLink->atts();
	$bsmlseq->addBsmlLink( $attr->{'rel'}, $attr->{'href'} );
    }
 
    $twig->purge_up_to( $seq );
  }

sub seqPairAlignmentHandler
  {
     my ($twig, $seq_aln) = @_;

     # add a new BsmlSeqPairAlignment object to the bsmlDoc

     my $bsmlaln = $bsmlDoc->{'BsmlSeqPairAlignments'}[$bsmlDoc->addBsmlSeqPairAlignment()];
     
     # add the BsmlSeqPairAlignment element's attributes

     my $attr = $seq_aln->atts();

     foreach my $key ( keys( %{$attr} ) )
       {
	 $bsmlaln->addattr( $key, $attr->{$key} );
       }

     #if an id has been specified, add the alignment pair to the general object lookups
     if( my $id = $bsmlaln->returnattr('id')){
       BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $bsmlaln );}


     #if refseq and compseq are defined, add the alignment to the alignment lookups
     if( (my $refseq = $bsmlaln->returnattr('refseq')) && (my $compseq = $bsmlaln->returnattr('compseq')))
       {
	 BSML::BsmlDoc::BsmlSetAlignmentLookup( $refseq, $compseq, $bsmlaln );
       }
   
     # add Bsml Attribute elements to the BsmlSeqPairAlignment 

     foreach my $BsmlAttr ( $seq_aln->children( 'Attribute' ) )
      {
	my $attr = $BsmlAttr->atts();
	$bsmlaln->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
      }

     foreach my $BsmlLink ( $seq_aln->children( 'Link' ) )
     {
	 my $attr = $BsmlLink->atts();
	 $bsmlaln->addBsmlLink( $attr->{'rel'}, $attr->{'href'} );
     }
     
     foreach my $seq_run ( $seq_aln->children('Seq-pair-run') )
       {
	 my $bsmlseqrun = $bsmlaln->returnBsmlSeqPairRunR( $bsmlaln->addBsmlSeqPairRun() ); 

	 my $attr = $seq_run->atts();
	 foreach my $key ( keys( %{$attr} ) ){
	   $bsmlseqrun->addattr( $key, $attr->{$key} );
	 }

	 foreach my $BsmlAttr ( $seq_run->children( 'Attribute' ) ){
	   my $attr = $BsmlAttr->atts();
	   $bsmlseqrun->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
	 }

	 foreach my $BsmlLink ( $seq_run->children( 'Link' ) )
	   {
	     my $attr = $BsmlLink->atts();
	     $bsmlseqrun->addBsmlLink( $attr->{'rel'}, $attr->{'href'} );
	   }
       }  
     
     foreach my $crossRef( $seq_aln->children( 'Cross-reference' ))
     {

	my $attr = $crossRef->atts();

	my $xref = $bsmlaln->returnBsmlCrossReferenceR( $bsmlaln->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) )
        {
	    $xref->addattr( $key, $attr->{$key});
        }
    }

     # Purge the twig rooted at the SeqPairAlignment Element
     $twig->purge_up_to( $seq_aln );
  }

sub analysisHandler
  {
    my ($twig, $analysis) = @_;
    my $bsml_analysis = $bsmlDoc->{'BsmlAnalyses'}[$bsmlDoc->addBsmlAnalysis()];
    my $attr = $analysis->atts();
    
    foreach my $key ( keys( %{$attr} ) )
      {
	$bsml_analysis->addattr( $key, $attr->{$key} );
      }

    #if an id has been specified, add the analysis to the general object lookups
    if( my $id = $bsml_analysis->returnattr('id')){
      BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $bsml_analysis );}

    foreach my $BsmlAttr ( $analysis->children( 'Attribute' ) )
      {
	my $attr = $BsmlAttr->atts();
	$bsml_analysis->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
      }

    foreach my $BsmlLink ( $analysis->children( 'Link' ) )
      {
	my $attr = $BsmlLink->atts();
	$bsml_analysis->addBsmlLink( $attr->{'rel'}, $attr->{'href'} );
      }

    foreach my $crossRef( $analysis->children( 'Cross-reference' )){
	
	my $attr = $crossRef->atts();
	
	my $xref = $bsml_analysis->returnBsmlCrossReferenceR( $bsml_analysis->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) ) {
	    $xref->addattr( $key, $attr->{$key});
	}
    }

  }

sub multiAlignmentTableHandler
{
    my ($twig, $mTable) = @_;
    my $bsml_mTable = $bsmlDoc->returnBsmlMultipleAlignmentTableR( $bsmlDoc->addBsmlMultipleAlignmentTable() );
  
    addBsmlAttrLinks( $twig, $mTable, $bsml_mTable );

    foreach my $alnSum ( $mTable->children( 'Alignment-summary' ) )
    {
	alignmentSummaryHandler( $twig, $alnSum, $bsml_mTable );
    }

    foreach my $pAlns ( $mTable->children( 'Pairwise-alignments' ))
    {
	pairwiseAlignmentsHandler( $twig, $pAlns, $bsml_mTable );
    }

    foreach my $alnSeq ( $mTable->children( 'Sequence-alignment' ))
    {
	sequenceAlignmentHandler( $twig, $alnSeq, $bsml_mTable );
    }

    foreach my $crossRef( $mTable->children( 'Cross-reference' ))
    {

	my $attr = $crossRef->atts();

	my $xref = $bsml_mTable->returnBsmlCrossReferenceR( $bsml_mTable->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) )
        {
	    $xref->addattr( $key, $attr->{$key});
        }
    }     

    # Purge the twig rooted at the multi-alignment table
    $twig->purge_up_to( $mTable );
}

sub alignmentSummaryHandler
{
    my ($twig, $alnSum, $bsmlMTable ) = @_;
    my $alignmentSummary = $bsmlMTable->returnBsmlAlignmentSummaryR( $bsmlMTable->addBsmlAlignmentSummary() );

    addBsmlAttrLinks( $twig, $alnSum, $alignmentSummary );

    foreach my $alnSeq ( $alnSum->children( 'Aligned-sequence' ) )
    {
	alignedSequenceHandler( $twig, $alnSeq, $alignmentSummary );
    }    
}

sub alignedSequenceHandler
{
    my ($twig, $alnSeq, $bsmlAlnSum ) = @_;
    my $alignedSequence = $bsmlAlnSum->returnBsmlAlignedSequenceR( $bsmlAlnSum->addBsmlAlignedSequence() );

    addBsmlAttrLinks( $twig, $alnSeq, $alignedSequence );       
}

sub pairwiseAlignmentsHandler
{
    my ($twig, $pAlns, $bsmlMTable ) = @_;

    my $bsmlpAlns = $bsmlMTable->returnBsmlPairwiseAlignmentsR( $bsmlMTable->addBsmlPairwiseAlignments() );
    addBsmlAttrLinks( $twig, $pAlns, $bsmlpAlns );       

    foreach my $alnPair ( $pAlns->children( 'Aligned-pair' ) )
    {
	alignedPairHandler( $twig, $alnPair, $bsmlpAlns );
    }    
}

sub alignedPairHandler
{
    my ($twig, $alnP, $bsmlpAlns ) = @_;

    my $bsmlAlnP = $bsmlpAlns->returnBsmlAlignedPairR( $bsmlpAlns->addBsmlAlignedPair() );
    addBsmlAttrLinks( $twig, $alnP, $bsmlAlnP );  
}

sub sequenceAlignmentHandler
{
    my ($twig, $seqAln, $bsmlMTable) = @_;
    my $bsmlSeqAln = $bsmlMTable->returnBsmlSequenceAlignmentR( $bsmlMTable->addBsmlSequenceAlignment() );

    addBsmlAttrLinks( $twig, $seqAln, $bsmlSeqAln );

    foreach my $seqDat ( $seqAln->children( 'Sequence-data' ))
    {
	sequenceDataHandler( $twig, $seqDat, $bsmlSeqAln );
    }

    foreach my $alnCon ( $seqAln->children( 'Alignment-consensus' ))
    {
	alignmentConsensusHandler( $twig, $alnCon, $bsmlSeqAln );
    }
}

sub sequenceDataHandler
{
    my( $twig, $seqDat, $bsmlSeqAln ) = @_;
    my $bsmlSeqDat = $bsmlSeqAln->returnBsmlSequenceDataR( $bsmlSeqAln->addBsmlSequenceData() );

    addBsmlAttrLinks( $twig, $seqDat, $bsmlSeqDat );

    $bsmlSeqDat->addSequenceAlignmentData( $seqDat->text() );
}

sub alignmentConsensusHandler
{

}

######################################################################
##
#  Support foruse of data associated with Genome elements and
#  children.


sub genomeHandler
{
    my ($twig, $genome) = @_;
    my $bsmlGenome = $bsmlDoc->returnBsmlGenomeR( $bsmlDoc->addBsmlGenome() );

    addBsmlAttrLinks( $twig, $genome, $bsmlGenome );

    foreach my $organism ( $genome->children( 'Organism' ))
    {


	organismHandler( $twig, $organism, $bsmlGenome );
    }

    foreach my $crossRef( $genome->children( 'Cross-reference' ))
    {
	
	my $attr = $crossRef->atts();
	
	my $xref = $bsmlGenome->returnBsmlCrossReferenceR( $bsmlGenome->addBsmlCrossReference() );

	foreach my $key ( keys( %{$attr} ) )
        {
	    $xref->addattr( $key, $attr->{$key});
        }
    }

    $twig->purge_up_to( $genome );
}

sub organismHandler
{
    my ($twig, $organism, $bsmlGenome) = @_;

    my $bsmlOrganism = $bsmlGenome->returnBsmlOrganismR( $bsmlGenome->addBsmlOrganism() );
    
    addBsmlAttrLinks( $twig, $organism, $bsmlOrganism );

    foreach my $strain ($organism->children( 'Strain' ))
    {
	strainHandler( $twig, $strain, $bsmlOrganism );
    }   
}

sub strainHandler
{
    my ($twig, $strain, $bsmlOrganism ) = @_;

    my $bsmlStrain = $bsmlOrganism->returnBsmlStrainR( $bsmlOrganism->addBsmlStrain() );
   
    addBsmlAttrLinks( $twig, $strain, $bsmlStrain );
}

# Generic function to add attributes, Bsml Attributes, and 
# Bsml links to generic objects

sub addBsmlAttrLinks
{
    my ($twig, $elem, $bsmlObj ) = @_;

    my $attr = $elem->atts();
    
    foreach my $key ( keys( %{$attr} ) )
    {
	$bsmlObj->addattr( $key, $attr->{$key} );
    }

    foreach my $BsmlAttr ( $elem->children( 'Attribute' ) )
    {
	my $attr = $BsmlAttr->atts();
	$bsmlObj->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
    }

    foreach my $BsmlLink ( $elem->children( 'Link' ) )
    {
	my $attr = $BsmlLink->atts();
	$bsmlObj->addBsmlLink( $attr->{'rel'}, $attr->{'href'} );
    }    
}

1
