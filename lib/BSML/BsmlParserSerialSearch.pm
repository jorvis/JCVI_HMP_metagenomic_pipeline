package BSML::BsmlParserSerialSearch;

# Directed parser for serial output of BSML search encodings. 

use strict;
use warnings;
use XML::Twig;
BEGIN {
use BSML::BsmlFeature;
use BSML::BsmlSequence;
use BSML::BsmlAnalysis;
use BSML::BsmlMultipleAlignmentTable;
use BSML::BsmlGenome;
use BSML::BsmlSeqPairAlignment;
use BSML::Logger;
}

use Data::Dumper;

my $logger = BSML::Logger::get_logger("Logger::BSML");

#-------------------------------------------------------------------------------------------------------------------
#
#
# Clients must define callback routines which are called during parsing. Callbacks may be defined
# for Sequence, SeqPairAlignment, and Analysis objects. These functions are called with an object 
# reference during the parse as each element is encountered. Twigs are rooted at Sequence, SeqPairAlignment, Feature, or Analysis elements
# for maximum efficiency. 
#
# modification: Jay Sundaram 2003-12-07
#               1) introduction of log4perl logging
#               2) MultipleAlignmentCallBack function
#
#
#
#
#
#
#
#-------------------------------------------------------------------------------------------------------------------



sub new
  {
    my $class = shift;
    my (%args) = @_; 
    my $self = {};
    bless $self, $class;

    $self->{'Ignore'} = {};

    $logger->debug("Instantiating BsmlParserSerialSearch") if($logger->is_debug());

    if( $args{'AlignmentCallBack'} ){
	$self->{'AlignmentCallBack'} = $args{'AlignmentCallBack'};
	
	$self->{'Roots'}->{'Seq-pair-alignment'} = sub {my ($twig, $seq_aln ) = @_;
							$logger->debug("Starting alignment call back") if($logger->is_debug());
							my $bsml_aln = seqPairAlignmentHandler($twig,$seq_aln);
							$self->{'AlignmentCallBack'}( $bsml_aln );
						    };
    }

    if( $args{'MultipleAlignmentCallBack'} ){
		
	$self->{'MultipleAlignmentCallBack'} = $args{'MultipleAlignmentCallBack'};

	$self->{'Roots'}->{'Multiple-alignment-table'} = sub {my ($twig, $seq_aln ) = @_;
							      $logger->debug("Starting multiple alignment call back") if($logger->is_debug());
							      my $bsml_aln = multipleAlignmentHandler($twig,$seq_aln);
							      $self->{'MultipleAlignmentCallBack'}( $bsml_aln );
							  };
    }


    if( $args{'AnalysisCallBack'} ){
	$self->{'AnalysisCallBack'} = $args{'AnalysisCallBack'};
	
	$self->{'Roots'}->{'Analysis'} = sub {my ($twig, $analysis) = @_;
					      $logger->debug("Starting analysis call back") if($logger->is_debug());
					      my $bsml_analysis = analysisHandler( $twig, $analysis );
					      $self->{'AnalysisCallBack'}( $bsml_analysis );
					  };
    }

    if( $args{'SequenceCallBack'} ){
	$self->{'SequenceCallBack'} = $args{'SequenceCallBack'};
	if( defined( $args{'ReadFeatureTables'} ) ){
	    if( $args{'ReadFeatureTables'} == 1 ){
		$self->{'Roots'}->{'Sequence'} = sub {my ($twig, $sequence) = @_;
						      $logger->debug("Starting sequence handler and reading feature tables") if($logger->is_debug());
						      my $bsml_sequence = sequenceHandler( $twig, $sequence );
						      $self->{'SequenceCallBack'}( $bsml_sequence );
					  };
		    }
	    else{
			$self->{'Roots'}->{'Sequence'} = sub {my ($twig, $sequence) = @_;
							      $logger->debug("Starting sequence handler and skipping reading feature tables") if($logger->is_debug());
							      if( defined($args{'SequenceClassFilter'}) && $args{'SequenceClassFilter'} ne "" && $args{'SequenceClassFilter'} ne $sequence->atts()->{'class'}){
								  $logger->debug("Skipping parse of sequence with class ".$sequence->atts()->{'class'}) if($logger->is_debug());
							      }
							      else{
								  my $bsml_sequence = minsequenceHandler( $twig, $sequence );
								  $self->{'SequenceCallBack'}( $bsml_sequence );
							      }
					  }; 
			$self->{'Ignore'}->{'Feature-tables'}=1;
		    }
	}
	else{
	    $self->{'Roots'}->{'Sequence'} = sub {my ($twig, $sequence) = @_;
						  $logger->debug("Reading sequences and subfeatures") if($logger->is_debug());
						  my $bsml_sequence = sequenceHandler( $twig, $sequence );
						  $self->{'SequenceCallBack'}( $bsml_sequence );
					  };
	}
    }

    if( $args{'FeatureCallBack'} ){
	$self->{'FeatureCallBack'} = $args{'FeatureCallBack'};
	
	$self->{'Roots'}->{'Feature'} = sub {my ($twig, $feature) = @_;
					     $logger->debug("Starting feature call back") if($logger->is_debug());
					     my $bsml_feature = featureHandler( $twig, $feature );
					     $self->{'FeatureCallBack'}( $bsml_feature );
					 };
    }

    if( $args{'SeqFeatureCallBack'} ){
	$self->{'SeqFeatureCallBack'} = $args{'SeqFeatureCallBack'};
	
	$self->{'Roots'}->{'Feature'} = sub {my ($twig, $feature) = @_;
					     $logger->debug("Starting feature call back") if($logger->is_debug());
					     my $bsml_feature = SeqFeatureHandler( $twig, $feature );
					     $self->{'SeqFeatureCallBack'}( $bsml_feature );
					 };
    }

    if( $args{'GenomeCallBack'} ){
	$self->{'GenomeCallBack'} = $args{'GenomeCallBack'};
	
	$self->{'Roots'}->{'Genome'} = sub {my ($twig, $genome ) = @_;
					    $logger->debug("Starting genome call back") if($logger->is_debug());
					    my $bsml_genome = genomeHandler($twig,$genome);
					    $self->{'GenomeCallBack'}( $bsml_genome );
					};
    }
	
    return $self;
  }

sub parse  {

    $logger->debug("") if($logger->is_debug());

    my $self = shift;
    my ( $fileOrHandle ) = @_;

    my $bsml_logger = BSML::Logger::get_logger( "Bsml" );

    # Set "rooted" twig handlers on Sequence, SeqPairAlignment, and Analysis elements. 
    my $twig = new XML::Twig( TwigRoots => $self->{'Roots'},
			      ignore_elts => $self->{'Ignore'});

    # parse will die if an xml syntax error is encountered or if
    # there is an io problem

    if( $fileOrHandle ) {
		if (ref($fileOrHandle) && ($fileOrHandle->isa("IO::Handle") || $fileOrHandle->isa("GLOB"))) {

			$logger->debug("twig->parse") if($logger->is_debug());
			$twig->parse( $fileOrHandle );
		} 
		else {

			# This is not a file handle- is an actual file.

			if ($fileOrHandle =~ /\.bsml\.(gz|gzip)$/) {

				# This file has been zipped

				my ($fh);
				
				open ($fh, "<:gzip", $fileOrHandle) || $logger->logdie("Could not open zipped file '$fileOrHandle': $!");

				$logger->debug("twig->parsefile $fileOrHandle") if($logger->is_debug());

				$twig->parse( $fh );

			} 
			else {

				# Regular unzipped file

				$twig->parsefile( $fileOrHandle );

			}

		}
    } 
    else {
		$twig->parse( \*STDIN );
    }
    $logger->debug("twig->dispose") if($logger->is_debug());
    $twig->dispose();  
}

# handler for SeqPairAlignment elements. Returns a reference to a BSML::SeqPairAlignment object.

sub seqPairAlignmentHandler
  {
     my ($twig, $seq_aln ) = @_;

     # create a new BsmlSeqPairAlignment object

     my $bsmlaln = new BSML::BsmlSeqPairAlignment;
     
     # add the BsmlSeqPairAlignment element's attributes

     my $attr = $seq_aln->atts();

     foreach my $key ( keys( %{$attr} ) )
       {
	 $bsmlaln->addattr( $key, $attr->{$key} );
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
	 $bsmlaln->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
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
	     $bsmlseqrun->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
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
     $twig->purge();

     return $bsmlaln;

  }

# handler for analysis elements. Returns a reference to a BSML::Analysis object.

sub analysisHandler
  {
    my ($twig, $analysis) = @_;
    my $bsml_analysis = new BSML::BsmlAnalysis;
    my $attr = $analysis->atts();
    
    foreach my $key ( keys( %{$attr} ) )
      {
	$bsml_analysis->addattr( $key, $attr->{$key} );
      }

    foreach my $BsmlAttr ( $analysis->children( 'Attribute' ) )
      {
	my $attr = $BsmlAttr->atts();
	$bsml_analysis->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
      }

    foreach my $BsmlLink ( $analysis->children( 'Link' ) )
      {
	my $attr = $BsmlLink->atts();
	$bsml_analysis->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
      }


    foreach my $crossRef( $analysis->children( 'Cross-reference' )){
	
	my $attr = $crossRef->atts();
	
	my $xref = $bsml_analysis->returnBsmlCrossReferenceR( $bsml_analysis->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) ) {
	    $xref->addattr( $key, $attr->{$key});
	}
    }



    $twig->purge();
    return $bsml_analysis;
  }

# handler for sequence elements. Returns a reference to a BSML::Sequence object.

sub sequenceHandler
  {
    my ($twig, $seq) = @_;


    $logger->debug("Entered sequenceHandler") if $logger->is_debug;

    # add a new Sequence object to the bsmlDoc

    my $bsmlseq = new BSML::BsmlSequence;
    
    # add the sequence element's attributes

    my $attr = $seq->atts();


    $logger->debug("Processing all //Sequence/\@attributes") if $logger->is_debug;

    foreach my $key ( keys( %{$attr} ) ) {
	$bsmlseq->addattr( $key, $attr->{$key} );
    }

    # add Sequence level Bsml Attribute elements 


    $logger->debug("Processing all //Sequence/Attribute elements") if $logger->is_debug;


    foreach my $BsmlAttr ( $seq->children( 'Attribute' ) ) {
	my $attr = $BsmlAttr->atts();
	$bsmlseq->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
    }
    
    # add any Bsml Links
    
    $logger->debug("Processing all //Sequence/Link elements") if $logger->is_debug;

    foreach my $BsmlLink ( $seq->children( 'Link' ) )
      {

	my $attr = $BsmlLink->atts();

	#
	# editor:   sundaram@tigr.org
	# date:     2005-08-18
	# bgzcase:  2053
	# URL:      http://serval.tigr.org:8080/bugzilla/show_bug.cgi?id=2053
	# comment:  The bsml2chado.pl script retrieves BSML <Link> info, however pls note that it is
	#           only interested in the following attributes: //Link/@rel and //Link/@href
	#
	$bsmlseq->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );

      }

    # add raw sequence data if found 


    $logger->debug("Processing all //Sequence/Seq-data elements") if $logger->is_debug;


    my $seqDat = $seq->first_child( 'Seq-data' );
    
    if( $seqDat ){
      $bsmlseq->addBsmlSeqData( $seqDat->text() );
    }

    # add appended sequence data if found


    $logger->debug("Processing all //Sequence/Seq-data-import elements") if $logger->is_debug;

    my $seqDatImport = $seq->first_child( 'Seq-data-import' );
    


    $logger->debug("Processing all //Sequence/Seq-data-import/\@attributes") if $logger->is_debug;


    if( $seqDatImport )       {
	my $attr = $seqDatImport->atts();
	$bsmlseq->addBsmlSeqDataImport( $attr->{'format'}, $attr->{'source'}, $attr->{'id'},$attr->{'identifier'});
    }
    
    # add numbering information


    $logger->debug("Processing all //Sequence/Numbering/ elements") if $logger->is_debug;

    my $numbering = $seq->first_child( 'Numbering' );

    if( $numbering ) {
	
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


	$logger->debug("Processing //Sequence/Numbering/Attribute elements") if $logger->is_debug;

	foreach my $BsmlAttr ( $numbering->children( 'Attribute' ) ) {
	    my $attr = $BsmlAttr->atts();
	    $bsmlNumbering->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
	}
    }

    # add Feature Tables with Feature, Reference, and Feature-group Elements




    my $BsmlFTables = $seq->first_child( 'Feature-tables' );

    if( $BsmlFTables ) {

	$logger->debug("Processing first_child //Feature-table/ element") if $logger->is_debug;
	
	foreach my $BsmlFTable ( $BsmlFTables->children( 'Feature-table' ) ) {


#	    $logger->debug("BsmlTable: " . Dumper $BsmlFTable) if $logger->is_debug;

	    
	    my $table = $bsmlseq->{'BsmlFeatureTables'}[$bsmlseq->addBsmlFeatureTable()];
	    
	    my $attr = $BsmlFTable->atts();
	    
	    
	    $logger->debug("Processing //Feature-table/\@attibutes") if $logger->is_debug;
	    
	    foreach my $key ( keys( %{$attr} ) ){
		$table->addattr( $key, $attr->{$key} );
	    }
	    

	    $logger->debug("Processing //Feature-table/Link elements") if $logger->is_debug;


	    foreach my $BsmlLink ( $BsmlFTable->children( 'Link' ) ) {
		my $attr = $BsmlLink->atts();
		$table->addBsmlLink( $attr->{'title'}, $attr->{'href'}, $attr->{'role'} );
	    }


	    $logger->debug("Processing //Feature-table/Reference/ elements") if $logger->is_debug;

	    foreach my $BsmlRef ($BsmlFTable->children( 'Reference' ) ) {
		
		my $ref = $table->{'BsmlReferences'}[$table->addBsmlReference()];
		my $attr = $BsmlRef->atts();
		
		foreach my $key ( keys( %{$attr} ) ) {
		    $ref->addattr( $key, $attr->{$key} );
		}
		
		foreach my $BsmlAuthor ($BsmlRef->children( 'RefAuthors' )) {
		    $ref->addBsmlRefAuthors( $BsmlAuthor->text() );
		}
		
		foreach my $BsmlRefTitle ($BsmlRef->children( 'RefTitle' )) {
		    $ref->addBsmlRefTitle( $BsmlRefTitle->text() );
		}
		
		foreach my $BsmlRefJournal ($BsmlRef->children( 'RefJournal' )) {
		    $ref->addBsmlRefJournal( $BsmlRefJournal->text() );
		}

		foreach my $BsmlLink ( $BsmlRef->children( 'Link' ) ) {
		    my $attr = $BsmlLink->atts();
		    $ref->addBsmlLink( $attr->{'title'}, $attr->{'href'}, $attr->{'role'} );
		}
	    }
	    

	    $logger->debug("Attempting to process //Feature-table/Feature/ elements") if $logger->is_debug;

	    foreach my $BsmlFeature ($BsmlFTable->children( 'Feature' )){
		
		my $feat = $table->{'BsmlFeatures'}[$table->addBsmlFeature()];
		my $attr = $BsmlFeature->atts();
		

		$logger->debug("Attempting to process //Feature-table/Feature/\@attributes") if $logger->is_debug;


		foreach my $key ( keys( %{$attr} ) ) {
		    $feat->addattr( $key, $attr->{$key} );
		}
		

		$logger->debug("Attempting to process //Feature-table/Feature/Cross-reference elements") if $logger->is_debug;

		foreach my $crossRef( $BsmlFeature->children( 'Cross-reference' )) {
		    
		    my $attr = $crossRef->atts();
		    
		    my $xref = $feat->returnBsmlCrossReferenceR( $feat->addBsmlCrossReference );
		    
		    foreach my $key ( keys( %{$attr} ) ) {

			$xref->addattr( $key, $attr->{$key});
		    }
		}
		
		# add Feature level Bsml Attribute elements 
		
		$logger->debug("Attempting to process //Feature-table/Feature/Attribute elements") if $logger->is_debug;

		foreach my $BsmlAttr ( $BsmlFeature->children( 'Attribute' ) ) {
		    my $attr = $BsmlAttr->atts();
		    $feat->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
		}
		
		$logger->debug("Attempting to process //Feature-table/Feature/Qualifier elements") if $logger->is_debug;

		foreach my $BsmlQualifier ($BsmlFeature->children( 'Qualifier' )) {
		    my $attr = $BsmlQualifier->atts();
		    $feat->addBsmlQualifier( $attr->{'value-type'} , $attr->{'value'} ); 
		}


		$logger->debug("Attempting to process //Feature-table/Feature/Interval-loc elements") if $logger->is_debug;

		foreach my $BsmlIntervalLoc ($BsmlFeature->children( 'Interval-loc' )) {

		    my $attr = $BsmlIntervalLoc->atts();
		    if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };
		    
		    $feat->addBsmlIntervalLoc( $attr->{'startpos'} , $attr->{'endpos'}, $attr->{'complement'} ); 
		}


		$logger->debug("Attempting to process //Feature-table/Feature/Site-loc elements") if $logger->is_debug;

		foreach my $BsmlSiteLoc ($BsmlFeature->children( 'Site-loc' )) {
		    my $attr = $BsmlSiteLoc->atts();
		    
		    if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };
		    $feat->addBsmlSiteLoc( $attr->{'sitepos'} , $attr->{'complement'}, $attr->{'class'} ); 
		}



		$logger->debug("Attempting to process //Feature-table/Feature/Link elements") if $logger->is_debug;

		foreach my $BsmlLink ( $BsmlFeature->children( 'Link' ) ) {
		    my $attr = $BsmlLink->atts();
		    $feat->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
		}


		$logger->debug("Attempting to process //Feature-table/Feature/Attribute-list elements") if $logger->is_debug;

		foreach my $attlist( $BsmlFeature->children( 'Attribute-list' )){
		    
		    my $listref = [];

		    foreach my $BsmlAttr ( $attlist->children( 'Attribute' ) ) {

			my $attr = $BsmlAttr->atts();
			
			push (@{$listref}, $attr);
		    }
 
		    $feat->addBsmlAttributeList( $listref );

		    $logger->debug("called feat->addBsmlAttributeList") if $logger->is_debug;

		}
	    }
	}


	$logger->debug("Attempting to process all //Feature-table/Feature-group elements") if $logger->is_debug;

	foreach my $BsmlFGroup  ($BsmlFTables->children( 'Feature-group' )) {

	    my $group = $bsmlseq->{'BsmlFeatureGroups'}[$bsmlseq->addBsmlFeatureGroup()];
	    
	    my $attr = $BsmlFGroup->atts();
	    
	    foreach my $key ( keys( %{$attr} ) )
	      {
		$group->addattr( $key, $attr->{$key} );
	      }

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
	    
	  }
	
      }


    
    $logger->debug("Attempting to process all //Sequence/Cross-reference/ elements") if $logger->is_debug;

    foreach my $crossRef( $seq->children( 'Cross-reference' )) {
	
	my $attr = $crossRef->atts();
	
	my $xref = $bsmlseq->returnBsmlCrossReferenceR( $bsmlseq->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) ) {
	    $xref->addattr( $key, $attr->{$key});
        }
    }
    
    
    $logger->debug("Attempting to process //Sequence/Attribute-list elements") if $logger->is_debug;
    
    foreach my $attlist( $seq->children( 'Attribute-list' )){
	
	my $listref = [];

	foreach my $BsmlAttr ( $attlist->children( 'Attribute' ) ) {

	    my $attr = $BsmlAttr->atts();
	    
	    push (@{$listref}, $attr);
	}
	
	$bsmlseq->addBsmlAttributeList( $listref );
    }
  
    $twig->purge;
    return $bsmlseq;
  }

# Feature tables and Feature groups are not parsed with the minisequenceHandler

sub minsequenceHandler
  {

    my ($twig, $seq, $detect) = @_;

    # add a new Sequence object to the bsmlDoc

    my $bsmlseq = new BSML::BsmlSequence;
    
    # add the sequence element's attributes

    my $attr = $seq->atts();

    foreach my $key ( keys( %{$attr} ) )
      {
	$bsmlseq->addattr( $key, $attr->{$key} );
      }

    # add Sequence level Bsml Attribute elements 

    foreach my $BsmlAttr ( $seq->children( 'Attribute' ) )
      {
	my $attr = $BsmlAttr->atts();
	$bsmlseq->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
      }

    # add any Bsml Links

    foreach my $BsmlLink ( $seq->children( 'Link' ) )
      {

	my $attr = $BsmlLink->atts();

	$bsmlseq->addBsmlLink( $attr->{'title'}, $attr->{'href'} );

	#
	# editor:   sundaram@tigr.org
	# date:     2005-08-18
	# bgzcase:  2053
	# URL:      http://serval.tigr.org:8080/bugzilla/show_bug.cgi?id=2053
	# comment:  The bsml2chado.pl script retrieves BSML <Link> info, however pls note that it is
	#           only interested in the following attributes: //Link/@rel and //Link/@href
	#
	$bsmlseq->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'});
	


      }

    # add raw sequence data if found 

    my $seqDat = $seq->first_child( 'Seq-data' );
    
    if( $seqDat ){
      $bsmlseq->addBsmlSeqData( $seqDat->text() );
    }

    # add appended sequence data if found

    my $seqDatImport = $seq->first_child( 'Seq-data-import' );
    
    if( $seqDatImport )
      {
	my $attr = $seqDatImport->atts();
	$bsmlseq->addBsmlSeqDataImport( $attr->{'format'}, $attr->{'source'}, $attr->{'id'},$attr->{'identifier'});
      }

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

	foreach my $BsmlAttr ( $numbering->children( 'Attribute' ) )
	{
	    my $attr = $BsmlAttr->atts();
	    $bsmlNumbering->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
	}
    }

      $logger->debug("Attempting to process all //Sequence/Cross-reference/ elements") if $logger->is_debug;

    foreach my $crossRef( $seq->children( 'Cross-reference' )) {
	
	my $attr = $crossRef->atts();
	
	my $xref = $bsmlseq->returnBsmlCrossReferenceR( $bsmlseq->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) ) {
	    $xref->addattr( $key, $attr->{$key});
        }
    }
    
    
    $logger->debug("Attempting to process //Sequence/Attribute-list elements") if $logger->is_debug;
    
    foreach my $attlist( $seq->children( 'Attribute-list' )){
	
	my $listref = [];

	foreach my $BsmlAttr ( $attlist->children( 'Attribute' ) ) {

	    my $attr = $BsmlAttr->atts();
	    
	    push (@{$listref}, $attr);
	}
	
	$bsmlseq->addBsmlAttributeList( $listref );
    }
    
    $twig->purge;
    return $bsmlseq;
}



# returns a BsmlFeature object

sub featureHandler
{
    my ($twig, $BsmlFeature) = @_;

    # print "executing featurHandler\n";

    my $bsmlseq = new BSML::BsmlSequence;
    my $feat = new BSML::BsmlFeature;

    my $attr = $BsmlFeature->atts();

    foreach my $key ( keys( %{$attr} ) )
    {
	$feat->addattr( $key, $attr->{$key} );
    }

    # add Feature level Bsml Attribute elements 

    foreach my $BsmlAttr ( $BsmlFeature->children( 'Attribute' ) )
      {
	my $attr = $BsmlAttr->atts();
	$feat->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
      }

    foreach my $BsmlQualifier ($BsmlFeature->children( 'Qualifier' ))
    {
	my $attr = $BsmlQualifier->atts();
	$feat->addBsmlQualifier( $attr->{'value-type'} , $attr->{'value'} ); 
    }

    foreach my $BsmlIntervalLoc ($BsmlFeature->children( 'Interval-loc' ))
    {
	my $attr = $BsmlIntervalLoc->atts();
	if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };
	
	$feat->addBsmlIntervalLoc( $attr->{'startpos'} , $attr->{'endpos'}, $attr->{'complement'} ); 
    }

    foreach my $BsmlSiteLoc ($BsmlFeature->children( 'Site-loc' ))
    {
	my $attr = $BsmlSiteLoc->atts();
	
	if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };
	$feat->addBsmlSiteLoc( $attr->{'sitepos'} , $attr->{'complement'}, $attr->{'class'} ); 
    }
    
    foreach my $BsmlLink ( $BsmlFeature->children( 'Link' ) )
    {
	my $attr = $BsmlLink->atts();
	$feat->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
    }


    foreach my $crossRef( $BsmlFeature->children( 'Cross-reference' )){
	
	my $attr = $crossRef->atts();
	
	my $xref = $feat->returnBsmlCrossReferenceR( $feat->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) ) {
	    $xref->addattr( $key, $attr->{$key});
	}
    }


    foreach my $attlist( $BsmlFeature->children( 'Attribute-list' )){
		     
	my $listref = [];
	
	foreach my $BsmlAttr ( $attlist->children( 'Attribute' ) ) {
	    
	    my $attr = $BsmlAttr->atts();
			  
	    push (@{$listref}, $attr);
	}

	$feat->addBsmlAttributeList( $listref );
    }   
   
    $twig->purge;

    return $feat;
}

# [$bsml_feature_ref, $bsml_sequence_ref] SeqFeatureHandler( $twig, $BsmlFeature )

# A handler that returns a reference to a list containg a bsml feature 
# reference and its parent sequence. This needs to be used in conjunction
# with the sequence handler in order to insure the parse tree containing 
# the parental elements is created. XML twig calls the feature callback each
# time the parser has fully parsed the subtree rooted at the feature element.
# This subtree is pruned (purged) from the Sequence parse tree after the call back
# has successfully completed. 

sub SeqFeatureHandler
{
    my ($twig, $BsmlFeature) = @_;

    # Retreive the parent sequence object - *NOTE this must be used with the sequence handler
    # or equivalent that looks at the parent nodes else (FeatureTable->FeatureTables->Sequence) 
    # will  be null.

    my $twigSequence = (($BsmlFeature->parent())->parent())->parent();
   

    ##################################################
    # Extract data from the sequence object, not looking at the feature tables. This is 
    # equivalent to the minSequenceHandler.

    my $bsmlseq = new BSML::BsmlSequence;

    # add the sequence element's attributes

    my $attr = $twigSequence->atts();

    foreach my $key ( keys( %{$attr} ) )
      {
	$bsmlseq->addattr( $key, $attr->{$key} );
    }

    # add Sequence level Bsml Attribute elements 

    foreach my $BsmlAttr ( $twigSequence->children( 'Attribute' ) )
    {
	my $attr = $BsmlAttr->atts();
	$bsmlseq->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
    }
    
    # add any Bsml Links
    
    foreach my $BsmlLink ( $twigSequence->children( 'Link' ) )
    {

	my $attr = $BsmlLink->atts();
	$bsmlseq->addBsmlLink( $attr->{'title'}, $attr->{'href'} );


	#
	# editor:   sundaram@tigr.org
	# date:     2005-08-18
	# bgzcase:  2053
	# URL:      http://serval.tigr.org:8080/bugzilla/show_bug.cgi?id=2053
	# comment:  The bsml2chado.pl script retrieves BSML <Link> info, however pls note that it is
	#           only interested in the following attributes: //Link/@rel and //Link/@href
	#
	$bsmlseq->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );

    }

    # add raw sequence data if found 
    
    my $seqDat = $twigSequence->first_child( 'Seq-data' );
    
    if( $seqDat ){
	$bsmlseq->addBsmlSeqData( $seqDat->text() );
    }

    # add appended sequence data if found
    
    my $seqDatImport = $twigSequence->first_child( 'Seq-data-import' );
    
    if( $seqDatImport )
      {
	  my $attr = $seqDatImport->atts();
	  $bsmlseq->addBsmlSeqDataImport( $attr->{'format'}, $attr->{'source'}, $attr->{'id'},$attr->{'identifier'});
      }
    
    my $numbering = $twigSequence->first_child( 'Numbering' );
    
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

	foreach my $BsmlAttr ( $numbering->children( 'Attribute' ) )
	{
	    my $attr = $BsmlAttr->atts();
	    $bsmlNumbering->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
	}
    }
    
    ###################################################
    #

    # parse the feature and construct a Bsml Feature object

    my $feat = new BSML::BsmlFeature;
    $attr = $BsmlFeature->atts();

    foreach my $key ( keys( %{$attr} ) )
    {
	$feat->addattr( $key, $attr->{$key} );
    }

    # add Feature level Bsml Attribute elements 

    foreach my $BsmlAttr ( $BsmlFeature->children( 'Attribute' ) )
      {
	my $attr = $BsmlAttr->atts();
	$feat->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
      }

    foreach my $BsmlQualifier ($BsmlFeature->children( 'Qualifier' ))
    {
	my $attr = $BsmlQualifier->atts();
	$feat->addBsmlQualifier( $attr->{'value-type'} , $attr->{'value'} ); 
    }

    foreach my $BsmlIntervalLoc ($BsmlFeature->children( 'Interval-loc' ))
    {
	my $attr = $BsmlIntervalLoc->atts();
	if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };
	
	$feat->addBsmlIntervalLoc( $attr->{'startpos'} , $attr->{'endpos'}, $attr->{'complement'} ); 
    }

    foreach my $BsmlSiteLoc ($BsmlFeature->children( 'Site-loc' ))
    {
	my $attr = $BsmlSiteLoc->atts();
	
	if( !( $attr->{'complement'} ) ){ $attr->{ 'complement' } = 0 };
	$feat->addBsmlSiteLoc( $attr->{'sitepos'} , $attr->{'complement'}, $attr->{'class'} ); 
    }
    
    foreach my $BsmlLink ( $BsmlFeature->children( 'Link' ) )
    {
	my $attr = $BsmlLink->atts();
	$feat->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
    }


    foreach my $crossRef( $BsmlFeature->children( 'Cross-reference' )){
	
	my $attr = $crossRef->atts();
	
	my $xref = $feat->returnBsmlCrossReferenceR( $feat->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) ) {
	    $xref->addattr( $key, $attr->{$key});
	}
    }


    foreach my $attlist( $BsmlFeature->children( 'Attribute-list' )){
		     
	my $listref = [];
	
	foreach my $BsmlAttr ( $attlist->children( 'Attribute' ) ) {
	    
	    my $attr = $BsmlAttr->atts();
			  
	    push (@{$listref}, $attr);
	}

	$feat->addBsmlAttributeList( $listref );
    }   
   

    # prune the subtree rooted at the feature element
    $twig->purge;


    # return a reference to a list containing references to the feature and parental
    # sequence objects

    my $listref = [$feat, $bsmlseq];
    return $listref;
}

















######################################################################
##
#  Support foruse of data associated with Genome elements and
#  children.


sub genomeHandler
{
    my ($twig, $genome) = @_;
    my $bsmlGenome = new BSML::BsmlGenome;


    $logger->debug("Entered genomeHandler") if $logger->is_debug;

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

	foreach my $BsmlAttr ( $crossRef->children( 'Attribute' ) ) {
	    my $attr = $BsmlAttr->atts();
	    $xref->addBsmlAttr( $attr->{'name'}, $attr->{'content'} );
	}


    }

    $twig->purge_up_to( $genome );


    $logger->debug("Leaving genomeHandler") if $logger->is_debug;

    return $bsmlGenome;
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
	$bsmlObj->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
    }    
}



#------------------------------------------------------------------------------
# multipleAlignmentHandler()
#
# Handler for multipleAlignment elements.
# Returns a reference to a BSML::SeqPairAlignment object.
#
#
#
# <Multiple-alignment-table molecule-type="protein">
#      <Alignment-summary seq-format="msf" seq-type="protein">
#           <Aligned-sequence length="1096" seqnum="1" name="lma2.204.m00029_protein:1"></Aligned-sequence>
#           <Aligned-sequence length="1096" seqnum="2" name="lma2.206.m00009_protein:2"></Aligned-sequence>
#      </Alignment-summary>
#      <Sequence-alignment sequences="1:2:">
#           <Sequence-data seq-name="lma2.204.m00029_protein:1">lma2.204.m00029_protein      MHSQRPCAGY NARCEHDRSV LPSLMSQQDE FFTGSSPRTA SVHSEARPHL 
#           </Sequence-data>
#           <Sequence-data seq-name="lma2.206.m00009_protein:2">lma2.206.m00009_protein      .MSLMPGYPV PRDQWVESAE ACAACSKRFT FFAFKENCPC CGRLFCSSCL 
#           </Sequence-data>
#           <Alignment-consensus></Alignment-consensus>
#      </Sequence-alignment>
# </Multiple-alignment-table>
#
#
#
#
#
#
#-------------------------------------------------------------------------------
sub multipleAlignmentHandler {

    $logger->debug("Entered multipleAlignmentHandler") if $logger->is_debug;


    my ($twig, $seq_aln ) = @_;

    $logger->logdie("twig was not defined")    if (!defined($twig));
    $logger->logdie("seq_aln was not defined") if (!defined($seq_aln));

    
    #
    # create a new BsmlMultipleAlignmentTable object
    #
    my $bsml_multi_aln = new BSML::BsmlMultipleAlignmentTable;
    

    #
    # This will retrieve all of the XML attributes associated to the  <Multiple-alignment-table> element
    #
    # e.g. <Multiple-alignment-table molecule-type="protein">
    # 
    my $attr = $seq_aln->atts();

    #
    # add all of the <Multiple-alignment-table> attributes to the BsmlMultipleAlignmentTable object
    #
    foreach my $key (keys %{$attr}){
	$bsml_multi_aln->addattr( $key, $attr->{$key} );
    }

    foreach my $BsmlLink ( $seq_aln->children( 'Link' ) ) {
	my $attr = $BsmlLink->atts();
	$bsml_multi_aln->addBsmlLink( $attr->{'rel'}, $attr->{'href'}, $attr->{'role'} );
    }
    
    ## extract values of BSML Attribute elements associated with BSML Multiple-alignment-table elements
    foreach my $BsmlAttr ( ($seq_aln->children('Attribute') )){
	my $attr  = $BsmlAttr->atts();
	$bsml_multi_aln->addBsmlAttr($attr->{'name'}, $attr->{'content'});
    }

    #
    # Process the <Alignment-summary> section:
    # e.g.
    #
    #    <Alignment-summary seq-format="msf" seq-type="protein">
    #         <Aligned-sequence length="1096" seqnum="1" name="lma2.204.m00029_protein:1"></Aligned-sequence>
    #         <Aligned-sequence length="1096" seqnum="2" name="lma2.206.m00009_protein:2"></Aligned-sequence>
    #    </Alignment-summary>
    #


    #
    # This will retrieve all of the <Alignment-summary> elements
    # 
    foreach my $alignment_summary ($seq_aln->children('Alignment-summary')){

	#
	# Creates and adds to BsmlAlignmentSummary object
	#
	my $bsml_alignment_summary = $bsml_multi_aln->returnBsmlAlignmentSummaryR( $bsml_multi_aln->addBsmlAlignmentSummary() ); 
	
	#
	# this will retrieve all of the XML attributes associated to the <Alignment-summary> element
	# e.g.  "seq-format" and "seq-type"
	my $attr = $alignment_summary->atts();


	#
	# this will add <Alignment-summary> XML attributes to the  BsmlAlignmentSummary object
	#
	foreach my $key (keys %{$attr}){
	    $bsml_alignment_summary->addattr( $key, $attr->{$key} );
	}

	#
	# this will retrieve all <Attribute> elements and add them to the BsmlAlignmentSummary object
	#
	foreach my $BsmlAttr ( ($alignment_summary->children('Attribute') )){
	    my $attr  = $BsmlAttr->atts();
	    $bsml_alignment_summary->addBsmlAttr($attr->{'name'}, $attr->{'content'});
	}

	#
	# this will retrieve all <Link> elements and add them to the BsmlAlignmentSummary object
	#
	foreach my $BsmlLink ($alignment_summary->children('Link')){
	    my $attr = $BsmlLink->atts();
	    $bsml_alignment_summary->addBsmlLink($attr->{'ref'}, $attr->{'href'}, $attr->{'role'});
	}


	#
	# this will retrieve all the <Aligned-sequence> elements
	#
	foreach my $aligned_sequence ($alignment_summary->children('Aligned-sequence')){

	    #
	    # Creates and adds to BsmlAlignedSequence object
	    #
	    my $bsml_aligned_sequence = $bsml_alignment_summary->returnBsmlAlignedSequenceR($bsml_alignment_summary->addBsmlAlignedSequence());

	    #
	    # this will retrieve the all of the XML attributes associated to the <Aligned-sequence> element
	    # e.g.  "length", "seqnum", and "name" 
	    my $attr = $aligned_sequence->atts();	

	    #
	    # this will add <Aligned-sequence> XML attributes to the BsmlAlignedSequence object
	    #
	    foreach my $key (keys %{$attr}){
		$bsml_aligned_sequence->addattr( $key, $attr->{$key} );
	    }
	    
	    #
	    # this will retrieve all <Attribute> elements and add them to the BsmlAlignedSequence object
	    #
	    foreach my $BsmlAttr ( ($aligned_sequence->children('Attribute') )){
		my $attr  = $BsmlAttr->atts();
		$bsml_aligned_sequence->addBsmlAttr($attr->{'name'}, $attr->{'content'});
	    }

	    #
	    # this will retrieve all <Link> elements and add them to the BsmlAlignedSequence object
	    #	    
	    foreach my $BsmlLink ($aligned_sequence->children('Link')){
		my $attr = $BsmlLink->atts();
		$bsml_aligned_sequence->addBsmlLink($attr->{'ref'}, $attr->{'href'}, $attr->{'role'});
	    }
	}
    }

    #
    # Parsing the <Sequence-alignment> section e.g.:
    #
    # <Sequence-alignment sequences="1:2:">
    #      <Sequence-data seq-name="lma2.204.m00029_protein:1">lma2.204.m00029_protein      MHSQRPCAGY NARCEHDRSV LPSLMSQQDE FFTGSSPRTA SVHSEARPHL 
    #                                                          lma2.204.m00029_protein      QERFPFATAA SSLPEHSDGA NLYSGAGATN TRGNFKVAVR VRPPLHRELH 
    #                                                          lma2.204.m00029_protein      GYRPFVDVVQ IVPEHPNSIT LCDALDTEDG RGAVYSRQSY TFDRVYAADA 
    #      </Sequence-data>
    #      <Sequence-data seq-name="lma2.206.m00009_protein:2">lma2.206.m00009_protein      .MSLMPGYPV PRDQWVESAE ACAACSKRFT FFAFKENCPC CGRLFCSSCL 
    #                                                          lma2.206.m00009_protein      SAQCTLFPTA PPKAVCLDCF RKAQDWRLSQ LEQQQQQQAT TQADVG...A 
    #                                                          lma2.206.m00009_protein      AAAPLSTSME VLESKLSALE EEFDRAKANA R.HLREENDS LIDLLAAKDS 
    #      </Sequence-data>
    #      <Alignment-consensus></Alignment-consensus>
    # </Sequence-alignment>
    #
    #
    #

    #
    # this will retrieve all the <Sequence-alignment> elements
    #
    foreach my $sequence_alignment ($seq_aln->children('Sequence-alignment')){


	#
	# Creates and adds to BsmlSequenceAlignment object
	#
	my $bsml_sequence_alignment = $bsml_multi_aln->returnBsmlSequenceAlignmentR( $bsml_multi_aln->addBsmlSequenceAlignment() ); 
	
	
	#
	# this will retrieve the "sequences" XML attribute of the <Sequence-alignment> element
	#
	my $attr = $sequence_alignment->atts();

	#
	# this will add <Sequence-alignment> XML attributes to the BsmlSequenceAlignment object
	#
	foreach my $key (keys %{$attr}){
	    $bsml_sequence_alignment->addattr( $key, $attr->{$key} );
	}
		
	#
	# this will retrieve <Attribute> elements and add them to the BsmlSequenceAlignment object
	#
	foreach my $BsmlAttr ( ($sequence_alignment->children('Attribute') )){
	    my $attr  = $BsmlAttr->atts();
	    $bsml_sequence_alignment->addBsmlAttr($attr->{'name'}, $attr->{'content'});
	}
	
	#
	# this will retrieve <Link> elements and add them to the BsmlSequenceAlignemnt object
	#
	foreach my $BsmlLink ($sequence_alignment->children('Link')){
	    my $attr = $BsmlLink->atts();
	    $bsml_sequence_alignment->addBsmlLink($attr->{'ref'}, $attr->{'href'}, $attr->{'role'});
	}
	

	
	#
	# this will retrieve all the <Sequence-data> elements
	#
	foreach my $sequence_data ($sequence_alignment->children('Sequence-data')){
	 
	    
	    #
	    # Creates and adds to BsmlSequenceData object
	    #
	    my $bsml_sequence_data = $bsml_sequence_alignment->returnBsmlSequenceDataR($bsml_sequence_alignment->addBsmlSequenceData());

	
	    #
	    # this will retrieve all of the XML attributes associated to the <Sequence-data> element
	    # e.g. "seq-name"
	    #
	    my $attr = $sequence_data->atts();	

	    # add raw sequence data if found 
	    
	    my $seqDat = $sequence_data->text();
	    if( $seqDat ){
		$bsml_sequence_data->addSequenceAlignmentData( $seqDat );
	    }

	    #
	    # this will add <Sequence-data> XML attributes to the BsmlSequenceData object
	    #
	    foreach my $key (keys %{$attr}){
		$bsml_sequence_data->addattr( $key, $attr->{$key} );
	    }
	    
	    #
	    # this will retrieve all <Attribute> elements and add them to the BsmlSequenceData object
	    #
	    foreach my $BsmlAttr ( ($sequence_data->children('Attribute') )){
		my $attr  = $BsmlAttr->atts();
		$bsml_sequence_data->addBsmlAttr($attr->{'name'}, $attr->{'content'});
	    }
	    
	    #
	    # this will retrieve all <Link> elements and add them to the BsmlSequenceData object
	    #
	    foreach my $BsmlLink ($sequence_data->children('Link')){
		my $attr = $BsmlLink->atts();
		$bsml_sequence_data->addBsmlLink($attr->{'ref'}, $attr->{'href'}, $attr->{'role'});
	    }

#	    my $sequence_data_pcdata = $bsml_sequence_data->returnSequenceAlignmentData($bsml_sequence_data->addSequenceAlignmentData());

	}
	
	if (0){
	    #
	    # for the time being we are not using the <Alignment-consensus> since is derivable
	    #



	    #
	    # this will retrieve all the <Alignment-consensus> elements
	    #
	    foreach my $alignment_consensus ($sequence_alignment->children('alignment-consensus')){
		
		
		#
		# Creates and adds to BsmlAlignmentConsensus object
		#
		my $bsml_alignment_consensus = $bsml_multi_aln->returnBsmlAlignmentConsensusR($bsml_multi_aln->addBsmlAlignmentConsensus());
		
		
		#
		# this will retrieve the XML attributes associated to the <Alignment-consensus> element
		#
		my $attr = $alignment_consensus->atts();	
		
		#
		# this will add <Alignment-consensus> XML attributes to the BsmlAlignmentConsensus object
		#
		foreach my $key (keys %{$attr}){
		    $bsml_alignment_consensus->addattr( $key, $attr->{$key} );
		}
		
		
		#
		# this will retrieve all <Attribute> elements and add them to the BsmlAlignmentConsensus object
		#
		foreach my $BsmlAttr ( ($alignment_consensus->children('Attribute') )){
		    my $attr  = $BsmlAttr->atts();
		    $bsml_alignment_consensus->addBsmlAttr($attr->{'name'}, $attr->{'content'});
		}
		
		#
		# this will retrieve all <Link> elements and add  them to the BsmlAlignmentConsensus object
		#
		foreach my $BsmlLink ($alignment_consensus->children('Link')){
		    my $attr = $BsmlLink->atts();
		    $bsml_alignment_consensus->addBsmlLink($attr->{'ref'}, $attr->{'href'}, $attr->{'role'});
		}
	    }
	}
    }


    foreach my $crossRef( $seq_aln->children( 'Cross-reference' ))
    {

	my $attr = $crossRef->atts();

	my $xref = $bsml_multi_aln->returnBsmlCrossReferenceR( $bsml_multi_aln->addBsmlCrossReference );
	
	foreach my $key ( keys( %{$attr} ) )
        {
	    $xref->addattr( $key, $attr->{$key});
        }
    }     
    $twig->purge();

    return $bsml_multi_aln;
    
}
1
