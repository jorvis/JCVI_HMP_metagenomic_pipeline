package BSML::BsmlBuilder;
@ISA = qw( BSML::BsmlDoc );

=head1 NAME

  BsmlBuilder.pm - class to facilitate the creation of BSML documents using a similar
  interface to that provided by the Labbook Bsml APIs Java class BsmlBuilder.

=head1 VERSION

  This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS

use BsmlBuilder;

my $seqdat = 'agctagctagctagctagctagct';
my $doc = new BsmlBuilder;
my $seq = $doc->createAndAddSequence( '_bsml001', 'test_basic_sequence', length( $seqdat ), 'dna' );
my $seq2 = $doc->createAndAddExtendedSequenceN( id => '_bsml002', title => 'test_extended_sequence', length => '24', molecule => 'dna', topology => 'linear' );
$doc->createAndAddSeqData( $seq, $seqdat );
$doc->write( 'output_file.xml' );

=cut

use strict;
use warnings;
BEGIN {
use BSML::BsmlDoc;
use BSML::Logger;
}
use Data::Dumper;

my $logger = BSML::Logger::get_logger("Logger::BSML");


# class variable used to create unique element identifiers
my $elem_id = 0;

=item $builder->createAndAddSequence( $id, $title, $length, $molecule, $class )

B<Description:> Creates a simple Bsml sequence containing minimum attributes

B<Parameters:> $id - the identifier of the Bsml element (unique at the document level)
  $title - the title of the sequence
  $length - the sequence length
  $molecule - the type of molecule from the controled vocabulary ( mol-not-set, dna, rna, aa, na, other-mol ) the default is mol-net-set

B<Returns:> a sequence object reference

=cut 

sub createAndAddSequence
{
  my $self = shift;
  my ( $id, $title, $length, $molecule, $class ) = @_;

  if( !($id) )
    {
	$id = $self->getUID();
    }

  if( !($molecule) )
    {
      $molecule = 'mol-not-set';
    }
  else
    {
      if( !(($molecule eq 'mol-not-set') || ($molecule eq 'dna') || ($molecule eq 'rna') || ($molecule eq 'aa') || ($molecule eq 'na') || ($molecule eq 'other-mol')) )
	{
	  $molecule = 'mol-not-set';
	}
    }

  my $seq = $self->returnBsmlSequenceR( $self->addBsmlSequence() );
  $seq->setattr( 'id', $id );
  $seq->setattr( 'title', $title );
  $seq->setattr( 'molecule', $molecule );

  if( defined($length) && !($length eq '') ){ $seq->setattr( 'length', $length ); }
  if( defined($class)  && !($class eq '') ){ $seq->setattr( 'class', $class ); }

  #Put a reference to the new sequence in the document namespace

  BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $seq );

  #return a reference to the new sequence

  return $seq;
}

=item $builder->createAndAddExtendedSequence( $id, $title, $length, $molecule, $locus, $dbsource, $icAcckey, $topology, $strand, $class )

B<Description:> Add a new sequence to the document with extended attributes

B<Parameters:> 
  $id - document wide unique identifier for the sequence
  $title - text description of sequence
  $length - integer value representing the total length of the sequence
  $molecule - type of molecule ( mol-not-set, dna, rna, aa, na, other-mol )
  $locus - sequence name
  $dbsource - database source of the sequence
  $icAcckey - internationl collboration accession number
  $topology - molecule shape (top-not-set, linear, circular, tandem, top-other )
  $strand - (std-not-set, ss, ds, mixed, std-other)

B<Returns:>
  A reference to the created sequence

=cut

sub createAndAddExtendedSequence
  {
    my $self = shift;
    my ( $id, $title, $length, $molecule, $locus, $dbsource, $icAcckey, $topology, $strand, $class ) = @_;

    my $seq = $self->createAndAddSequence( $id, $title, $length, $molecule );

    if( $locus    ){ $seq->setattr('locus', $locus);}
    if( $dbsource ){ $seq->setattr('dbsource', $dbsource);}
    if( $icAcckey ){ $seq->setattr('ic-acckey', $icAcckey);}
    if( defined($class)  && !($class eq '') ){ $seq->setattr( 'class', $class ); }

    if( $topology ){ 
      if( $topology eq 'top-not-set' || $topology eq 'linear' || $topology eq 'circular' || $topology eq 'tandem' || $topology eq 'top-other' ){
	$seq->setattr('topology', $topology);}
      else{
	$topology = 'top-not-set';
	$seq->setattr('topology', $topology);}
    }

    if( $strand ){ 
      if( $strand eq 'std-not-set' || $strand eq 'ss' || $strand eq 'ds' || $strand eq 'mixed' || $strand eq 'std-other' ){
	$seq->setattr('strand', $strand);}
      else{
	$strand = 'std-not-set';
	$seq->setattr('strand', $strand);}
    }
		
    #return a reference to the created sequence

    return $seq;
  }

=item $builder->createAndAddExtendedSequenceN( key => value )

B<Description:>
  Add a new sequence to the document with named extended attributes

B<Parameters:> 
  Looks for the following keys corresponding to the descriptions given above.
  {id, title, length, molecule, locus, dbsource, icAcckey, topology, strand}

B<Returns:>
  A reference to the created sequence

=cut

sub createAndAddExtendedSequenceN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddExtendedSequence( $args{'id'}, $args{'title'}, $args{'length'}, $args{'molecule'}, $args{'locus'}, $args{'dbsource'}, $args{'icAcckey'}, $args{'topology'}, $args{'strand'}, $args{'class'});

  }

sub createAndAddFeatureTables
  {
    #not used in the Bsml API
  }

sub createAndAddFeatureTable
  {
    my $self = shift;
    my ( $seq, $id, $title, $class ) = @_;

    if( !($id) )
      {
	  $id = $self->getUID();
      }

    if( ref($seq) eq 'BSML::BsmlSequence' )
      {
	my $FTable = $seq->returnBsmlFeatureTableR( $seq->addBsmlFeatureTable() );
	$FTable->setattr( 'id', $id );
	if( $title ){ $FTable->setattr('title', $title); }
	if( $class ){ $FTable->setattr('class', $class); }

	#Put a reference to the new feature table in the document namespace
	BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $FTable );

	return $FTable;
      }
    else
      {
	my $sequences = $self->returnBsmlSequenceListR();

	foreach my $seqR ( @{$sequences} )
	  {
	    if( $seqR->returnattr('id') eq $seq )
	      {
		my $FTable = $seqR->returnBsmlFeatureTableR( $seqR->addBsmlFeatureTable() );
		$FTable->setattr( 'id', $id );
		if( $title ){ $FTable->setattr('title', $title); }
		if( $class ){ $FTable->setattr('class', $class); }
		
		#Put a reference to the new feature table in the document namespace
		BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $FTable );

		return $FTable;
	      }
	  }
      }
  }

sub createAndAddFeatureTableN
  {
    my $self = shift;

    my %args = @_;

    return $self->createAndAddFeatureTable( $args{'seq'}, $args{'id'}, $args{'title'}, $args{'class'} );
  }

sub createAndAddReference
  {
    my $self = shift;
    my ( $FTable, $refID, $refAuthors, $refTitle, $refJournal, $dbxref ) = @_;

    if( !($refID) )
      {
	$refID = "Bsml"."$elem_id";
	$elem_id++;
      }

     if( ref($FTable) eq 'BSML::BsmlFeatureTable' )
      {
	my $rref = $FTable->returnBsmlReferenceR( $FTable->addBsmlReference() );
	$rref->setattr( 'id', $refID );

	if( defined($refAuthors) && !($refAuthors eq '') ){ $rref->addBsmlRefAuthors( $refAuthors ); }
	if( defined($refTitle) && !($refTitle eq '')){ $rref->addBsmlRefTitle( $refTitle ); }
	if( defined($refJournal) && !($refJournal eq '')){ $rref->addBsmlRefJournal( $refJournal ); }

	if( defined($dbxref) && !($dbxref eq '' )){ $rref->setattr( 'dbxref', $dbxref ); }

	#Put a reference to the new Bsml reference in the document namespace
	BSML::BsmlDoc::BsmlSetDocumentLookup( $refID, $rref );

	return $rref;
      }
    else
      {
	my $seqs =  $self->returnBsmlSequenceListR();

	foreach my $seq ( @{$seqs} )
	  {
	    foreach my $FeatureTable ( @{$seq->returnBsmlFeatureTableListR()} )
	      {
		if( $FeatureTable->returnattr( 'id' )){
		  #if( $FeatureTable->returnattr('id') eq $FTable )
		  {
		    my $rref = $FeatureTable->returnBsmlReferenceR( $FeatureTable->addBsmlReference() );
		    
		    $rref->setattr( 'id', $refID );

		    if( defined($refAuthors) && !($refAuthors eq '') ){ $rref->addBsmlRefAuthors( $refAuthors ); }
		    if( defined($refTitle) && !($refTitle eq '')){ $rref->addBsmlRefTitle( $refTitle ); }
		    if( defined($refJournal) && !($refJournal eq '')){ $rref->addBsmlRefJournal( $refJournal ); }

		    if( defined($dbxref) && !($dbxref eq '' )){ $rref->setattr( 'dbxref', $dbxref ); }

		    #Put a reference to the new Bsml reference in the document namespace
		    BSML::BsmlDoc::BsmlSetDocumentLookup( $refID, $rref );

		    return $rref;
		  }
		}
	      }
	  }
      }
  }

sub createAndAddReferenceN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddReferenceN( $args{'FTable'}, $args{'refID'}, $args{'refAuthors'}, $args{'refTitle'}, $args{'refJournal'}, $args{'dbxref'} );
  }

sub createAndAddFeature
  {
    my $self = shift;
    my ( $FTable, $id, $title, $class, $comment, $displayAuto ) = @_;


     if( !($id) )
      {
	  $id = $self->getUID();
      }

      if( ref($FTable) eq 'BSML::BsmlFeatureTable' )
      {
	my $fref = $FTable->returnBsmlFeatureR( $FTable->addBsmlFeature() );

	$fref->setattr( 'id', $id );

	if( defined($title) && !($title eq '') ){ $fref->setattr( 'title', $title ); }
	if( defined($class) && !($class eq '')){ $fref->setattr('class', $class); }
	if( defined($comment) && !($comment eq '')){ $fref->setattr('comment', $comment); }

	if( defined($displayAuto) && !($displayAuto eq '' )){ $fref->setattr( 'display-auto', $displayAuto ); }

	#Put a reference to the new Bsml feature in the document namespace
	BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $fref );

	return $fref;
      }
    else
      {
	my $seqs =  $self->returnBsmlSequenceListR();

	foreach my $seq ( @{$seqs} )
	  {
	    foreach my $FeatureTable (@{$seq->returnBsmlFeatureTableListR()})
	      {
		if( $FeatureTable->returnattr( 'id' ) eq $FTable )
		  {
		    my $fref = $FeatureTable->returnBsmlFeatureR( $FeatureTable->addBsmlFeature() );

		    
		    $fref->setattr( 'id', $id );

		    if( defined($title) && !($title eq '') ){ $fref->setattr( 'title', $title ); }
		    if( defined($class) && !($class eq '')){ $fref->setattr('class', $class); }
		    if( defined($comment) && !($comment eq '')){ $fref->setattr('comment', $comment); }

		    if( defined($displayAuto) && !($displayAuto eq '' )){ $fref->setattr( 'display-auto', $displayAuto ); }

		    #Put a reference to the new Bsml feature in the document namespace
		    BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $fref );

		    return $fref;
		  }
	      }
	  }
	}
  }

sub createAndAddFeatureN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddFeature( $args{'FTable'}, $args{'id'}, $args{'$title'}, $args{'class'}, $args{'comment'}, $args{'displayAuto'} );
  }

# see bug 2307
# Added $startopen and $endopen support
sub createAndAddFeatureWithLoc
  {
    my $self = shift;
    my ( $FTable, $id, $title, $class, $comment, $displayAuto, $start, $end, $complement, $startopen, $endopen ) = @_;

    my $feature = $self->createAndAddFeature( $FTable, $id, $title, $class, $comment, $displayAuto );
    

#     if (($start == 0) or ($end == 0)){
# 	$logger->fatal("id '$id' start '$start' end '$end'");
#     }


    if( (defined($start)) && (defined($end)) && ($start ne "") && ($end ne "") )
      {
	if( $start == $end )
	  {
	    #add a site position to the feature
	    $feature->addBsmlSiteLoc( $start, $complement );
	  }
	else
	  {
	    #add an interval location to the feature
	    $feature->addBsmlIntervalLoc( $start, $end, $complement, $startopen, $endopen );
	  }
      }
    
    return $feature;     
  }

# see bug 2307
# Added $startopen and $endopen support
sub createAndAddFeatureWithLocN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddFeatureWithLoc( $args{'FTable'}, $args{'id'}, $args{'title'}, $args{'class'}, $args{'comment'}, $args{'displayAuto'}, $args{'start'}, $args{'end'}, $args{'complement'}, $args{'startopen'}, $args{'endopen'} )
  }

# see bug 2307
# Added $startopen and $endopen support
sub createAndAddIntervalLoc
  {
    my $self = shift;
    my ( $feature, $start, $end, $complement, $startopen, $endopen ) = @_;

    if( ref($feature) eq 'BSML::BsmlFeature' )
      {
	$feature->addBsmlIntervalLoc( $start, $end, $complement, $startopen, $endopen );
      }
    else
      {
	my $seqs =  $self->returnBsmlSequenceListR();

	foreach my $seq ( @{$seqs} )
	  {
	    foreach my $FeatureTable ( @{$seq->returnBsmlFeatureTableListR()} )
	      {
		foreach my $rFeature( @{$FeatureTable->returnBsmlFeatureListR()} )
		  {
		    if( $rFeature->returnattr( 'id' ) eq $feature )
		      {
			$rFeature->addBsmlIntervalLoc( $start, $end, $complement, $startopen, $endopen );
		      }
		  }
	      }
	  }
      }

    return $feature;
  }

# see bug 2307
# Added $startopen and $endopen support
sub createAndAddIntervalLocN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddIntervalLoc( $args{'feature'}, $args{'start'}, $args{'end'}, $args{'complement'}, $args{'startopen'}, $args{'endopen'});
  }

sub createAndAddSiteLoc
  {
    my $self = shift;
    my ( $feature, $site, $complement, $class ) = @_;

    if( ref($feature) eq 'BSML::BsmlFeature' )
      {
	$feature->addBsmlSiteLoc( $site, $complement, $class );
      }
    else
      {
	my $seqs =  $self->returnBsmlSequenceListR();

	foreach my $seq ( @{$seqs} )
	  {
	    foreach my $FeatureTable ( @{$seq->returnBsmlFeatureTableListR()} )
	      {
		foreach my $rFeature( @{$FeatureTable->returnBsmlFeatureListR()} )
		  {
		    if( $rFeature->returnattr( 'id' ) eq $feature )
		      {
			$rFeature->addBsmlSiteLoc( $site, $complement, $class );
		      }
		  }
	      }
	  }
      }

    return $feature;
  }

sub createAndAddSiteLocN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddSiteLoc( $args{'feature'}, $args{'site'}, $args{'complement'} );
  }

sub createAndAddQualifier
  {
    my $self = shift;
    my ( $feature, $valuetype, $value ) = @_;

    if( ref($feature) eq 'BSML::BsmlFeature' )
      {
	$feature->addBsmlQualifier( $valuetype, $value );
      }
    else
      {
	my $seqs =  $self->returnBsmlSequenceListR();

	foreach my $seq ( @{$seqs} )
	  {
	    foreach my $FeatureTable ( @{$seq->returnBsmlFeatureTableListR()} )
	      {
		foreach my $rFeature( @{$FeatureTable->returnBsmlFeatureListR()} )
		  {
		    if( $rFeature->returnattr( 'id' ) eq $feature )
		      {
			$rFeature->addBsmlQualifier( $valuetype, $value );
		      }
		  }
	      }
	  }
      }

    return $feature;
  }

sub createAndAddQualifierN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddQualifier( $args{'feature'}, $args{'valuetype'}, $args{'value'} )
  }

sub createAndAddBsmlAttribute
  {
    my $self = shift;
    my ($elem, $key, $value ) = @_;

    $elem->addBsmlAttr( $key, $value );
  }

sub createAndAddBsmlAttributes
  {
    my $self = shift;
    my ($elem, %atts ) = @_;

    for my $key (keys %atts) {
        $elem->addBsmlAttr( $key, $atts{$key} );
    }
  }

sub createAndAddBsmlAttributeN
  {
    my $self = shift;
    my %args = @_;

    my $elem = $args{'elem'};
    $elem->addBsmlAttr( $args{'key'}, $args{'value'} ); 
  }

sub createAndAddAttribute
  {
    my $self = shift;
    my ($elem, $key, $value) = @_;

    $elem->addattr( $key, $value );
  }

sub createAndAddAttributeN
  {
    my $self = shift;
    my %args = @_;

    my $elem = $args{'elem'};
    $elem->addattr( $args{'key'}, $args{'value'} );
  }

sub createAndAddLink
  {
    my $self = shift;
    my ($elem, $rel, $href, $role) = @_;

    $elem->addBsmlLink( $rel, $href, $role );
  }

sub createAndAddLinkN
  {
    my $self = shift;
    my %args = @_;
    my $elem = $args{'elem'};
    $elem->addBsmlLink( $args{'title'}, $args{'href'} );
  }

sub createAndAddSeqData
  {
    my $self = shift;
    my ( $seq, $seqdat ) = @_;

    if( ref($seq) eq 'BSML::BsmlSequence' )
      {
	$seq->addBsmlSeqData( $seqdat );	
	return $seq;
      }
    else
      {
	my $sequences = $self->returnBsmlSequenceListR();
	
	foreach my $seqR ( @{$sequences} )
	  {
	    if( $seqR->returnattr('id') eq $seq )
	      {
		$seqR->addBsmlSeqData( $seqdat );
		
		return $seqR;
	      }
	  }
      }
  }



sub createAndAddSeqDataN
  {
    my $self = shift;
    my %args = @_;
    
    return $self->createAndAddSeqData( $args{'seq'}, $args{'seqdat'} );
  }

sub createAndAddSeqDataImport
  {
    my $self = shift;
    my ( $seq, $format, $source, $id, $identifier ) = @_;

    if( !($id) )
    {
	$id = $self->getUID();
    }
    
    if( ref($seq) eq 'BSML::BsmlSequence' )
      {
	$seq->addBsmlSeqDataImport( $format, $source, $id, $identifier );	
	return $seq;
      }
    else
      {
	my $sequences = $self->returnBsmlSequenceListR();
	
	foreach my $seqR ( @{$sequences} )
	  {
	    if( $seqR->returnattr('id') eq $seq )
	      {
		$seqR->addBsmlSeqDataImport( $format, $source, $id, $identifier );
		
		return $seqR;
	      }
	  }
      }
  }

sub createAndAddSeqDataImportN
  {
    my $self = shift;
    my %args = @_;
    
    return $self->createAndAddSeqDataImport( $args{'seq'}, $args{'format'}, $args{'source'}, $args{'id'}, $args{'identifier'} );
  }


# Add a pairwise alignment object to the document.  

sub createAndAddSequencePairAlignment
{
    my $self = shift;
    my %args = @_;

    ## the lines below commented out to prevent this method from checking to see if an alignment
    ##  object already exists for this refseq->compseq combination. Tools such as AAT can have multiple
    ##  separate alignment chains between the same two sequences.

#    #determine if a sequence pair alignment for the query and dbmatch already exists in the document
#
#    my $alignment_pair_list = BSML::BsmlDoc::BsmlReturnAlignmentLookup( "$args{'refseq'}", "$args{'compseq'}" );
#
#    if( $alignment_pair_list && !( $args{'force'} )){
#
#	# if coordinates specifying the search window for the reference sequence are provided, the object layer will
#	# search for an alignment object having the same coordinates and return it if found. Otherwise a new alignment
#	# object will be created. This functionality is intended to facilitate Doug's blast use case.
#
#	if( ($args{'refstart'} ne "") && ($args{'refend'} ne "") )
#	{
#	    foreach my $aln ( @{$alignment_pair_list} )
#	    {
#		if( $args{'refstart'} eq $aln->returnattr( 'refstart' ) && $args{'refend'} eq $aln->returnattr( 'refend' ) )
#		{
#		    return $aln;
#		}   
#	    }
#	}
#	else
#	{
#	    # Else return the first alignment pair to retain legacy compatitibility
#
#	    return $alignment_pair_list->[0];
#	}
#    }

    #no alignment pair matches, add a new alignment pair
    
    #check to see if sequences exist in the BsmlDoc, if not add them with basic attributes
    
    if( !( $self->returnBsmlSequenceByIDR( "$args{'refseq'}")) ){
	$self->createAndAddSequence( "$args{'refseq'}", "$args{'refseq'}", $args{'reflength'}, '',"$args{'class'}");}
    
    if( !( $self->returnBsmlSequenceByIDR( "$args{'compseq'}")) ){
	$self->createAndAddSequence( "$args{'compseq'}", "$args{'compseq'}", $args{'complength'}, '', "$args{'class'}" );}
    
    my $alignment_pair = $self->returnBsmlSeqPairAlignmentR( $self->addBsmlSeqPairAlignment() );
    
    $alignment_pair->setattr( 'refseq', "$args{'refseq'}" );
    $alignment_pair->setattr( 'compseq', "$args{'compseq'}" );
    
    BSML::BsmlDoc::BsmlSetAlignmentLookup( "$args{'refseq'}", "$args{'compseq'}", $alignment_pair );
    
    $alignment_pair->setattr( 'refxref', $args{'refxref'});
    $alignment_pair->setattr( 'refstart', $args{'refstart'} );
    $alignment_pair->setattr( 'refend', $args{'refend'} );
    $alignment_pair->setattr( 'reflength', $args{'reflength'} );
    $alignment_pair->setattr( 'compxref', $args{'compxref'});
    $alignment_pair->setattr( 'compstart', $args{'compstart'} );
    $alignment_pair->setattr( 'compend', $args{'compend'} );
    $alignment_pair->setattr( 'complength', $args{'complength'} );
    $alignment_pair->setattr( 'method', $args{'method'} );

    ## add a class if the user passed it
    if ($args{class}) {
        $alignment_pair->setattr( 'class', $args{class} );
    }

    return $alignment_pair;
}


sub createAndAddSequencePairRun
{
    my $self = shift;
    my %args = @_;

    my $alignment_pair = $args{'alignment_pair'};

    if( ref( $alignment_pair) eq 'BSML::BsmlSeqPairAlignment' ) 
    {
	#add a new BsmlSeqPairRun to the alignment pair and return
	my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );
	
	$seq_run->setattr( 'refpos', $args{'refpos'} );
	$seq_run->setattr( 'runlength', $args{'runlength'} );
	$seq_run->setattr( 'refcomplement', $args{'refcomplement'});
	
	$seq_run->setattr( 'comppos', $args{'comppos'} );
	$seq_run->setattr( 'comprunlength', $args{'comprunlength'} );
	$seq_run->setattr( 'compcomplement', $args{'compcomplement'} );
	
	$seq_run->setattr( 'runscore', $args{'runscore'} );
	$seq_run->setattr( 'runprob', $args{'runprob'} );

	$seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} );
	$seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} );
	$seq_run->addBsmlAttr( 'chain_number', $args{'chain_number'} );
	$seq_run->addBsmlAttr( 'segment_number', $args{'segment_number'} );
	$seq_run->addBsmlAttr( 'p_value', $args{'p_value'} );

	# The seq pair run object is returned these client specific attributes should be 
	# appended by the client. Potential the first five above, as well...

	$seq_run->addBsmlAttr( 'PEffect_Cluster_Id', $args{'PEffect_Cluster_Id'} );
	$seq_run->addBsmlAttr( 'PEffect_Cluster_Gap_Count', $args{'PEffect_Cluster_Gap_Count'} );
	$seq_run->addBsmlAttr( 'PEffect_Cluster_Gene_Count', $args{'PEffect_Cluster_Gene_Count'} );

	$seq_run->addBsmlAttr( 'Mummer_percent_coverage_ref', $args{'Mummer_percent_coverge_ref'} );
	$seq_run->addBsmlAttr( 'Mummer_percent_coverage_comp', $args{'Mummer_percent_coverage_comp'} );

    #If the user passes a class, add it.
    if($args{'class'}) {
        $seq_run->setattr('class', $args{'class'});
    }
	
	return $seq_run;
    }
}

# Utility to add pairwise alignments specified in the btab format. 

sub createAndAddBtabLineN
  {
    my $self = shift;
    my %args = @_;

    #determine if the query name and the dbmatch name are a unique pair in the document

    my $alignment_pair_list = BSML::BsmlDoc::BsmlReturnAlignmentLookup( "$args{'query_name'}", "$args{'dbmatch_accession'}" );

    my $alignment_pair = '';
    if( $alignment_pair_list )
    {
	$alignment_pair = $alignment_pair_list->[0];
    }

    if( $alignment_pair  )
	  {
	    #add a new BsmlSeqPairRun to the alignment pair and return
	    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

	    if( $args{'start_query'} > $args{'stop_query'} )
	      {
		$seq_run->setattr( 'refpos', $args{'stop_query'} );
		$seq_run->setattr( 'runlength', $args{'start_query'} - $args{'stop_query'} + 1 );
		$seq_run->setattr( 'refcomplement', 1 );
	      }
	    else
	      {
		$seq_run->setattr( 'refpos', $args{'start_query'} );
		$seq_run->setattr( 'runlength', $args{'stop_query'} - $args{'start_query'} + 1 );
		$seq_run->setattr( 'refcomplement', 0 );
	      }

	    #the database sequence is always 5' to 3'

	    $seq_run->setattr( 'comppos', $args{'start_hit'} );
	    $seq_run->setattr( 'comprunlength', $args{'stop_hit'} - $args{'start_hit'} + 1 );
	    $seq_run->setattr( 'compcomplement', 0 );

	    $seq_run->setattr( 'runscore', $args{'bit_score'} );
	    $seq_run->setattr( 'runprob', $args{'e_value'} );

	    $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} );
	    $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} );
	    $seq_run->addBsmlAttr( 'chain_number', $args{'chain_number'} );
	    $seq_run->addBsmlAttr( 'segment_number', $args{'segment_number'} );
	    $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} );

	    return $alignment_pair;
	  }

    #no alignment pair matches, add a new alignment pair and sequence run

    #check to see if sequences exist in the BsmlDoc, if not add them with basic attributes

    if( !( $self->returnBsmlSequenceByIDR( "$args{'query_name'}")) ){
      $self->createAndAddSequence( "$args{'query_name'}", "$args{'query_name'}", $args{'query_length'}, 'aa', "$args{'class'}" );}

    if( !( $self->returnBsmlSequenceByIDR( "$args{'dbmatch_accession'}")) ){
      $self->createAndAddSequence( "$args{'dbmatch_accession'}", "$args{'dbmatch_accession'}", '', 'aa', "$args{'class'}" );}

    $alignment_pair = $self->returnBsmlSeqPairAlignmentR( $self->addBsmlSeqPairAlignment() );
    

    $alignment_pair->setattr( 'refseq', "$args{'query_name'}" );
    $alignment_pair->setattr( 'compseq', "$args{'dbmatch_accession'}" );

    BSML::BsmlDoc::BsmlSetAlignmentLookup( "$args{'query_name'}", "$args{'dbmatch_accession'}", $alignment_pair );

    $alignment_pair->setattr( 'refxref', ':'.$args{'query_name'});
    $alignment_pair->setattr( 'refstart', 0 );
    $alignment_pair->setattr( 'refend', $args{'query_length'} - 1 );
    $alignment_pair->setattr( 'reflength', $args{'query_length'} );

    $alignment_pair->setattr( 'method', $args{'blast_program'} );

    $alignment_pair->setattr( 'compxref', $args{'search_database'}.':'.$args{'dbmatch_accession'} );

    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

    if( $args{'start_query'} > $args{'stop_query'} )
      {
	$seq_run->setattr( 'refpos', $args{'stop_query'} );
	$seq_run->setattr( 'runlength', $args{'start_query'} - $args{'stop_query'} + 1 );
	$seq_run->setattr( 'refcomplement', 1 );
      }
    else
      {
	$seq_run->setattr( 'refpos', $args{'start_query'} );
	$seq_run->setattr( 'runlength', $args{'stop_query'} - $args{'start_query'} + 1 );
	$seq_run->setattr( 'refcomplement', 0 );
      }

    #the database sequence is always 5' to 3'
    
    $seq_run->setattr( 'comppos', $args{'start_hit'} );
    $seq_run->setattr( 'comprunlength', ($args{'stop_hit'} - $args{'start_hit'} + 1));
    $seq_run->setattr( 'compcomplement', 0 );
    
    $seq_run->setattr( 'runscore', $args{'bit_score'} );
    $seq_run->setattr( 'runprob', $args{'e_value'} );

    $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} );
    $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} );
    $seq_run->addBsmlAttr( 'chain_number', $args{'chain_number'} );
    $seq_run->addBsmlAttr( 'segment_number', $args{'segment_number'} );
    $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} );

    return $alignment_pair;
  }

sub createAndAddFeatureGroup
  {
    my $self = shift;
    my ($seq, $id, $groupset) = @_;

    if( !($id) )
      {
	$id = $self->getUID();
      }
    
    my $FeatureGroup = $seq->returnBsmlFeatureGroupR( $seq->addBsmlFeatureGroup() );
    $FeatureGroup->setattr( 'id', $id ); 

    if( ($groupset) )
      {
	$FeatureGroup->setattr( 'group-set', $groupset );
	BSML::BsmlDoc::BsmlSetFeatureGroupLookup( $groupset, $FeatureGroup );
      }

    BSML::BsmlDoc::BsmlSetDocumentLookup( $id, $FeatureGroup );

    return $FeatureGroup;
  }

sub createAndAddFeatureGroupN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddFeatureGroupN( $args{'seq'}, $args{'id'}, $args{'groupset'} );
  }

sub createAndAddFeatureGroupMember
  {
    my $self = shift;
    my ($FeatureGroup, $featref, $feattype, $grouptype, $cdata) = @_;

    $FeatureGroup->addBsmlFeatureGroupMember( $featref, $feattype, $grouptype, $cdata );
    
    return $FeatureGroup;
  }

sub createAndAddFeatureGroupMemberN
  {
    my $self = shift;
    my %args = @_;

    return $self->createAndAddFeatureGroupMember( $args{'featref'}, $args{'feattype'}, $args{'grouptype'}, $args{'cdata'} );
  }

sub createAndAddBsmlAttr
  {
    my $self = shift;
    my ($elm, $key, $value) = @_;

    $elm->addBsmlAttr( $key, $value );
    return $elm;
  }

sub createAndAddAttr
  {
    my $self = shift;
    my ($elm, $key, $value) = @_;

    $elm->addattr( $key, $value );
    return $elm;
  }

# First attempts at encoding search parameters in a Bsml Analysis element. Better use of the 
# Bsml DTD may be required. 
sub createAndAddAnalysis
{
    my $self = shift;
    my %args = @_;

    my $analysis = $self->returnBsmlAnalysisR( $self->addBsmlAnalysis() );

    $analysis->addattr( 'id', $args{'id'} );

    $analysis->addBsmlAttr( 'algorithm', $args{'algorithm'} );
    $analysis->addBsmlAttr( 'description', $args{'description'} );
    $analysis->addBsmlAttr( 'name', $args{'name'} );
    $analysis->addBsmlAttr( 'program', $args{'program'} );
    $analysis->addBsmlAttr( 'programversion', $args{'programversion'} );
    $analysis->addBsmlAttr( 'sourcename', $args{'sourcename'} );
    $analysis->addBsmlAttr( 'sourceuri', $args{'sourceuri'} );
    $analysis->addBsmlAttr( 'sourceversion', $args{'sourceversion'} );
    $analysis->addBsmlAttr( 'queryfeature_id', $args{'queryfeature_id'} );
    $analysis->addBsmlAttr( 'timeexecuted', $args{'timeexecuted'} );

    $analysis->addBsmlLink( $args{'bsml_link_relation'}, $args{'bsml_link_url'} );

    return $analysis;
}

# Support for multiple sequence alignments

sub createAndAddMultipleAlignmentTable
{
    my $self = shift;
    my %args = @_;

    my $multipleAlignmentTable = $self->returnBsmlMultipleAlignmentTableR( $self->addBsmlMultipleAlignmentTable() );

    $multipleAlignmentTable->addattr( 'id', $args{'id'} );
    $multipleAlignmentTable->addattr( 'molecule-type', $args{'molecule-type'} );
    $multipleAlignmentTable->addattr( 'seq-format', $args{'seq-format'} );

    if($args{'class'}) {
        $multipleAlignmentTable->addattr( 'class', $args{'class'} );
    }
    
    return $multipleAlignmentTable;
}

sub createAndAddAlignmentSummary
{
    my $self = shift;
    my %args = @_;

    my $table = $args{'multipleAlignmentTable'};
    my $seqType = $args{'seq-type'};
    my $seqFormat = $args{'seq-format'};
    
    my $alignmentSummary = $table->returnBsmlAlignmentSummaryR( $table->addBsmlAlignmentSummary() );

    $alignmentSummary->addattr( 'seq-type', $seqType );
    $alignmentSummary->addattr( 'seq-format', $seqFormat );

    return $alignmentSummary;
}

sub createAndAddAlignedSequence
{
    my $self = shift;
    my %args = @_;

    my $alignmentSummary = $args{'alignmentSummary'};
 
    my $alignedSequence = $alignmentSummary->returnBsmlAlignedSequenceR( $alignmentSummary->addBsmlAlignedSequence() );
    $alignedSequence->addattr( 'seqref', $args{'seqref'} );
    $alignedSequence->addattr( 'start', $args{'start'} );
    $alignedSequence->addattr( 'on-complement', $args{'on-complement'} );
    $alignedSequence->addattr( 'translated', $args{'translated'} );
    $alignedSequence->addattr( 'frame', $args{'frame'} );
    $alignedSequence->addattr( 'trans-table', $args{'trans-table'} );
    $alignedSequence->addattr( 'seqnum', $args{'seqnum'} );
    $alignedSequence->addattr( 'name', $args{'name'} );
    $alignedSequence->addattr( 'length', $args{'length'} );

    return $alignedSequence;    
}

sub createAndAddPairwiseAlignments
{
    my $self = shift;
    my %args = @_;

    my $table = $args{'multipleAlignmentTable'};
    
    my $PairwiseAlns = $table->returnBsmlPairwiseAlignmentsR( $table->addBsmlPairwiseAlignments() );

    return $PairwiseAlns;
}

sub createAndAddAlignedPair
{
    my $self = shift;
    my %args = @_;

    my $PairwiseAlns = $args{'pairwiseAlignments'};
    my $pair = $PairwiseAlns->returnBsmlAlignedPairR( $PairwiseAlns->addBsmlAlignedPair() );

    $pair->addattr( 'seqnum1', $args{'seqnum1'} );
    $pair->addattr( 'seqnum2', $args{'seqnum2'} );
    $pair->addattr( 'score', $args{'score'} );

    return $pair;
}

sub createAndAddSequenceAlignment
{
    my $self = shift;
    my %args = @_;

    my $table = $args{'multipleAlignmentTable'};

    my $seqAln = $table->returnBsmlSequenceAlignmentR( $table->addBsmlSequenceAlignment() );

    $seqAln->addattr( 'sequences', $args{'sequences'} );

    return $seqAln;
}

sub createAndAddSequenceData
{
    my $self = shift;
    my %args = @_;
    
    my $seqAln = $args{'sequenceAlignment'};

    my $seqDat = $seqAln->returnBsmlSequenceDataR( $seqAln->addBsmlSequenceData() );

    $seqDat->addattr( 'seq-name', $args{'seq-name'} );
    $seqDat->addSequenceAlignmentData( $args{'seq-data'} );

    return $seqDat;
}

sub createAndAddGenome
{
    my $self = shift;
    my %args = @_;

    my $genome = $self->returnBsmlGenomeR( $self->addBsmlGenome() );

       
    #
    # editor:    sundaram@tigr.org
    # date:      2005-08-17
    # bgzcase:   2051
    # URL:       http://serval.tigr.org:8080/bugzilla/show_bug.cgi?id=2051
    # comment:   The <Sequence> will now be explicitly linked with the <Genome>.
    #	         This hack matches what we did for the Cross-reference id assignment.

    if (defined($args{'id'})){
	
	my $id = $args{'id'};
	my $idfinal;

	if ($id =~ /^\d+/){
	    $idfinal = '_' . $id;
	}

	$genome->addattr( 'id', $idfinal) 
    }

    elsif( !($args{'id'}) ) {
	$args{'id'} = "_$elem_id";
	$elem_id++;
	$genome->addattr( 'id', $args{'id'} );		 
    }


    $genome->addattr( 'autosomal-chromosome-count', $args{'autosomal-chromosome-count'} );
    $genome->addattr( 'sex-chromosome-count', $args{'sex-chromosome-count'} );
    $genome->addattr( 'ploidy-count', $args{'ploidy-count'} );
    $genome->addattr( 'distinct-chromosome-count', $args{'distinct-chromosome-count'} );
    $genome->addattr( 'total-chromosome-count', $args{'total-chromosome-count'} );
    
    return $genome;
}

sub createAndAddOrganism
{
    my $self = shift;
    my %args = @_;

    my $genome = $args{ 'genome' };

    my $organism = $genome->returnBsmlOrganismR( $genome->addBsmlOrganism() );

    $organism->addattr( 'id', $args{'id'} );
    $organism->addattr( 'genus', $args{'genus'} );
    $organism->addattr( 'species', $args{'species'} );
    $organism->addattr( 'taxon-num', $args{'taxon-num'} );
    $organism->addattr( 'taxonomy', $args{'taxonomy'} );
    $organism->addattr( 'url', $args{'url'} );

    return $organism;    
}

sub createAndAddStrain
{
    my $self = shift;
    my %args = @_;

    my $organism = $args{'organism'};
    my $strain_name = $args{'name'};

    #-------------------------------------------------------------------------
    # editor:  sundaram@tigr.org
    # date:    2005-10-04
    # comment: The database and source_database information should each be 
    #          stored separately in a <Cross-reference> element.
    #
    #my $database = $args{'database'};
    #my $source_database = $args{'source_database'};

    my $strain = $organism->returnBsmlStrainR( $organism->addBsmlStrain() );

    $strain->addBsmlAttr( 'name', $strain_name );
    #$strain->addBsmlAttr( 'database', $database );
    #$strain->addBsmlAttr( 'source_database', $source_database );

    return $strain;

}

sub createAndAddSegmentSet
{
    my $self = shift;
    my %args = @_;

    my $segmentSetType = $args{'segment_set_type'};

    my $segmentSet = $self->returnBsmlSegmentSetR( $self->addBsmlSegmentSet() );

    $segmentSet->addattr( 'seg-set-type', $segmentSetType );

    return $segmentSet;
}

sub createAndAddSegment
{
    my $self = shift;
    my %args = @_;

    my $segmentSet = $args{'segment_set'};
    my $segSourceType = $args{'seg_source_type'};
    my $segSource = $args{'seg_source'};
    my $segId = $args{'seg_id'};
    my $segURL = $args{'seg_url'};
    my $segRole = $args{'seg_role'};
    my $segStart = $args{'seg_start'};
    my $segEnd = $args{'seg_end'};
    my $segOnComplement = $args{'seg_on_complement'};

    my $segment = $segmentSet->returnBsmlSegmentR( $segmentSet->addBsmlSegment() );

    $segment->addattr( 'seg-source-type', $segSourceType );
    $segment->addattr( 'seg-source', $segSource );
    $segment->addattr( 'seg-id', $segId );
    $segment->addattr( 'seg-url', $segURL );
    $segment->addattr( 'seg-role', $segRole );
    $segment->addattr( 'seg-start', $segStart );
    $segment->addattr( 'seg-end', $segEnd );
    $segment->addattr( 'seg-on-complement', $segOnComplement );

    return $segment;
}

sub createAndAddNumbering
{
    my $self = shift;
    my %args = @_;

    my $numbering = $args{'seq'}->addBsmlNumbering();
    
    $numbering->addattr( 'seqref', $args{'seqref'} );
    $numbering->addattr( 'use-numbering', $args{'use_numbering'} );
    $numbering->addattr( 'type', $args{'type'} );
    $numbering->addattr( 'units', $args{'units'} );
    $numbering->addattr( 'a', $args{'a'} );
    $numbering->addattr( 'b', $args{'b'} );
    $numbering->addattr( 'dec-places', $args{'dec_places'} );
    $numbering->addattr( 'refnum', $args{'refnum'} );
    $numbering->addattr( 'has-zero', $args{'has_zero'} );
    $numbering->addattr( 'ascending', $args{'ascending'} );
    $numbering->addattr( 'names', $args{'names'} );
    $numbering->addattr(' from-aligns', $args{'from_aligns'} );
    $numbering->addattr( 'aligns', $args{'aligns'} );
    
    return $numbering;
}

sub createAndAddCrossReference
{
  my $self = shift;
  my %args = @_;

  my $parent = $args{ 'parent' };

  my $xref = $parent->returnBsmlCrossReferenceR( $parent->addBsmlCrossReference );


   if (defined($args{'id'})){
      
       my $id = $args{'id'};
       my $idfinal;
       if ($id =~ /^\d+/){
	   $idfinal = '_' . $id;
       }
       $xref->addattr( 'id', $idfinal) 
   }
  elsif( !($args{'id'}) ) {
      $args{'id'} = "_$elem_id";
      $elem_id++;
      $xref->addattr( 'id', $args{'id'} );		 
  }



  $xref->addattr( 'context', $args{'context'} );
  $xref->addattr( 'database', $args{'database'} );
  $xref->addattr( 'identifier', $args{'identifier'} );
  $xref->addattr( 'identifier-type', $args{'identifier-type'} );
  $xref->addattr( 'title', $args{'title'} );
  $xref->addattr( 'behavior', $args{'behavior'} );
  $xref->addattr( 'href', $args{'href'} );
  $xref->addattr( 'role', $args{'role'} );
  
  return $xref;

}

sub createAndAddCrossReferencesByParse {
    ## this looks for recognized header formats and creates appropriate
    ##  Cross-Reference elements within the passed Sequence.
    my $self = shift;
    my %args = @_;
    
    my $seq = $args{'sequence'};
    my $str = $args{'string'} || '';
    
    ## PANDA OMNI lines like: OMNI|NTL01TW0001|AAO44098.1|28476008|
    if ($str =~ /^OMNI\|(.+)\|(.+)\|(.+)\|$/) {
        $self->createAndAddCrossReference( parent => $seq, database => 'OMNI', 'identifier-type' => 'version', identifier => $2 );
        $self->createAndAddCrossReference( parent => $seq, database => 'OMNI', 'identifier-type' => 'gi', identifier => $3 );
    ## PANDA Swiss-Prot lines like: SP|Q88BK3|DNAA_PSESM
    } elsif ($str =~ /^SP\|(.+)\|(.+)$/) {
        $self->createAndAddCrossReference( parent => $seq, database => 'SP', 'identifier-type' => 'accession', identifier => $1 );
        $self->createAndAddCrossReference( parent => $seq, database => 'SP', 'identifier-type' => 'entry_name', identifier => $2 );
    ## PANDA GenBank lines like: GB|AAN14420.1|23193191|AF319579
    } elsif ($str =~ /^GB\|(.+)\|(.+)\|(.+)$/) {
        $self->createAndAddCrossReference( parent => $seq, database => 'GB', 'identifier-type' => 'version', identifier => $1 );
        $self->createAndAddCrossReference( parent => $seq, database => 'GB', 'identifier-type' => 'gi', identifier => $2 );
    ## PANDA Protein Research Foundation lines like: PRF|2203412A|1586345|2203412A
    } elsif ($str =~ /^PRF\|(.+)\|.+\|.+$/) {
        $self->createAndAddCrossReference( parent => $seq, database => 'PRF', 'identifier-type' => 'name', identifier => $1 );
    ## PANDA Brookhaven Protein Data Bank lines like: >PDB|1I3Q_A|14278324|1I3Q_A
    } elsif ($str =~ /^PDB\|(.+)\|(.+)\|.+$/) {
        $self->createAndAddCrossReference( parent => $seq, database => 'PDB', 'identifier-type' => 'version', identifier => $1 );
        $self->createAndAddCrossReference( parent => $seq, database => 'PDB', 'identifier-type' => 'gi', identifier => $2 );
    ## PANDA NBRF PIR lines like: PIR|T30335|T30335
    } elsif ($str =~ /^PIR\|(.+)\|.+$/) {
        $self->createAndAddCrossReference( parent => $seq, database => 'PIR', 'identifier-type' => 'entry', identifier => $1 );
    }
}

#
# Added 2004-11-02 Bugzilla case 1808
#

sub createAndAddAttributeListMember {
    
    my $self = shift;
    my ($AttributeList, $name, $content) = @_;
    
    $AttributeList->addBsmlAttributeListMember( $name, $content );
    
    return $AttributeList;
}

sub createAndAddAttributeListMemberN {


    my $self = shift;
    my %args = @_;
    
    return $self->createAndAddAttributeListMember( $args{'name'}, $args{'content'} );
}

sub createAndAddBsmlAttributeList {

    my $self = shift;
    my ( $elem, $id ) = @_;

    if (!defined($id)){
	$id = "Bsml"."$elem_id";
	$elem_id++;
    }

    $elem->addBsmlAttrlist();
}




1;
