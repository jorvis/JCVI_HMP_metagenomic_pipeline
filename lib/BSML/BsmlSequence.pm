package BSML::BsmlSequence;

@ISA = qw( BSML::BsmlElement );


=head1 NAME

  BsmlSequence.pm - Bsml API Object representing the Bsml Sequence Element 

=head1 VERSION
  
  This document refers to version 1.0 of the BSML Object Layer

=head1 Description

=head2 Overview

  The BsmlSequence class allows sequence data <Seq-dat> and feature tables to be added and manipulated.

=head2 Constructor and initialization

  Typically a BsmlSequence is created by the BsmlDoc object it is contained within and manipulated
  as a reference.

  my $doc = new BsmlDoc;
  my $seq = $doc->returnBsmlSequenceR( $doc->addBsmlSequence() );

=head2 Class and object methods

=over 4

=cut

use BSML::BsmlCrossReference;
use BSML::BsmlNumbering;
use BSML::BsmlElement;
use BSML::BsmlFeatureGroup;
use BSML::BsmlFeatureTable;
use BSML::Logger;
use BSML::Indexer::Fasta;
use strict;
use warnings;
use XML::Writer;
use Data::Dumper;
use File::Basename;

my $logger = BSML::Logger::get_logger("Logger::BSML");


# a bsml sequence stores raw sequence data and a list of feature tables

sub new
  {

      $logger->debug("") if $logger->is_debug;


    my $class = shift;
    my $self = {};
    bless $self, $class;
    
    $self->init();
    return $self;
  }

sub init
  {

      $logger->debug("") if $logger->is_debug;
    my $self = shift;

    $self->{ 'attr' } = {};
    $self->{ 'BsmlAttr' } = {};
    $self->{ 'BsmlFeatureTables' } = [];
    $self->{ 'BsmlFeatureGroups' } = [];
    $self->{ 'BsmlSeqData' } = undef;
    $self->{ 'BsmlSeqDataImport' } = {};
    $self->{ 'BsmlLink' } = [];
    $self->{ 'BsmlNumbering' } = undef;
    $self->{ 'BsmlCrossReference' } = [];
    $self->{ 'BsmlAttributeList' } = undef;
    
  }

=item $seq->addBsmlFeatureTable()

B<Description:> adds a feature table to the sequence object

B<Parameters:> None

B<Returns:> The index of the added feature table

=cut 

sub addBsmlFeatureTable
  {

      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    push( @{$self->{'BsmlFeatureTables'}}, new BSML::BsmlFeatureTable );

    my $index = @{$self->{'BsmlFeatureTables'}} - 1;
    return $index;
  }

=item $seq->dropBsmlFeatureTable( $index )

B<Description:> deletes a feature table from the sequence object

B<Parameters:> The index of the feature table to be deleted

B<Returns:> None

=cut 

sub dropBsmlFeatureTable
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlFeatureTables'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlFeatureTables'}[$i] );
	  }
      }

    $self->{'BsmlFeatureTables'} = \@newlist;    
  }

=item $seq->returnBsmlFeatureTableListR()

B<Description:> Return a list of references to all the feature table objects contained in the document.

B<Parameters:> None

B<Returns:> a list of BsmlFeatureTable object references

=cut 

sub returnBsmlFeatureTableListR
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    return $self->{'BsmlFeatureTables'};
  }

=item $seq->returnBsmlFeatureTableR( $index )

B<Description:> Return a reference to a feature table object given its index

B<Parameters:> ($index) - the feature table index returned from addBsmlFeatureTable (position of the table in the reference list)

B<Returns:> a BsmlFeatureTable object reference

=cut

sub returnBsmlFeatureTableR
  {

      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlFeatureTables'}[$index];
  }

=item $seq->addBsmlSeqData( $seq_string )

B<Description:> add a sequence string <Seq_dat> to the object

B<Parameters:> ( $seq_string ) - string containing raw sequence data

B<Returns:> None

=cut

sub addBsmlSeqData
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    ($self->{'BsmlSeqData'}) = @_; 
  }

=item $seq->setBsmlSeqData( $seq_string )

B<Description:> same as addBsmlSeqData - maintained to make methods consistent with the rest of the API

B<Parameters:> ( $seq_string ) - string containing raw sequence data

B<Returns:> None

=cut

sub setBsmlSeqData
  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($seq) = @_;

    $self->addBsmlSeqData( $seq );
  }

=item $seq->dropBsmlSeqData()

B<Description:> delete the raw sequence <Seq-dat> from the sequence object

B<Parameters:> None

B<Returns:> None

=cut

sub dropBsmlSeqData
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
   
    $self->{'BsmlSeqData'} = undef;
  }

=item $seq->returnSeqData

B<Description:> return a string containing the raw sequence from the sequence object

B<Parameters:> None

B<Returns:> a string containing raw sequence data

=cut

sub returnSeqData
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    return $self->{'BsmlSeqData'};
  }

sub addBsmlSeqDataImport
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($format, $source, $id, $identifier) = @_;
    my($fasta_file,$fasta_dir) = fileparse($source);  
      if(-e $fasta_file){
	  #index fasta file if it exists
	  my $indexer = BSML::Indexer::Fasta->new($fasta_file, $fasta_dir);
	  
	  # Check the health of the various indices for the data file.
	  my @check = $indexer->check_indices;
	  
	  # Create the indices if necessary...
	  if ($check[0] == 1) { $indexer->index_entries };
	  if ($check[1] == 1) { $indexer->index_headers };
    }  
    $self->{'BsmlSeqDataImport'}->{'format'} = $format;
    $self->{'BsmlSeqDataImport'}->{'source'} = $source;

    $self->{'BsmlSeqDataImport'}->{'identifier'} = $identifier;

    $self->{'BsmlSeqDataImport'}->{'id'} = $id;
  }

sub setBsmlSeqDataImport
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($format, $source, $id, $identifier) = @_;

    $self->addBsmlSeqDataImport( $format, $source, $id, $identifier );
  }

sub dropBsmlSeqDataImport
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    $self->{'BsmlSeqDataImport'} = {};
  }

sub returnBsmlSeqDataImport
  {
      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    return $self->{'BsmlSeqDataImport'};
  }

=item $seq->addBsmlFeatureGroup()

B<Description:> add a feature group to the sequence

B<Parameters:> None

B<Returns:> the index of the added feature group

=cut

sub addBsmlFeatureGroup
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    push( @{$self->{'BsmlFeatureGroups'}}, new BSML::BsmlFeatureGroup );

    my $index = @{$self->{'BsmlFeatureGroups'}} - 1;

    #In order to retain the relationship of genes to the assembly on which they
    #are contained in a memory efficient manner consistent with the document level
    #lookups, the sequence id is embedded in each feature group. 

    $self->{'BsmlFeatureGroups'}->[$index]->{'ParentSequenceId'} = $self->returnattr( 'id' );

    return $index;

  }

=item $seq->dropBsmlFeatureGroup( $index )

B<Description:> drop a feature group from the sequence

B<Parameters:> $index - the index of the feature group to be deleted

B<Returns:> None

=cut

sub dropBsmlFeatureGroup
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($index) = @_;

    my @newlist;

    for( my $i=0; $i<length(@{$self->{'BsmlFeatureGroups'}}); $i++ )
      {
	if( $i != $index )
	  {
	    push( @newlist, $self->{'BsmlFeatureGroups'}[$i] );
	  }
      }

    $self->{'BsmlFeatureGroups'} = \@newlist;    
  }

=item $seq->returnBsmlFeatureGroupListR()

B<Desciption:> returns the list of feature group references associated with the sequence

B<Parameters:> None

B<Returns:> a list of feature group references

=cut

sub returnBsmlFeatureGroupListR
  {
      $logger->debug("") if $logger->is_debug;
      
      my $self = shift;
      return $self->{'BsmlFeatureGroups'};
  }

=item $seq->returnFeatureGroupR( $index )

B<Description:> returns the feature group reference of index $index

B<Parameters:> $index - the index of the feature group to be returned

B<Returns:> a reference to a feature group

=cut

sub returnBsmlFeatureGroupR
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    my ($index) = @_;
    return $self->{'BsmlFeatureGroups'}[$index];
  } 


sub addBsmlNumbering
{
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    $self->{'BsmlNumbering'} = new BSML::BsmlNumbering;

    return $self->{'BsmlNumbering'};
}

sub returnBsmlNumberingR
{
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    return $self->{'BsmlNumbering'};
}

sub dropBsmlNumbering
{
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    $self->{'BsmlNumbering'} = undef;
}


sub write  {

      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    my $writer = shift;
    
    $writer->startTag( "Sequence", %{$self->{'attr'}} );
    

    foreach my $bsmlattr (sort (keys( %{$self->{ 'BsmlAttr'}})))
    {
	foreach my $value (@{$self->{'BsmlAttr'}->{$bsmlattr}}){

	    $writer->startTag( "Attribute", 'name' => $bsmlattr, 'content' => $value );
	    $writer->endTag( "Attribute" );
	}
    }

    if( my $tcount = @{$self->{'BsmlFeatureTables'}} > 0 ) {

	$writer->startTag( "Feature-tables" );

	foreach my $tbl ( @{$self->{'BsmlFeatureTables'}} ) {

	    $tbl->write( $writer );
	}
	
	if( my $gcount = @{$self->{'BsmlFeatureGroups'}} > 0 )  {
	    
	    foreach my $grp ( @{$self->{'BsmlFeatureGroups'}} ){
		$grp->write( $writer );
	    }
	}
	    
	$writer->endTag( "Feature-tables" );
    }

    # either imbedded or linked sequence data is expected, not both
    
    if( $self->{'BsmlSeqData'} )  {
	$writer->startTag( "Seq-data" );
	$writer->characters( $self->{'BsmlSeqData'} );
	$writer->endTag( "Seq-data" );
    }
    else {
	
	if( $self->{'BsmlSeqDataImport'}->{'source'} ) {
	    $writer->startTag( 'Seq-data-import', %{$self->{'BsmlSeqDataImport'}} );
	    $writer->endTag( 'Seq-data-import' );
	}
    }

    if( $self->{'BsmlNumbering'} ) {

	$self->{'BsmlNumbering'}->write( $writer );
    }

    foreach my $xref (@{$self->{'BsmlCrossReference'}}){
	$xref->write( $writer );
    }
    
      foreach my $link (@{$self->{'BsmlLink'}})  {
	  $writer->startTag( "Link", %{$link} );
	  $writer->endTag( "Link" );
      }

      foreach my $listref (@{$self->{'BsmlAttributeList'}}){
	  $writer->startTag( "Attribute-list");
	    foreach my $hash ( @{$listref} ){ 
		    $writer->startTag( "Attribute",  'name' => $hash->{'name'}, 'content' => $hash->{'content'} );
		    $writer->endTag( "Attribute" );
		}
	   $writer->endTag( "Attribute-list" );
      }
    
    $writer->endTag( "Sequence" );
    
}



sub detectBsmlFeatureTables {
    
    $logger->debug("") if $logger->is_debug;
    
    my $self = shift;

    if (scalar(@{$self->{'BsmlFeatureTables'}}) > 1 )  {
	return 1;
    }
    else {
	return 0;
    }

}


sub subSequence
  {
      my $seq = shift; 
      my ($start, $stop, $complement) = @_;
      
      my $logger = BSML::Logger::get_logger();

    $logger->debug("Entered subSequence") if $logger->is_debug;

    # Hack to check if the substring is on the reverse complement. 
    # Note some of the clients are not setting
    # the complement attributes correctly. 

    if( ($start > $stop) )
    {
	($start, $stop) = ($stop, $start);
	$complement = 1;
    }

    # get the sequence data if it is already in memory. (inline sequence objects)

    my $seqdat = $seq->returnSeqData();

    # otherwise pull in the data from an external file

    if( !($seqdat) )
    {
	my $seqimpt = $seq->returnBsmlSeqDataImport();
	
	# no external file data either; this sequence has no sequence data

	if ( (!($seqimpt->{'format'})) )
	{
	    $logger->warn("The specified sequence (ID=\", $seq->returnattr('id'), \") has no Seq-data or Seq-data-import");
	    return undef;
	}

	if( $seqimpt->{'format'} eq 'fasta' )
	{
	    $seqdat = parse_multi_fasta( $seqimpt->{'source'}, $seqimpt->{'identifier'} );
	}
	
	# Both interal and external lookups for sequence data failed. 
	# Log and error and return the empty string.

	if( !($seqdat) )
	{
	    $logger->warn("No sequence to retrieve for ource: $seqimpt->{'source'} Id: $seqimpt->{'identifier'}");
	    return undef;
	}
	
	#Store the imported sequence in the object layer so further file
	#io is not necessary
	$seqdat =~ s/\s//g;
	$seq->addBsmlSeqData( $seqdat );

      }

    #return the whole sequence if '-1' is passed as the start coordinate

    if( $start == -1 )
      {
	  if( $complement == '0' ){
	      return $seqdat;}
	  else{
	      return reverse_complement( $seqdat );}
      }

    # pull out the substring of interest and return it

    if( $seqdat )
      {
	  if( ($start) > length($seqdat) )
	  {
	      # Log error condition and return the empty string
	      $logger->logdie("BsmlReader::subSequence() - Error: Out of Bounds Substring Access. SeqId($seq) Start($start) Stop($stop) Complement($complement)\n");
	      return '';
	  }
	  
	  # return 5' to 3' data
	  if( $complement == 0 )
	  {
	      my $rseq = substr( $seqdat, $start, ($stop-$start+1));
	      return $rseq;
	  }
	  # return 5' to 3' on rev complement
	  else
	  {
	      return reverse_complement(substr( $seqdat, $start, ($stop-$start+1)));
	  }
      }
}

sub parse_multi_fasta {
  my $fasta_file = shift;
  my $specified_header = shift;

  # extract the directory from the filepath
  my $pos = -1;
  my $posMax = -1;

  while (($pos = index( $fasta_file, "\/", $pos)) > -1 )
  {
      $posMax = $pos;
      $pos++;
  }

  my $fasta_dir = substr( $fasta_file, 0, $posMax ); 
  
  # try to load the sequence using a pregenerated index through. BSML::FastaIndex

  my $indexer = BSML::Indexer::Fasta->new($fasta_file, $fasta_dir);
  
  # Check the health of the various indices for the data file.
  my @check = $indexer->check_indices;

  # Create the indices if necessary...
  if ($check[0] == 1) { $indexer->index_entries };
  if ($check[1] == 1) { $indexer->index_headers };

  # Get the name of the header index file.
  my $header_index = $indexer->header_index;
  my $entry_index = $indexer->entry_index;
 



  my (%h, %e);
  my $c = tie %h, 'CDB_File', $header_index or die "tie failed: $!\n";
  $c = tie %e, 'CDB_File', $entry_index or die "tie failed: $!\n";

  my $entry = $h{$specified_header};

  my ($offset, $length) = split( ',', $e{$entry} );

  # if the sequence could not be loaded with an index, perform a grep operation on the
  # FASTA file.

  open (IN, $fasta_file) or die "Unable to open $fasta_file due to $!";

  if( $offset )
  {
      seek( IN, $offset, 0 );
  }

  my $line = <IN>;
  my $seq_ref = [];
  while(defined($line)) {
    unless($line =~ /^>([^\s]+)/) {
      $line = <IN>;
    } else {
      my $header = $1;

      if($specified_header eq $header) {
	while(defined($line=<IN>) and $line !~ /^>/ ) {
	  next if($line =~/^\s+$/);                   #skip blank lines
	  chomp($line);
	  push(@$seq_ref, $line);
	}
	last;   #seq found, terminating fasta_file parasing
      } else { $line = <IN>; };  #wrong seq, keep looking
    }
    }
  close IN;
  
  my $final_seq = join("", @$seq_ref);
  
  return $final_seq; 
}

sub reverse_complement {

    my $seq = shift;

    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    $seq =~ tr/ATGCatgc/TACGtacg/;
    $seq = reverse($seq);
    return $seq;
}

1;

