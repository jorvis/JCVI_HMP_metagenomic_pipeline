package BSML::GenePredictionBsml;

=head1 NAME

    GenePredictionBsml.pm - Bsml document creation module specific to gene prediction algorithms

=head1 SYNOPSIS

my $geneBsmlDoc = new BSML::GenePredictionBsml($algorithm, $sourcename, $programversion);

$geneBsmlDoc->addGene($gene);                #Gene.pm object
$geneBsmlDoc->setFasta($gene->getSeqId, $fastaFile);  #Sequence id the gene is located on
                                             #and the fasta file name.

$geneBsmlDoc->writeBsml($outputFile);
    

=head1 USAGE

    This module facilitates the create of a gene prediction BSML document and uses BSML::BsmlBuilder.
    
    It is suggested (although not enforced) for each gene to have a transcript, CDS, polypeptide,
    and exons.  Each gene will be represented as a feature element in the resulting Bsml and each feature
    that is associated with the gene will also be represented as a feature.  If the groups were created
    for the gene, this will be represented as Feature-groups and each feature related to a certain group
    will be come a member of that group.  By using different groups with the same gene we can represent
    differentially expressed genes.  

=head1 CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use warnings;
use BSML::BsmlBuilder;
use Data::Dumper;
use Chado::Gene;

=item my $genePredBsml = new BSML::GenePredictionBsml($algorithm,$sourcename,$programversion);

B<Description:> Creates a new Gene Prediction Bsml builder.

B<Parameters:> $alrogithm - The method of analysis (ex. fgenesh, genscan)
               $sourcename - The path to the source directory of the analysis (component)

B<Returns:> A GenePredictionBsml object

=cut

sub new {
    my $class = shift;
    my $self = { doc       => new BSML::BsmlBuilder(),
                 seq       => [],        #Will hold the sequence id's.
                 algorithm => shift,     #For the analysis (ex. fgenesh, augustus)
                 seq_files => {},        #holds the paths of the files parsed for seq info
                 sourcename => shift,    #for the analysis
                 programversion => shift,  #for the analysis
               };


    bless $self, $class;
    return $self;
}

=item $genePredBsml->addGene($gene);

B<Description:> Adds a Gene object to the document.

B<Parameters:> $gene - a Gene.pm object (see Gene.pm for usage)

B<Returns:> Zero if the gene was not added, otherwise the gene id.

=cut

sub addGene {
    my ($self, $gene, $score) = @_;
    my $retval = 0;

    if($gene->{'seq'}) {
        $self->{'genes'}->{$gene->{'id'}}->{'gene'} = $gene;
        $retval = $gene->{'id'};
    } else {
        warn "$gene->{id} was not added to the bsml document because it is not".
            " located on a sequence.";
    }
    
    $self->{'genes'}->{$gene->{'id'}}->{'score'} = $score
        if($score);

    return $retval;
}

=item $genePredBsml->removeGene($geneId);

B<Description:> Removes a gene from the document.

B<Parameters:> $geneId - The id of the gene to be removed.

B<Returns:> Zero if the gene was not removed, otherwise the gene id.

=cut

sub removeGene {
    my ($self, $geneID) = @_;
    my $retval = $geneID;
    delete($self->{'genes'}->{$geneID}) or
        $retval = 0;
    return $retval;
}

=item $genePredBsml->addSequence($seqId, $fastaFile);

B<Description:> Adds a Sequence to the document.  If parameter
    $seqId is null (or empty), will add every sequence present
    in the fasta file.

B<Parameters:> 

B<Returns:> 

=cut

sub addSequence {
    my ($self, $seqId, $fastaFile) = @_;

    ## has this file been parsed?
    if (! defined $self->{seq_files}->{$fastaFile} ) {
        ## parse the entire file, storing info on each seq
        
        ## support input fasta files that have been gzipped
        if (! -e $fastaFile && -e $fastaFile.".gz") {
            $fastaFile .= ".gz";
        }

        my $infh;

        if ($fastaFile =~ /\.gz$/) {
            open($infh, "<:gzip", $fastaFile)
              || warn("Could not open $fastaFile in addSequence ($!)");
        } else {
            #Open the fasta file.
            open($infh, "< $fastaFile") or
            warn("Could not open $fastaFile in addSequence ($!)");
        }

        ## strip off .gz extension (if present) for the seq-data-import
        $fastaFile =~ s/\.gz//;

        while(my $line = <$infh>) {
            if($line =~ /^>((\S+).*)/) {
                my $ident = $2;
                $self->{'seqs'}->{$ident}->{'fasta'} = $fastaFile;
                $self->{'seqs'}->{$ident}->{'defline'} = $1;
            }
        }
        
        ## just save that it was parsed
        $self->{seq_files}->{$fastaFile} = 1;
    }

    ## check whether it has been found already and fail if not
    if ( ! $seqId || defined $self->{'seqs'}->{$seqId} ) {
        return 1;
    } else {
        die( "failed to find $seqId in file $fastaFile" );
    }
}


=item $genePredBsml->setFasta($seqId, $fastaFile);

B<Description:> See addSequence subroutine.  This 
    function calls that.
    $genePredBsml->addSequence($seqId, $fastaFile);

B<Parameters:> $seqId - id of the sequence. 
               $fastaFile - The file to link to the sequence

B<Returns:> Zero if it could not open the fasta file.

=cut

sub setFasta {
    my ($self, $seqId, $fastaFile) = @_;
    return $self->addSequence($seqId, $fastaFile);
}

=item $genePredBsml->writeBsml($outputFile);

B<Description:> Writes the bsml to the specified file.

B<Parameters:> $seqId - Ouput file name.

B<Returns:> Nothing.

=cut

sub writeBsml {
    my ($self, $outFile, $dtd, $gzip) = @_;

    #Add sequences to the document
    $self->_writeSeqsToDoc();

    #Add genes
    $self->_writeGenesToDoc();

    #Finish Doc
    $self->_finishDoc();

    #And write it all
    $self->{'doc'}->write($outFile, $dtd, $gzip);
    
}

###################PRIVATE METHODS##############################
#or they would be if perl could do that...
sub _writeSeqsToDoc {
    my $self = shift;
    my $doc = $self->{'doc'};

    foreach my $seqId (keys %{$self->{'seqs'}}) {
        unless ( exists($self->{'seqs'}->{$seqId}->{'elem'}) ) {
            $self->{'seqs'}->{$seqId}->{'elem'} = 
                $doc->createAndAddSequence($seqId, undef, undef, 'dna', 'assembly');
            $self->{'seqs'}->{$seqId}->{'elem'}->addBsmlLink('analysis', 
                                                             '#'.$self->{'algorithm'}.'_analysis', 
                                                             'input_of');
            if( exists($self->{'seqs'}->{$seqId}->{'fasta'}) ) {
                $doc->createAndAddSeqDataImport($self->{'seqs'}->{$seqId}->{'elem'}, 'fasta',
                                                $self->{'seqs'}->{$seqId}->{'fasta'}, '', $seqId);
            }

            $doc->createAndAddBsmlAttribute($self->{'seqs'}->{$seqId}->{'elem'}, 'defline',
                                            $self->{'seqs'}->{$seqId}->{'defline'} )
                if(exists($self->{'seqs'}->{$seqId}->{'defline'}));
        }

    }
}

sub _writeGenesToDoc {
    my $self = shift;
    my $doc = $self->{'doc'};

    foreach my $geneId (keys %{$self->{'genes'}}) {
        my $gene = $self->{'genes'}->{$geneId}->{'gene'};

        #If the sequence this gene is on was already created then don't recreate it.
        unless(exists($self->{'seqs'}->{$gene->getSeqId}->{'elem'})) {
            $self->{'seqs'}->{$gene->getSeqId}->{'elem'} = 
                $doc->createAndAddSequence($gene->getSeqId, undef, undef, 'dna', 'assembly');
            $self->{'seqs'}->{$gene->getSeqId}->{'elem'}->addBsmlLink('analysis', 
                                                            '#'.$self->{'algorithm'}.'_analysis', 
                                                            'input_of');
            $doc->createAndAddSeqDataImport($self->{'seqs'}->{$gene->getSeqId}->{'elem'}, 'fasta',
                                            $self->{'seqs'}->{$gene->getSeqId}->{'fasta'}, '', $gene->getSeqId)
                if(exists($self->{'seqs'}->{$gene->getSeqId}->{'fasta'}));

            $doc->createAndAddBsmlAttribute($self->{'seqs'}->{$gene->getSeqId}->{'elem'}, 'defline',
                                            $self->{'seqs'}->{$gene->getSeqId}->{'defline'} )
                if(exists($self->{'seqs'}->{$gene->getSeqId}->{'defline'}));
            
        }

        $self->{'featureTables'}->{$gene->getSeqId}->{'elem'} = 
            $doc->createAndAddFeatureTable($self->{'seqs'}->{$gene->getSeqId}->{'elem'})
            unless(exists($self->{'featureTables'}->{$gene->getSeqId}->{'elem'}));

        my $seqId = $gene->getSeqId;
        my $ft = $self->{'featureTables'}->{$seqId}->{'elem'};
        my $geneFeat = $doc->createAndAddFeature($ft, $gene->getId, '', 'gene');

        #Make sure that the start is always less than the stop
        ($gene->{'start'},$gene->{'stop'}) = ($gene->{'stop'},$gene->{'start'}) if( $gene->{'start'} > $gene->{'stop'} ); 

        $doc->createAndAddIntervalLoc($geneFeat, $gene->{'start'}, $gene->{'stop'}, $gene->{'strand'});
        $geneFeat->addBsmlLink('analysis', '#'.$self->{'algorithm'}."_analysis", 'computed_by');
                
        
        foreach my $groupId (keys %{$gene->{'groups'}}) {
            
            my $features = $gene->getGroup($groupId);
            my $fg = $doc->createAndAddFeatureGroup( $self->{'seqs'}->{$gene->getSeqId}->{'elem'},
                                                     '', $gene->getId);
            $fg->addBsmlFeatureGroupMember( $gene->getId, 'gene' );
            
            foreach my $featObj (@{$features}) {
                my $feat = $doc->createAndAddFeature($ft, $featObj->id, '', $featObj->type);
                $feat->addBsmlLink('analysis', '#'.$self->{'algorithm'}."_analysis", 'computed_by');

                #Make sure the start is less than the stop
                if( $featObj->start > $featObj->stop ) {
                    my $tmp = $featObj->start;
                    $featObj->start($featObj->stop);
                    $featObj->stop($tmp);
                }

                $feat->addBsmlIntervalLoc($featObj->start,
                                          $featObj->stop,
                                          $featObj->strand);
                if($featObj->attribute) {
                    foreach my $scoreType( keys %{$featObj->attribute} ) {
                        $doc->createAndAddBsmlAttribute( 
                                                        $feat, 
                                                        $scoreType, 
                                                        $featObj->attribute->{$scoreType},
                                                       );
                    }
                }

                $fg->addBsmlFeatureGroupMember( $featObj->id, $featObj->type );
                
            }
        }
    }
}

sub _finishDoc {
    my $self = shift;
    
    #Find the sourcename if it hasn't been passed in already
    unless ($self->{'sourcename'}) {

        my $fastaFile = $self->{'seqs'}->{(keys %{$self->{'seqs'}})[0]}->{'fasta'};
        $self->{'sourcename'} = $1 if($fastaFile =~ m|^(.*/\d+_[^/]+)/|);
        warn("Please pass sourcename in at constructor.  We'll let it slide for now".
             " but only if this parses, which doesn't mean it's right.");

    } 
    warn("Could not get sourcename") unless ($self->{'sourcename'});

    #Fill version if it hasn't been passed
    $self->{'programversion'} = 'current' unless ($self->{'programversion'});

    #Make sure algorithm and input have been added.
    if($self->{'algorithm'}) {

        #Write analysis
        $self->{'doc'}->createAndAddAnalysis( 'id' => $self->{'algorithm'}."_analysis",
                                              'sourcename' => $self->{'sourcename'},
                                              'program'    => $self->{'algorithm'},
                                              'algorithm'  => $self->{'algorithm'},
                                              'programversion'    => $self->{'programversion'},
                                             );

    } else {
        warn "The analysis section for the bsml doc was not written because".
            " either input or algorithm was not passed to the GenePredictionBsml constructor.".
            " See GenePredictionBsml.pm documentation for usage details.";
    }
    
}
1;
