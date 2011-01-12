#!/usr/local/bin/perl

#use lib "/share/apps/ergatis/jboss-test/lib";
#use lib "/usr/local/devel/ANNOTATION/mg-annotation/testing/smurphy-20090825-DO_NOT_USE_FOR_PRODUCTION/current/annotation_tool";

use strict;
use warnings;
use Carp;
use CAMERA::Polypeptide;
use CAMERA::AnnotationData::Polypeptide;
use CAMERA::AnnotationRules::PredictedProtein;
use Data::Dumper;

use Getopt::Long;

$| = 1;

my %options = ();
GetOptions (\%options,
            'input|i=s',
            'output|o=s',
            'synonyms|s=s',
            'help|h') || pod2usage();

unless($options{'input'}) {
    confess "must specify value for --input";
}
unless($options{'output'}) {
    confess "must specify value for --output";
}
unless (-e $options{'input'}) {
    confess "input file '$options{input}' doesn't exist";
}

open (my $infh, $options{'input'}) || die "Couldn't open input table file '$options{input}' for reading: $!";
open (my $outfh, ">".$options{'output'}) || die "couldn't open output file '$options{output_file}' for writing: $!";
my $rules = new CAMERA::AnnotationRules::PredictedProtein($options{'synonyms'});

my $last_id = '';
my $pep = undef;
while (<$infh>) {
    chomp;
  
    my $annot = new CAMERA::AnnotationData::Polypeptide();
    $annot->from_string($_);

    ## we've hit the first or a new polypeptide id in the table
    if ($annot->{'id'} ne $last_id) {
        if (defined($pep)) {
            $rules->annotate($pep,$outfh);
            $pep = undef;
        }
        
        $pep = new CAMERA::Polypeptide( 'id' => $annot->{'id'} );
    }
    $pep->add_annotation($annot);
    $last_id = $annot->{'id'};
}
## do the last entry if needed
if (defined($pep)) {
    $rules->annotate($pep);
}
