#!/usr/local/bin/perl

use strict;
use warnings;
use Carp;

use CAMERA::Parser::BTAB; 
use CAMERA::Parser::HTAB; 
use CAMERA::Parser::ECTable; 
use CAMERA::Parser::TMHMMBSML; 
use CAMERA::Parser::LipoproteinMotifBSML; 
use CAMERA::Parser::Hypothetical;
use CAMERA::PolypeptideSet;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../lib";

$| = 1;

print "\@INC=@INC\n";

my %supported_type = (
                        "BTAB" => 1,
                        "HTAB" => 1,
                        "ECTable" => 1,
                        "TMHMMBSML" => 1,
                        "LipoproteinMotifBSML" => 1,
                        "Hypothetical" => 1,
                        "FRAG_HTAB" => 1,
                     );

my %options = ();
GetOptions (\%options,
            'input_file|i=s',
            'input_type|t=s',
            'output_file|o=s',
            'work_dir|d=s',
            'help|h') || pod2usage();

unless ($options{'input_file'}) {
    confess "must provide a value for flag --input_file";
}
unless ($options{'input_type'}) {
    confess "must provide a value for flag --input_type";
}
unless ($options{'output_file'}) {
    confess "must provide a value for flag --output_file";
}
unless ($options{'work_dir'}) {
    confess "must provide a value for flag --work_dir";
}

unless ($supported_type{$options{'input_type'}}) {
    confess "Input type '$options{input_type}' is not supported";
}
my $input_type = $options{'input_type'};

unless (-e $options{'input_file'}) {
    confess "input file '$options{input_file}' doesn't exist";
}
unless (-d $options{'work_dir'}) {
    confess "work dir '$options{work_dir}' doesn't exist";
}

## create a container for the polypeptides we're annotating
my $pep_set_ref = new CAMERA::PolypeptideSet('HASH' => {});

#print "Check1\n";

open (my $outfh, ">".$options{'output_file'}) || die "couldn't open output file '$options{output_file}' for writing: $!";

#print "Check2\n";


### adding code to handle the input_type FRAG_HTAB it will the HTAB parser but will pass the constructor with one extra param(nikhat)
my $frag =0;
if($input_type =~ /FRAG_HTAB/){
    $frag =1;
    $input_type = "HTAB";
}

#print "Check3\n";

my $parser;
my $parser_type = "CAMERA::Parser::$input_type";
$parser = new $parser_type('polypeptides' => $pep_set_ref, 'work_dir' => $options{'work_dir'},'frag' => $frag);

#print "Check4\n";
my $result = $parser->parse($options{'input_file'}, 1, $outfh);
#print "Check5\n";

if ($result) {
    print $result;
}
