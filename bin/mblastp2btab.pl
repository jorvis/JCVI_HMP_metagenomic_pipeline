#!/usr/bin/perl

################## TO DO
# fix- it's breaking output into 2 lines
# colunm 19:hit length - is this in the file??
# colunm 18:strand - confirm that strand info is always in () in query_id col
# column 17: blast frame - confirm format...currently is in 0/0 format.
# column 12: rounding issue
# TEST
# ################


=head1 NAME

mblastp2btab.pl - convert an mBlastP tab-delimited output file to BTAB format.

=head1 SYNOPSIS

USAGE: mblastp2btab.pl 
            --input=/path/to/input_file.raw
            --output=/path/to/output_file.btab
          [ --log=/path/to/logfile
            --debug=N
          ]

=head1 OPTIONS

B<--input,-i>
    The input mblastp tab-delimited output file from the mblastp suite.

B<--output,-o>
    The file to which the parsed output will be written.

B<--reference_fasta,-r> 
    The FASTA file used as the reference in the BLAST searches.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to parse a tab-delimited output file from the mblastp suite of programs into
btab format. 

=head1  INPUT

The input to this sequence is defined using the --input option.  This should point
to the  mblastp tab-delimited output file.

The --reference_fasta refers to the FASTA file your queries were searched against.  This is
necessary because mBlastP doesn't carry over any annotation in the FASTA descriptors so it is
indexed by this script.  The headers will be hashed in memory so it's better to run this script
on large, concatenated result files.

=head1  OUTPUT

The output is defined using the --output option.  The file created is tab-delimited and
composed of the following columns.

    1   query_name
    2   date
    3   query_length
    4   algorithm
    5   database_name
    6   hit_name
    7   qry_start
    8   qry_end
    9   hit_start
    10  hit_end
    11  percent_identity
    12  percent_similarity
    13  raw_score
    14  bit_score
    15  NULL
    16  hit_description
    17  blast_frame
    18  qry_strand (Plus | Minus)
    19  hit_length
    20  e_value
    21  p_value

=head1 TESTING

cd /usr/local/projects/dacc/jcvi_metagenomic_autoannotate/output/PGA/SRS011140

perl -I ~/svn/ergatis/lib/ /home/jorvis/git/JCVI_HMP_metagenomic_pipeline/bin/mblastp2btab.pl -i SRS011140.mblastp.custom2.out -o SRS011140.mblastp.custom2.btab -r /usr/local/projects/dacc/jcvi_metagenomic_autoannotate/data/uniref100.fasta

/usr/local/projects/dacc/jcvi_metagenomic_autoannotate/bin/camera_parse_annotation_results_to_text_table.pl --input_file=SRS011140.mblastp.custom2.btab --input_type=BTAB --output_file=SRS011140.mblastp.custom2.btab.parsed --work_dir=/usr/local/projects/dacc/jcvi_metagenomic_autoannotate/data/ 

=head1  CONTACT

    Kemi Abolude
    kabolude@som.umaryland.edu

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use bigint;
use Ergatis::Logger;
use Bio::SearchIO;

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
                          'reference_fasta|r=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## holds the identifiers from the reference FASTA file as keys and the rest of each
#   header line as the vlue.
my %annotations = ();

open(my $fasta_fh, $options{reference_fasta}) || die "failed to read reference FASTA file: $!";

while (<$fasta_fh>) {
    $annotations{$1} = $2 if ( /\>(\S+) (.+)/ );
}


## parse database name
my $database_name = "";
# open the input file for database name
open(IN, "<$options{input}") || $logger->logdie("can't read the input sequence: $!");
while (my $line = <IN>){
    next unless ($line =~ /Database/i);
	my @db_name = split /:/, $line;
	$database_name = $db_name[1];
	$database_name =~ s/^\s+//;
}
close IN;


## open the input file for parsing                                                                                                                                                                    
open(IN, "<$options{input}") || $logger->logdie("can't read the input sequence: $!");

## open the output file:                                                                                                                                                                              
open (OUT, ">$options{output}") || $logger->logdie("can't create output file for BLAST report: $!");


   # parse each line:
   while( <IN> ) {
      
       chomp;
      
       #skip mblastp header or footer
       next if ( /^Query Id\s+Reference Id/i ||
		 /^\s*Number of Sequences/i ||
		 /^\s*Reference Database/i);

       #skip whitespace
       next if ( /^\s*$/ );

       my @cols = split /\t/;
       my @x;

       ## there should be 12 elements in cols, if not, we have an unrecognized format.
       unless (scalar @cols == 12) {
            $logger->error("the following mblastp line was not recognized and could not be parsed:\n$_\n") 
                if ($logger->is_error);
            next;
        }

       #1: query_name
       my @q_name = split /\s+/,$cols[0];
       $x[0] = $q_name[0];

       #2: date
       ##info not in file, today's date for now
       my @time = localtime(time);
       $x[1] = join('/', $time[4]+1,$time[3],$time[5]+1900,);  
      
       #3: query_length
       $x[2] = $cols[9]-$cols[8]+1;  #double-check

       #4: algorithm
       $x[3] = "mblastp";

       #5: database_name
       $x[4] = $database_name;
       #database name will get parsed with whitespace if its path is long
       $x[4] =~ s/\s//g;
            
       #6: hit_name
       $x[5] = $cols[1];

       #7: qry_start
       $x[6] = $cols[8];

       #8: qry_end
       $x[7] = $cols[9];

       #9: hit_start
       $x[8] = $cols[10];
       
       #10: hit_end
       $x[9] = $cols[11];

       #11: percent_identity
       $x[10] = sprintf ("%.1f", $cols[4]);

       #12: percent_similarity
       my $similarity = ($cols[6] * 100)/ $cols[5];
       $x[11] = sprintf ("%.1f", $similarity);

       #13: raw_score
  
       #14: bit_score
       $x[13] = $cols[3];
  
       #15: NULL
       
       #16: hit_description
       $x[15] = $annotations{$cols[1]};
       
       #17: blast_frame
       $x[16] = $cols[7];
	   #  $x[16] = ( ($hsp->query->frame + 1) * $hsp->query->strand); #blast frame (1, 2, 3, -1, -2, -3).

       #18: qry_strand
       #assuming that strand information will always be last item in () Query_id colunm
       my $strand_descript = "null";
       my $query_strand =  substr $cols[0], -2, 1;
         if ($query_strand eq '+') {
             $strand_descript = "Plus";
         } elsif ($query_strand eq '-') {
             $strand_descript = "Minus";
         }
       $x[17] = $strand_descript;

       #19: hit_length
       #doesn't seem like there's a way to get this

       #20: e_value
       $x[19] = $cols[2];
       
       #21: p_value
       my $pvalue = &calculate_pvalue( $cols[2] );
       $x[20] = $pvalue;

       my $outline = join ("\t", @x);
       print OUT "$outline\n";
   }
close IN;
close OUT;

#See http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
sub calculate_pvalue {
    my $evalue = shift;

    my $estimate = 0.57721566490153;
    
    #my $p = 1 - (bexp( (-1*$evalue), 4 ) );
    my $pvalue = ( 1 - ( $estimate**(-1*$evalue) ) );
    if ($pvalue eq 'NaN'){
	return 0;
    }else{
	return $pvalue;
    }
}

sub check_parameters {
    my $options = shift;
    
    unless (-e $options{input}) {
        $logger->logdie("input option not passed or does not exist!");
        exit(1);
    }

    unless (defined $options{output}) {
        $logger->logdie("output option not passed");
        exit(1);
    }
    
    unless (defined $options{reference_fasta}) {
        $logger->logdie("reference_fasta option not passed");
        exit(1);
    }
}

