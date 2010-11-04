#!/usr/local/bin/perl

use lib (@INC, '/usr/local/devel/ANNOTATION/ard/testing_manual/lib/5.8.8');

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use DBI;
use Split_DB_File;
#require "autoAnnotate.data";

use vars qw(%isoType %ber_rank);
use Time::localtime;

$|++;

my %options = ();
my $results = GetOptions (\%options, 
                          'output_panda_db|o=s',
                          'panda_file|f=s',
                          'password_file|p=s',
                          'omni_db|d=s',
                          'debug=s',
                          'help|h') || die("Bad option getting");

################# GLOBALS AND CONSTANTS #################################
my ($acc,$name,$ec_num,$gene_sym,$species,$go_id,$role_id,$exp,$wgp,$cg) = (0..9);


my @parsed;
#########################################################################
my ($h,$m,$s) = (localtime->hour,localtime->min,localtime->sec);
print "Starting $h:$m:$s\n";

#establish database connection
my ($user, $pass) = ('access', 'access');
my $dbh = &connectToDatabase($user, $pass);

my $omniInfo;
$omniInfo = &getOmniInfo;

($h,$m,$s) = (localtime->hour,localtime->min,localtime->sec);



print "Omni done $h:$m:$s\n";
print "Got omni info\n";


&removeOldFiles('parsedPanda');

print "Running the cdbyank\n";
system('cdbyank -l AllGroup.niaa.cidx > /tmp/tmpIndices.txt');

print "Have the indices file\n";
system('cdbyank -F AllGroup.niaa.cidx < /tmp/tmpIndices.txt > /tmp/pandaTmp_header');

print "Done running cdbyank\n";
open(PANDA, "< /tmp/pandaTmp_header") or die("Can't open /tmp/pandaTmp_header");
open(OUT, " > parsedPanda_header") or die("Unable to open parsedPanda_header");
my %index;
print "About to tie index to pandaTest\n";
tie(%index, 'Split_DB_File', 'parsedPanda') or die("Can't tie hash");
print "Tied index to pandaTest";
my $count = 0;



while(<PANDA>) {
    my $header = &parseHeader($_);
    my $pandaAcc = $1 if(/^(\S+)/);    
    my $outVal = "";

    foreach (@{$header}) {
        $outVal .= $_ if($_);
        $outVal .= "\t";
    }    

    my $indVal = tell(OUT)." ".length($outVal);
    
    #adds values to the index
    $index{$header->[$acc]} = $indVal;
    $index{$pandaAcc} = $indVal;
    
    print OUT $outVal;
    
    $count++;

    if($count%1000 == 0) {
        print "\r$count";
    }
}
print "\n$count headers\n";
($h,$m,$s) = (localtime->hour,localtime->min,localtime->sec);
print "Finishing $h:$m:$s\n";

#close store
$dbh->disconnect();
close(OUT);
close(PANDA);
untie(%index);

sub removeOldFiles {
    my $file = shift;
    
    my $num = 0;

    while(-e $file.".$num") {
        $num++;
        unlink($file.".$num");
        print "\tremoved $file.$num\n";
    }

    if(-e $file."_header") {
        unlink($file."_header");
        print "\tremoved ${file}_header\n";
    }
}

sub parseHeader {
    my $header = shift;
    my $best = [];
    my $best_rank = 100;

    #All the accessions need to be queryable
    my $allAccs = [];

    #Panda is deleminited by this ^|^
    my @individualProts = split(/\^\|\^/,$header);

	#loop through header
    foreach my $protHead (@individualProts) {

        my $current = [];
        my $current_rank;
        my $desc;
        
        if($protHead =~ /^>?(\S+)\s+(.*)/) {
            $current->[$acc] = $1;
            $desc = $2;
        } else {
            print "\n\n:|$protHead|:\n\n";
            die("Could not parse description from $protHead");
        }
        
        #Start parsing the headers to find the best annotated header

        #This will match the part of the line that looks something like this:
        #             (exp=0; wgp=0; cg=1; genomes=0; )
        if ( $desc =~ 
             s/ \(exp(erimental)?=(\-?\d+); wgp=(\-?\d+); cg=(\-?\d+);[^\)]+\)//) {
            $current->[$exp] = $2;
            $current->[$wgp] = $3;
            $current->[$cg] = $4;
        }

        #We should also check from the database here I guess
        ### HERE ###
        # checkExperimentalData( $current->{'accession'}) unless($current->{'exp'});

        #If we already have an experimentally characterized match, 
        #don't bother with the ones that aren't.
        next if(defined($best->[$exp]) && 
                $best->[$exp] == 1 && 
                $current->[$exp] != 1);
        
        #Find the species if it looks like this
        #        {Chlorobium chlorochromatii CaD3;}
        $current->[$species] = $1 if ($desc =~ s/\{([^;\}]+);.*\}//);

        #Specific parsing depeding on originating database
        #Swissprot
        if($current->[$acc] =~ /^SP/) {

            #Get the ec_num(s)
            my @ec;
            while ($desc =~ s/\(EC ([^\)]+)\)//) { 
                push @ec, $1;
            }
            $current->[$ec_num] = join " ", @ec;

            # SwissProt has alternate names (synonyms in SP talk) in parens following the EC
            while ($desc =~ /(\s*\([^\)]*\)\s*\.?\s*)$/) { 
                $desc =~ s/\s*\([^\)]*\)\s*\.?\s*$//;
            }
            
            # bifunctional proteins have separate functions enclosed in brackets. kill them.
            $desc =~ s/\[(includes|contains)\:.*\]//i;
            $current->[$name] = $desc;
        #checks for omnium header
        } elsif($current->[$acc] =~ /^OMNI/) {
            $current->[$name] = $desc;
                        
            my $tmpAcc = $1 if($current->[$acc] =~ /^(OMNI\|[^\|]+)|/);
            foreach my $i (0..9) {
            	#if locus is part of the $omniInfo array add to current
                $current->[$i] = $omniInfo->{$tmpAcc}->[$i]
                    if($omniInfo->{$tmpAcc}->[$i]);
                                       
            } 
 			#return $best; 
 
            #This is basically for genbank headers (although I've seen some PDB
            #entries showing up that we don't deal with).  Not sure if they
            # *should* also be handled here, but they currently are.
        } else {
            
            # FRAMESHIFT and the like sometimes appear between the gene
            # symbol and the species. This needs to be removed
            $desc =~ s/(authentic )?frameshift//i;
            $desc =~ s/(authentic )?point mutation//i;

            if ($desc =~ /\(([a-z]{2,4}[A-Z\d\-]{,3})\)/) {
                my $tmpGS = $1;
                $current->[$gene_sym] = $1;
                $desc =~ s/\($tmpGS\)\s*//;
            }
            
            # get rid of trailing stuff in parens
            $desc =~ s/\(.*\)\s*$//;

            #Now this is the common name.  Who knew.
            $current->[$name] = $desc;

        }

        #Check for go (unless it's omnium because they should already have them)
        if($current->[$name] !~ /^OMNI/ && $current->[$ec_num]) {
            my $go_ids = &getGoIdFromEc($current->[$ec_num]);
            $current->[$go_id] = join(" ", @{$go_ids});
        }   
        
        #Give it a rank
        $current_rank = 100;

        &cleanAnnotation($current);                
        if(&swapBestHeader($best, $best_rank, $current, $current_rank)) {
            ($best, $best_rank) = ($current, $current_rank);
        }
    }
    
    return $best;  
}

#returns the best representative among header entries
sub swapBestHeader {
    my ($best, $best_rank, $current, $current_rank) = @_;
    
       
    my $retval = 0;
    
    if (! $best->[$name]) {
        $retval = 1;
    } else {

        if (# hits with rank < 3 are either experimental or highly curated, so info is to be trusted.
            # But we don't want conserved hypo unless necessary, so don't substitute for that.
            ($current_rank < 3 && $current_rank < $best_rank ) ||
            ($best->[$name] =~ /hypothetical|probable|putative|related|unknown/i &&
             $current->[$name] !~ /hypothetical|probable|putative|related|unknown/i) &&
             $current->[$name] !~ /hypothetical protein/ ) {
            
            $retval = 1;
           
        } elsif ($best_rank == $current_rank) {
        	 
            my ($bestc, $testc) = (0,0);
            $bestc++ if ($best->[$gene_sym]);
            $bestc++ if ($best->[$ec_num]);
            $testc++ if ($current->[$gene_sym]);
            $testc++ if ($current->[$ec_num]);
            
            #If there is a better score
            if ($testc > $bestc) {
                $retval = 1;
            }
        }
    }
    return $retval;
}




sub getGoIdFromEc {
    my $ec_nums = shift;
    my @ecs = split(/\s+/, $ec_nums);
    my @go_ids = ();
    foreach my $ec (@ecs) {
        my $query = "select go_id "
            . " FROM common..go_map "
            . " where db = \"EC\" " 
            . " AND identifier = \"$ec\" "
            . " AND go_id !=\"GO:0000004\" "
            . " AND go_id !=\"GO:0005554\"";

        my $sth = $dbh->prepare($query);
        $sth->execute;
        my $go;
        $sth->bind_columns(\$go);
        while($sth->fetch) {
            push(@go_ids, $go);
        }
        $sth->finish;
    }
    return \@go_ids;
}

#loads the omnium data into an associative array
sub getOmniInfo {
	 
    my $retval;
    
    my $query = "select i.locus, i.com_name, i.gene_sym, i.ec_num, t.genus, t.species"
        . " from omnium..nt_ident i, omnium..asm_feature a, omnium..asmbl_data d, omnium..taxon_link tl, omnium..taxon t "
#        . " where i.locus=\"NTL01HC5962\""
# 		. " and i.locus=a.locus"
        . " where i.locus=a.locus"
		. " and a.asmbl_id=d.asmbl_id"
		. " and d.db_data_id=tl.db_taxonl_id"
		. " and tl.taxon_uid=t.uid";
       
    my $sth = $dbh->prepare($query);
    print "\tExecuting the query: [$query]\n";
    $sth->execute();
    print "\tFetching the rows\n";
    my $arrayref;
    while($arrayref = $sth->fetch) {
        my $locus = "OMNI|".$arrayref->[0];
        $retval->{$locus}->[$acc] = $locus;
        $retval->{$locus}->[$name] = $arrayref->[1];
        $retval->{$locus}->[$gene_sym] = $arrayref->[2];
        $retval->{$locus}->[$ec_num] = $arrayref->[3];
        $retval->{$locus}->[$species] = $arrayref->[4]." ".$arrayref->[5];
    }
    $sth->finish();

    $query = "select i.locus, i.com_name, i.gene_sym, i.ec_num, t.genus, t.species"
        . " from omnium..ident i, omnium..asm_feature a, omnium..asmbl_data d, omnium..taxon_link tl, omnium..taxon t "
#        . " where i.locus='NTL01HC5962'"
# 		. " and i.locus=a.locus"
        . " where i.locus=a.locus"
		. " and a.asmbl_id=d.asmbl_id"
		. " and d.db_data_id=tl.db_taxonl_id"
		. " and tl.taxon_uid=t.uid";
    
    $sth = $dbh->prepare($query);
    print "\tExecuting the query: [$query]\n";
    $sth->execute();
    print "\tFetching the rows\n";
    
    while($arrayref = $sth->fetch) {
        my $locus = "OMNI|".$arrayref->[0];
        $retval->{$locus}->[$acc] = $locus;
        $retval->{$locus}->[$name] = $arrayref->[1];
        $retval->{$locus}->[$gene_sym] = $arrayref->[2];
        $retval->{$locus}->[$ec_num] = $arrayref->[3];
        $retval->{$locus}->[$species] = $arrayref->[4]." ".$arrayref->[5];
    }
    $sth->finish();

     $query = "select distinct locus, role_id"
         . " from omnium..role_link"
         . " where role_id != 156"
         . " and role_id != 185"
         . " and role_id != 270";

    $sth = $dbh->prepare($query);
    print "\tExecuting the query: [$query]\n";
    $sth->execute();
    my ($locus,$role);
    $sth->bind_columns(\$locus,\$role);
    while($sth->fetch) {
        my $newLoc = "OMNI|$locus";
        $retval->{$newLoc}->[$role_id] .= $role." ";
        
    }
    $sth->finish();


    $query = "select r.locus, r.go_id "
        . " from omnium..go_role_link r, omnium..go_evidence e"
        . " where r.go_id !=\"GO:0000004\" "
        . " AND r.go_id !=\"GO:0005554\" "
        . " and r.id=e.role_link_id"
        . " and e.ev_code != 'IEA'";
        #. " and r.locus='GSGORF0055'";
    
    $sth = $dbh->prepare($query);
    print "\tExecuting the query: [$query]\n";
    $sth->execute();
    my $go;
    $sth->bind_columns(\$locus,\$go);
    while($sth->fetch) {
        my $tmpLoc= "OMNI|$locus";
        $retval->{$tmpLoc}->[$go_id] .= $go." ";
    }
    $sth->finish();
    print "\tFinished gathering omnium info\n";
	
    return $retval;
}

sub connectToDatabase {
    my ($u, $p) = @_;
    my $retval = DBI->connect('dbi:Sybase:server=SYBTIGR; packetSize=8092', $u, $p)
        or die("Unable to connect to database".${DBI->errstr});
    return $retval;    
}




sub cleanAnnotation {
    my $annot = shift;
    
    # this subroutine should hold ALL error checking for annotation info
    my ($name, $gene_sym, $ec) = (\$annot->[$name], \$annot->[$gene_sym], \$annot->[$ec_num]);

    # this is to try to get rid of accessions that may be part of the name
    $$name =~ s/\b[A-Z]{2,}_?[A-Z]*[\d\.]{3,}[a-z]?\s*//g;
    $$name =~ s/\b[A-Z][a-z]{1,2}\d*\-?\d{3,}[a-z]?//g;
    
    # try to identify names that are in all caps and put them in lowercase
    if ($$name =~ /[A-Z]{4,}/) { $$name =~ tr/A-Z/a-z/ }
    # but uppercase these words
    if ($$name =~ /\b((nadp?h?|atp|v?i+v?|\w))\b/i) { 
        my $s = uc($1);
        $$name =~ s/$1/$s/;
    }
    
    # capitalize the first letter of protein names
    if ($$name =~ /\b([a-z]{3}([A-Z]|\d+))\b/) {
        my $gene = $1;
        my $new = ucfirst($gene);
        $$name =~ s/$gene/$new/;
    }

    # get rid of bad text
    $$name =~ s/[\:\;] *.* similar to .*//i;
    $$name =~ s/\[imported\]//;
    $$name =~ s/\"//g;
    $$name =~ s/\([^\)]*$//;
    $$name =~ s/\[[^\]]*$//;
    $$name =~ s/\{[^\}]*$//;

    # look for non-homology based info
    $$name =~ s/(authentic )?frameshift//i;
    $$name =~ s/(authentic )?point mutation//i;
    $$name =~ s/ (selenoprotein|selenocysteine\-? ?containing)//i;
    $$name =~ s/, ?(truncation|degenerate|pseudogene)//i;
    $$name =~ s/, ?interruption\-[nc]//i;
    $$name =~ s/(probable|putative|predicted|hypothetical) (lipoprotein|transmembrane protein)(, ?putative)?//i;
    $$name =~ s/similar to .*//i;
    $$name =~ s/alternate start at bp \d+\;//i;
    $$name =~ s/\s+/ /g;
    $$name =~ s/^\s+//;
    $$name =~ s/\s+$//;
    $$name =~ s/\s*\.$//;
    $$name =~ s/\bpossible\b/putative/;
    $$name =~ s/precursor//i;
    $$name =~ s/(chloroplast|mitochondrial)//i;

    # look for bad names
    if ($$name =~ /^orf/i ||
        $$name =~ /^\s*protein\s*$/i ||
        $$name =~ /^\s*putative\s*$/i ||
        $$name =~ /hypothetical[0-9\.kda ]* protein/i ||
        $$name =~ /^\s*hypothetical\s*$/i ||
        $$name =~ /^unknown( protein)?/i ||
        $$name =~ /^\s*$/ || 
        $$name =~ /\(\).*/ || 
        $$name =~ /^\s*\d+\s*$/ || 
        $$name =~ /^AGR/ ||
        $$name =~ /^CG(\d+) gene product/ ||
        $$name =~ /^RX/ ||
        $$name =~ /unnamed protein product/i ||
        $$name =~ /No significant matches/i ||
        $$name =~ /conserved domain protein/i ||
        $$name =~ /\bpredicted\b/i ||
        $$name =~ /\bstructure\b/i ||
        $$name =~ /\bCOG\d+/ ||
        $$name =~ /incomplete/i ||
        $$name =~ /\-like/i ||
        $$name =~ /uncharacteri(s|z)ed/i ||
        $$name =~ /\bsimilar\b/i ||
        $$name =~ /\bidentical\b/i
        ) {
        $$name = "conserved hypothetical protein";
    }

    if ($$name =~ /\bhomolog(ue)?\b/ && $$name !~ /putative/) { 
        $$name = "putative ".$$name; 
    }

    # if name is really long for some reason, truncate it
    if (length($$name)>255) {
        my $trunc_name = substr($$name, 0, 255);
        $$name = $trunc_name;
    }

    # Is there an ec hiding in the name?
    if ($$name =~ s/(\d+(\.(\d+|-)){0,2}((\.\d+){3}|\.-))//) {
        if (!$$ec) { 
            $$ec = $1;
        }
    }
    $$name =~ s/(\[|\()ec\:?\s*(\]|\))//i; # this could be a leftover from EC finding

    ## remove leading, trailing spaces:
    $$name =~ s/^ *//;
    $$name =~ s/\,*\.* *$//;

    # try to lowercase first letter
    if ($$name =~ /^[A-Z][a-z]{3,}/) {
        $$name = "\l$$name";
    }

    # if the data doesn't look good, delete it
    if ($$name !~ /\w+/ || $$name =~ /\w{,2}/) { 
        $$name = "";
    }

    if ($$gene_sym) {

        #leading and trailing spaces
        $$gene_sym =~ s/^ *(.*) *$/$1/;

        # TIGR indexing
        $$gene_sym =~ s/\-?\d+$//;

        # screen for bad data
        if ($$gene_sym !~ /\w+/) { 
            $$gene_sym = "";
        }
        if (length($$gene_sym) < 3) { 
            $$gene_sym = "";
        }
    }

    if ($$ec) {
        $$ec =~ s/^ *(.*) *$/$1/;
        if ($$ec =~ /^[^\d\.\-]+$/) {
            $$ec = "";
        }
    }
}
