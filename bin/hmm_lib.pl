#!/usr/local/bin/perl

sub prot_acc_link {
    use strict;
    my ( $acc, $mode, $out_r ) = @_;
    $mode = 1
      if ( $mode !~ /^\s*\d+\s*$/
        or $mode <= 0 );
    if ( $mode == 1 ) {
        if ( $acc =~ /^EGAD\|(\d+)/ ) {
            $$out_r =
"<a target=_blank href=\"http://wild/docs/tigr-scripts/egad_scripts/ht_report.spl?prot_id=$1\">$acc</a>";
        }
        elsif ( $acc =~ /^SP\|([^\|\/\:\s]+)/ ) {
            $$out_r =
"<a target=_blank href=\"http://ca.expasy.org/cgi-bin/niceprot.pl?$1\">$acc</a>";

#	    $$out_r = "<a target=_blank href=\"http://www.expasy.ch/cgi-bin/sprot-search-ac?$1\">$acc</a>";
        }
        elsif ( $acc =~ /^GP\|([^\|\/\:\s]+)/ ) {
            $$out_r =
"<a target=_blank href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Protein&cmd=Retrieve&dopt=GenPept&list_uids=$1\">$acc</a>";
        }
        elsif ( $acc =~ /^(?:GB|RF)\|([^\|\/\:\s]+)/ )
        {    # GB = Genbank,  RF = Refseq.
            $$out_r =
"<a target=_blank href=\"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=$1\">$acc</a>";
        }
        elsif ( $acc =~ /^PIR\|([^\|\/\:\s]+)/ ) {
            $$out_r =
"<a target=_blank href=\"http://www-nbrf.georgetown.edu/cgi-bin/pirwww/nbrfget?id=$1&xref=0&fmt=C\">$acc</a>";
        }
        elsif ( $acc =~ /^OMNI\|([^\|\/\:\s]+)/ ) {
            $$out_r =
"<a target=_blank href=\"http://www.tigr.org/tigr-scripts/CMR2/GenePage.spl?locus=$1\">$acc</a>";
        }
        elsif ( $acc =~ /^gi\|(\d+)/i ) {
            $$out_r =
"<a target=_blank href=\"http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=Protein&list_uids=$1&dopt=GenPept\">$acc</a>\n";
        }
        elsif ( $acc =~ /^PDB\|(\d\S*)\_\S*/i ) {
            $$out_r =
"<a target=_blank href=\"http://www.rcsb.org/pdb/cgi/explore.cgi?pdbId=$1\">$acc</a>\n";
        }
        else {
            $$out_r = "$acc";
        }
    }
}

sub filter_seq_by_cutoff {

#delete a sequence if its score is below cutoff
#this subroutine need accurate lend and rend saved in sequence and header info in STOCKHOLM format
    use strict;
    my ( $seq_r, $cutoff, $cutoff2 ) = @_;
    my %tmp;
    for my $n ( 0 .. $$seq_r{number} - 1 ) {
        my $acc = (
            $$seq_r{$n}{db} eq "EGAD"
            ? "$$seq_r{$n}{db}|$$seq_r{$n}{prot_id}"
            : "$$seq_r{$n}{db}|$$seq_r{$n}{0}"
        );
        push @{ $tmp{$acc} },
          {
            "seq_num" => $n,
            "lend"    => $$seq_r{$n}{lend},
            "rend"    => $$seq_r{$n}{rend}
          };
    }
    for my $acc ( keys %tmp ) {
        next
          if ( scalar( @{ $tmp{$acc} } ) == 1 );
        my @tmp1 = sort { $$a{lend} <=> $$b{lend} } @{ $tmp{$acc} };
        for my $n ( 0 .. @tmp1 - 1 ) {
            $$seq_r{ $tmp1[$n]{seq_num} }{domain_num} = $n + 1;
        }
    }
    for my $n ( 0 .. $$seq_r{number} - 1 ) {
        my $acc = (
            $$seq_r{$n}{db} eq "EGAD"
            ? "$$seq_r{$n}{db}|$$seq_r{$n}{prot_id}"
            : "$$seq_r{$n}{db}|$$seq_r{$n}{0}"
        );
        if (   $$seq_r{stockholm}{GF}{score}{$acc}{hmm} ne ''
            && $$seq_r{stockholm}{GF}{score}{$acc}{b} < $cutoff
            and $cutoff =~ /\d+/ )
        {
            delete $$seq_r{$n};
            next;
        }
        my $domain_num = $$seq_r{$n}{domain_num};
        if (    $$seq_r{stockholm}{GF}{score}{$acc}{hmm} ne ''
            and $domain_num <= 0
            and $$seq_r{stockholm}{GF}{score}{$acc}{b} < $cutoff2
            and $cutoff2 =~ /\d+/ )
        {
            delete $$seq_r{$n};
            next;
        }
        if (    $$seq_r{stockholm}{GF}{score}{$acc}{hmm} ne ''
            and $domain_num > 0
            and $$seq_r{stockholm}{GF}{score}{$acc}{$domain_num}{b} < $cutoff2
            and $cutoff2 =~ /\d+/ )
        {
            delete $$seq_r{$n};
            next;
        }
    }
    &shrink_seq($seq_r);
    &shrink_alignment($seq_r);
}

sub truncate_file {

    #truncate a file by keeping a number of lines on top or end of the file
    #if $place eq "end", the last $max_line lines will be keeped in file.
    use strict;
    my ( $file, $max_line, $place, $log_header, $log_r ) = @_;
    $log_header .= "--sub truncate_file";
    if ( -e $file ) {
        system "chmod 760 $file"
          if ( !-w $file or !-r $file );
        if (    -w $file
            and -r $file )
        {
            open( FH, "wc $file |" );
            my $tmp = <FH>;
            close FH;
            my $line = $1
              if ( $tmp =~ /(\d+)/ );
            return "success"
              if ( $line < $max_line );
            my $dir       = '';
            my $file_name = $file;
            if ( $file =~ /^(.+)\/([^\/]+)$/ ) {
                $dir       = $1;
                $file_name = $2;
            }
            my $tmp_file = &tmp_file( $dir, $file_name );
            if ( $place eq "top" ) {
                system "head -$max_line $file > $dir/$tmp_file";
                system "mv $dir/$tmp_file $file";
            }
            elsif ( $place eq "end" ) {
                system "tail -$max_line $file > $dir/$tmp_file";
                system "mv $dir/$tmp_file $file";
            }
            else {
                $$log_r .=
"***Error: $log_header--unknoen place $place, it should be top or end\n"
                  if ($log_r);
                return "error";
            }
        }
        else {
            $$log_r .=
"***Error: $log_header--do not have permission to change the file $file\n"
              if ($log_r);
            return "error";
        }
    }
    else {
        $$log_r .= "***Error: $log_header--can not find file $file\n"
          if ($log_r);
        return "error";
    }
    return "success";
}

sub print_cyc {
    use strict;
    my ( $file, $string, $overwrite, $to_STDOUT ) = @_;
    (
        $overwrite == 1
        ? open( CYC, ">$file" )
        : open( CYC, ">>$file" )
    );
    print CYC $string;
    print $string
      if ( $to_STDOUT == 1 );
    close CYC;
}

sub modify {
    use strict;
    my ( $hmm_acc, $hmm_name, $hmm_com_name, $hmm_r ) = @_;
    $$hmm_r =~ s/\nNAME\s+\S+\n/\nNAME  $hmm_acc\n/;
    if ( $$hmm_r =~ /\nDESC(.*)\n/ ) {
        $$hmm_r =~ s/\nDESC(.*)\n/\nDESC  $hmm_name: $hmm_com_name\n/;
    }
    else {
        $$hmm_r =~ s/\nLENG  \d+\n/\nDESC  $hmm_name: $hmm_com_name$&/;
    }
}

sub cat_hmm_lib {

    # cat new hmms to lib files
    use strict;
    my ( $hmm_new_r, $dir_ind, $dir_lib ) = @_;
    system "chmod 660 $dir_lib/*";

# set hmmconvert command - as of LDhmmer-1.5.3, hmmconvert by default converts models to binary format
# Because of this they also deprecated the -b option, and in fact the program fails if you specify -b.
# But on Solaris machines, we're using regular hmmer, and thus regular hmmconvert so we need the -b.
# Thus:
    my $arch = `uname -s`;
    chomp $arch;
    my $hmmconvert;
    if ( $arch != /linux/i ) {
        die "Your operating system is impure: It does not support LDhmmer\n";
    }
    $hmmconvert = "hmmconvert";
    for my $hmm_acc ( keys %$hmm_new_r ) {
        if ( -s "$dir_ind/$hmm_acc.HMM" ) {
            system
"$hmmconvert $dir_ind/$hmm_acc.HMM  $dir_ind/${hmm_acc}_bin.HMM  > $dir_ind/${hmm_acc}_bin.HMM.log";
            system "cat $dir_ind/$hmm_acc.HMM       >> $dir_lib/ALL_LIB.HMM";
            system
"$hmmconvert -A $dir_ind/${hmm_acc}_bin.HMM $dir_lib/ALL_LIB_bin.HMM";

     #	    system "cat $dir_ind/${hmm_acc}_bin.HMM >> $dir_lib/ALL_LIB_bin.HMM";
            system "cat $dir_ind/$hmm_acc.HMM       >> $dir_lib/TOGA_LIB.HMM"
              if ( $hmm_acc =~ /^TIGR/ );
            system
"$hmmconvert -A $dir_ind/${hmm_acc}_bin.HMM $dir_lib/TOGA_LIB_bin.HMM"
              if ( $hmm_acc =~ /^TIGR/ );

#	    system "cat $dir_ind/${hmm_acc}_bin.HMM >> $dir_lib/TOGA_LIB_bin.HMM" if($hmm_acc =~ /^TIGR/);
            system "cat $dir_ind/$hmm_acc.HMM       >> $dir_lib/PFAM_LIB.HMM"
              if ( $hmm_acc =~ /^PF/ );
            system
"$hmmconvert -A $dir_ind/${hmm_acc}_bin.HMM $dir_lib/PFAM_LIB_bin.HMM";

#	    system "cat $dir_ind/${hmm_acc}_bin.HMM >> $dir_lib/PFAM_LIB_bin.HMM" if($hmm_acc =~ /^PF/);
            unlink(
                "$dir_ind/${hmm_acc}_bin.HMM.log",
                "$dir_ind/${hmm_acc}_bin.HMM"
            );
        }
        if ( -s "$dir_ind/$hmm_acc.FRAG" ) {
            system
"$hmmconvert $dir_ind/$hmm_acc.FRAG $dir_ind/${hmm_acc}_bin.FRAG > $dir_ind/${hmm_acc}_bin.FRAG.log";
            system "cat $dir_ind/$hmm_acc.FRAG       >> $dir_lib/ALL_LIB.FRAG";

            # fix 1 for Logical Depth binary libraries constructions
            system
"$hmmconvert -A $dir_ind/${hmm_acc}_bin.FRAG $dir_lib/ALL_LIB_bin.FRAG";

  #             "cat $dir_ind/${hmm_acc}_bin.FRAG >> $dir_lib/ALL_LIB_bin.FRAG";
            system "cat $dir_ind/$hmm_acc.FRAG       >> $dir_lib/TOGA_LIB.FRAG"
              if ( $hmm_acc =~ /^TIGR/ );

            # fix 2
            system
"$hmmconvert -A $dir_ind/${hmm_acc}_bin.FRAG $dir_lib/TOGA_LIB_bin.FRAG"

  #             "cat $dir_ind/${hmm_acc}_bin.FRAG >> $dir_lib/TOGA_LIB_bin.FRAG"
              if ( $hmm_acc =~ /^TIGR/ );
            system "cat $dir_ind/$hmm_acc.FRAG       >> $dir_lib/PFAM_LIB.FRAG"
              if ( $hmm_acc =~ /^PF/ );

            # fix 3
            system
"$hmmconvert -A $dir_ind/${hmm_acc}_bin.FRAG $dir_lib/PFAM_LIB_bin.FRAG"

  #             "cat $dir_ind/${hmm_acc}_bin.FRAG >> $dir_lib/PFAM_LIB_bin.FRAG"
              if ( $hmm_acc =~ /^PF/ );
            unlink(
                "$dir_ind/${hmm_acc}_bin.FRAG.log",
                "$dir_ind/${hmm_acc}_bin.FRAG"
            );
        }
    }
    unlink(
        "$dir_lib/PFAM_LIB_bin.HMM.gsi", "$dir_lib/PFAM_LIB_bin.FRAG.gsi",
        "$dir_lib/TOGA_LIB_bin.HMM.gsi", "$dir_lib/TOGA_LIB_bin.FRAG.gsi",
        "$dir_lib/ALL_LIB_bin.HMM.gsi",  "$dir_lib/ALL_LIB_bin.FRAG.gsi"
    );

#     system "hmmindex $dir_lib/PFAM_LIB_bin.HMM  > $dir_lib/PFAM_LIB_bin.HMM.gsi.log"  if(-s "$dir_lib/PFAM_LIB_bin.HMM");
#     system "hmmindex $dir_lib/PFAM_LIB_bin.FRAG > $dir_lib/PFAM_LIB_bin.FRAG.gsi.log" if(-s "$dir_lib/PFAM_LIB_bin.FRAG");
#     system "hmmindex $dir_lib/TOGA_LIB_bin.HMM  > $dir_lib/TOGA_LIB_bin.HMM.gsi.log"  if(-s "$dir_lib/TOGA_LIB_bin.HMM");
#     system "hmmindex $dir_lib/TOGA_LIB_bin.FRAG > $dir_lib/TOGA_LIB_bin.FRAG.gsi.log" if(-s "$dir_lib/TOGA_LIB_bin.FRAG");
#     system "hmmindex $dir_lib/ALL_LIB_bin.HMM   > $dir_lib/ALL_LIB_bin.HMM.gsi.log"   if(-s "$dir_lib/ALL_LIB_bin.HMM");
#     system "hmmindex $dir_lib/ALL_LIB_bin.FRAG  > $dir_lib/ALL_LIB_bin.FRAG.gsi.log"  if(-s "$dir_lib/ALL_LIB_bin.FRAG");
}

sub build_lib {

# create XXX_LIB, XXX_LIB_bin, XXX_LIB_bin.YY.gsi files from individual HMM files
    use strict;
    my ( $family, $dir_ind, $dir_lib ) = @_;

# set hmmconvert command - as of LDhmmer-1.5.3, hmmconvert by default converts models to binary format
# Because of this they also deprecated the -b option, and in fact the program fails if you specify -b.
# But on Solaris machines, we're using regular hmmer, and thus regular hmmconvert so we need the -b.
# Thus:
    my $arch = `uname -s`;
    chomp $arch;
    my $hmmconvert;
    if ( $arch != /linux/i ) {
        die "Your operating system is impure: It does not support LDhmmer.\n";
    }
    $hmmconvert = "hmmconvert";
    if (   $family eq 'TIGR'
        || $family =~ /^ALL$/i )
    {
        unlink(
            "$dir_lib/TOGA_LIB.HMM",
            "$dir_lib/TOGA_LIB.FRAG",
            "$dir_lib/TOGA_LIB_bin.HMM",
            "$dir_lib/TOGA_LIB_bin.FRAG",
            "$dir_lib/TOGA_LIB_bin.HMM.log",
            "$dir_lib/TOGA_LIB_bin.FRAG.log",
            "$dir_lib/TOGA_LIB_bin.HMM.gsi.log",
            "$dir_lib/TOGA_LIB_bin.FRAG.gsi.log",
            "$dir_lib/TOGA_LIB_bin.HMM.gsi",
            "$dir_lib/TOGA_LIB_bin.FRAG.gsi"
        );
        while (<$dir_ind/TIGR*HMM>) {
            system
              "cat $_ >> $dir_lib/TOGA_LIB.HMM"; # if($_ =~ /^TIGR\d{5}\.HMM$/);
        }
        while (<$dir_ind/TIGR*FRAG>) {
            system "cat $_ >> $dir_lib/TOGA_LIB.FRAG"
              ;    # if($_ =~ /^TIGR\d{5}\.FRAG$/);
        }
        system
"$hmmconvert $dir_lib/TOGA_LIB.HMM  $dir_lib/TOGA_LIB_bin.HMM  > $dir_lib/TOGA_LIB_bin.HMM.log"
          if ( -s "$dir_lib/TOGA_LIB.HMM" );
        system
"$hmmconvert $dir_lib/TOGA_LIB.FRAG $dir_lib/TOGA_LIB_bin.FRAG > $dir_lib/TOGA_LIB_bin.FRAG.log"
          if ( -s "$dir_lib/TOGA_LIB.FRAG" );

# 	system "hmmindex $dir_lib/TOGA_LIB_bin.HMM  > $dir_lib/TOGA_LIB_bin.HMM.gsi.log"  if(-s "$dir_lib/TOGA_LIB_bin.HMM");
# 	system "hmmindex $dir_lib/TOGA_LIB_bin.FRAG > $dir_lib/TOGA_LIB_bin.FRAG.gsi.log" if(-s "$dir_lib/TOGA_LIB_bin.FRAG");
    }
    if (   $family eq 'PF'
        || $family =~ /^ALL$/i )
    {
        unlink(
            "$dir_lib/PFAM_LIB.HMM",
            "$dir_lib/PFAM_LIB.FRAG",
            "$dir_lib/PFAM_LIB_bin.HMM",
            "$dir_lib/PFAM_LIB_bin.FRAG",
            "$dir_lib/PFAM_LIB_bin.HMM.log",
            "$dir_lib/PFAM_LIB_bin.FRAG.log",
            "$dir_lib/PFAM_LIB_bin.HMM.gsi.log",
            "$dir_lib/PFAM_LIB_bin.FRAG.gsi.log",
            "$dir_lib/PFAM_LIB_bin.HMM.gsi",
            "$dir_lib/PFAM_LIB_bin.FRAG.gsi"
        );
        while (<$dir_ind/PF*HMM>) {
            system
              "cat $_ >> $dir_lib/PFAM_LIB.HMM";   # if($_ =~ /^PF\d{5}\.HMM$/);
        }
        while (<$dir_ind/PF*FRAG>) {
            system
              "cat $_ >> $dir_lib/PFAM_LIB.FRAG"; # if($_ =~ /^PF\d{5}\.FRAG$/);
        }
        system
"$hmmconvert $dir_lib/PFAM_LIB.HMM  $dir_lib/PFAM_LIB_bin.HMM  > $dir_lib/PFAM_LIB_bin.HMM.log"
          if ( -s "$dir_lib/PFAM_LIB.HMM" );
        system
"$hmmconvert $dir_lib/PFAM_LIB.FRAG $dir_lib/PFAM_LIB_bin.FRAG > $dir_lib/PFAM_LIB_bin.FRAG.log"
          if ( -s "$dir_lib/PFAM_LIB.FRAG" );

#	system "hmmindex $dir_lib/PFAM_LIB_bin.HMM  > $dir_lib/PFAM_LIB_bin.HMM.gsi.log"  if(-s "$dir_lib/PFAM_LIB_bin.HMM");
#	system "hmmindex $dir_lib/PFAM_LIB_bin.FRAG > $dir_lib/PFAM_LIB_bin.FRAG.gsi.log" if(-s "$dir_lib/PFAM_LIB_bin.FRAG");
    }
    if ( $family =~ /^ALL$/i ) {
        unlink(
            "$dir_lib/ALL_LIB.HMM",
            "$dir_lib/ALL_LIB.FRAG",
            "$dir_lib/ALL_LIB_bin.HMM",
            "$dir_lib/ALL_LIB_bin.FRAG",
            "$dir_lib/ALL_LIB_bin.HMM.log",
            "$dir_lib/ALL_LIB_bin.FRAG.log",
            "$dir_lib/ALL_LIB_bin.HMM.gsi.log",
            "$dir_lib/ALL_LIB_bin.FRAG.gsi.log",
            "$dir_lib/ALL_LIB_bin.HMM.gsi",
            "$dir_lib/ALL_LIB_bin.FRAG.gsi"
        );
        if (    -s "$dir_lib/TOGA_LIB.HMM"
            and -s "$dir_lib/PFAM_LIB.HMM" )
        {
            system
"cat $dir_lib/TOGA_LIB.HMM $dir_lib/PFAM_LIB.HMM >> $dir_lib/ALL_LIB.HMM";
        }
        elsif ( -s "$dir_lib/TOGA_LIB.HMM"
            and !-s "$dir_lib/PFAM_LIB.HMM" )
        {
            system "cp $dir_lib/TOGA_LIB.HMM $dir_lib/ALL_LIB.HMM";
        }
        elsif (!-s "$dir_lib/TOGA_LIB.HMM"
            and -s "$dir_lib/PFAM_LIB.HMM" )
        {
            system "cp $dir_lib/PFAM_LIB.HMM $dir_lib/ALL_LIB.HMM";
        }
        system
"$hmmconvert $dir_lib/ALL_LIB.HMM $dir_lib/ALL_LIB_bin.HMM > $dir_lib/ALL_LIB_bin.HMM.log"
          if ( -s "$dir_lib/ALL_LIB.HMM" );

#	system "hmmindex $dir_lib/ALL_LIB_bin.HMM > $dir_lib/ALL_LIB_bin.HMM.gsi.log" if(-s "$dir_lib/ALL_LIB_bin.HMM");
        if (    -s "$dir_lib/TOGA_LIB.FRAG"
            and -s "$dir_lib/PFAM_LIB.FRAG" )
        {
            system
"cat $dir_lib/TOGA_LIB.FRAG $dir_lib/PFAM_LIB.FRAG >> $dir_lib/ALL_LIB.FRAG";
        }
        elsif ( -s "$dir_lib/TOGA_LIB.FRAG" ) {
            system "cp $dir_lib/TOGA_LIB.FRAG $dir_lib/ALL_LIB.FRAG";
        }
        elsif ( -s "$dir_lib/PFAM_LIB.FRAG" ) {
            system "cp $dir_lib/PFAM_LIB.FRAG $dir_lib/ALL_LIB.FRAG";
        }
        system
"$hmmconvert $dir_lib/ALL_LIB.FRAG $dir_lib/ALL_LIB_bin.FRAG > $dir_lib/ALL_LIB_bin.FRAG.log"
          if ( -s "$dir_lib/ALL_LIB.FRAG" );

#	system "hmmindex $dir_lib/ALL_LIB_bin.FRAG > $dir_lib/ALL_LIB_bin.FRAG.gsi.log" if(-s "$dir_lib/ALL_LIB_bin.FRAG");
    }
}

sub get_process_id {

    #look for process ids defined by $process_name using system call "ps"
    use strict;
    my ( $proc_id_r, $proc_name, $log_header, $log_r ) = @_;
    open( FH, "ps -elf |" );
    while ( my $line = <FH> ) {
        chomp $line;

        #old matching not working
        if ( $line =~ /perl\s+(\S+\/)?\Q$proc_name\E/ ) {

            #	if($line =~ /\s+perl\s+(\S+\/)?\Q$proc_name\E(\s+|\Z)/) {
            my @tmp = split /\s+/, $line;
            $$proc_id_r{ $tmp[3] } = 1;
        }
    }
    close FH;
    for my $k ( keys %$proc_id_r ) {
        delete $$proc_id_r{$k}
          if ( $k !~ /^\d+$/ );
    }
    return "success";
}

sub load_hmm {
    use strict;

 # implement: need to get iso_id first if not assigned. this is the case for PF.
    my ( $db_proc, $db, $in_r, $log_header, $log_r ) = @_;
    my %status;
    &get_status_vocabulary( $db_proc, $db, \%status );

    # step 1: check data before loading
    for my $k ( keys %$in_r ) {
        next
          if ( ref $$in_r{$k} );
        $$in_r{$k} =~ s/^\s+|\s+$//g;
    }
    $log_header .= "--sub load_hmm";
    my $check = "success";

    # data that can not be null
    if ( $$in_r{iso_id} !~ /^\d+$/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--iso_id can not be null\n"
          if ($log_r);
    }
    if ( $$in_r{hmm_name} !~ /\S/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--hmm_name can not be null\n"
          if ($log_r);
    }
    if ( length( $$in_r{hmm_name} ) > 15 ) {
        $check = "error";
        $$log_r .=
          "***Error: $log_header--hmm_name must be shorter than 15 characters\n"
          if ($log_r);
    }
    if ( $$in_r{hmm_com_name} !~ /\S/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--hmm_com_name can not be null\n"
          if ($log_r);
    }
    if ( $$in_r{iso_type} !~ /^\w+$/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--iso_type can not be null\n"
          if ($log_r);
    }
    if ( $$in_r{hmm_type} !~ /^\w+$/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--hmm_type can not be null\n"
          if ($log_r);
    }
    if ( $$in_r{author} !~ /\S/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--author can not be null\n"
          if ($log_r);
    }
    if (   $$in_r{method_seed} !~ /\S/
        && $$in_r{hmm_id} !~ /^\d+$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--method_seed can not be null\n"
          if ($log_r);
    }
    if (   $$in_r{hmm} !~ /\S/
        && $$in_r{hmm_id} !~ /^\d+$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--hmm can not be null\n"
          if ($log_r);
    }
    if ( $$in_r{hmm_len} !~ /\S/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--hmm_len can not be null\n"
          if ($log_r);
    }
    if (   $$in_r{hmm_seed} !~ /\S/
        && $$in_r{seed_align_id} !~ /^\d+$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--hmm_seed can not be null\n"
          if ($log_r);
    }
    if ( $$in_r{is_current} !~ /\S/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--is_current can not be null\n"
          if ($log_r);
    }

    # check data format;
    if (   $$in_r{trusted_cutoff} =~ /\S/
        && $$in_r{trusted_cutoff} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for trusted_cutoff\n"
          if ($log_r);
    }
    if (   $$in_r{trusted_cutoff2} =~ /\S/
        && $$in_r{trusted_cutoff2} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for trusted_cutoff2\n"
          if ($log_r);
    }
    if (   $$in_r{noise_cutoff} =~ /\S/
        && $$in_r{noise_cutoff} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for noise_cutoff\n"
          if ($log_r);
    }
    if (   $$in_r{noise_cutoff2} =~ /\S/
        && $$in_r{noise_cutoff2} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for noise_cutoff2\n"
          if ($log_r);
    }
    if (   $$in_r{ec_num} =~ /\S/
        && $$in_r{ec_num} !~
        /^(\s*(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)\s*)*\s*$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for EC number\n"
          if ($log_r);
    }
    if ( $$in_r{is_current} !~ /^-?\d+$/ ) {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for is_current\n"
          if ($log_r);
    }
    if (   $$in_r{related_hmm} =~ /\S/
        && $$in_r{related_hmm} !~ /^(\s*(TIGR|PF)\d{5}\s*)*$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for related HMM\n"
          if ($log_r);
    }
    if (   $$in_r{hmm_len_check} =~ /^\d+$/
        && $$in_r{hmm_len} != $$in_r{hmm_len_check} )
    {
        $$log_r .=
"***Error: $log_header--hmm length from HMM($$in_r{hmm_len}) and BUILD($$in_r{hmm_len_check}) files are different\n"
          if ($log_r);
        $check = "error";
    }
    if (   $$in_r{num_seed} =~ /\S/
        && $$in_r{num_seed} !~ /^\d+$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for num_seed\n"
          if ($log_r);
    }
    if (   $$in_r{avg_score} =~ /\S/
        && $$in_r{avg_score} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for avg_score\n"
          if ($log_r);
    }
    if (   $$in_r{std_dev} =~ /\S/
        && $$in_r{std_dev} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for std_dev\n"
          if ($log_r);
    }
    if (   $$in_r{min_score} =~ /\S/
        && $$in_r{min_score} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for min_score\n"
          if ($log_r);
    }
    if (   $$in_r{max_score} =~ /\S/
        && $$in_r{max_score} !~ /^(-?\d+\.?\d*|\d*\.\d+)$/ )
    {
        $check = "error";
        $$log_r .= "***Error: $log_header--wrong format for max_score\n"
          if ($log_r);
    }

    #check whether iso_id is in isology
    my $ref =
      &do_query( $db_proc,
        "select iso_id from isology where iso_id = $$in_r{iso_id}",
        0, $log_header, $log_r );
    if ( @$ref != 1 ) {
        $$log_r .=
"***Error: $log_header--can not found the isology \#$$in_r{iso_id} in isology table\n"
          if ($log_r);
        return "error";
    }

    #check hmm id
    if ( $$in_r{hmm_id} =~ /^\d+$/ ) {
        my $ref =
          &do_query( $db_proc,
            "select hmm_acc from $db..hmm2 where id=$$in_r{hmm_id}",
            0, $log_header, $log_r );
        if ( @$ref == 0 ) {
            $log_r .=
"***Error: $log_header--can not find hmm id $$in_r{hmm_id}, assign new hmm id and hmm accession\n"
              if ($log_r);
            $$in_r{hmm_id}  = '';
            $$in_r{hmm_acc} = '';
        }
        elsif ($$in_r{hmm_acc} =~ /\S/
            && $$ref[0]{hmm_acc} ne $$in_r{hmm_acc} )
        {
            $log_r .=
"***Error: $log_header--hmm_acc $$in_r{hmm_acc} in INFO file is different from hmm_acc $$ref[0]{hmm_acc} found in DB as specified by hmm id $$in_r{hmm_id}, assign new hmm_id and hmm_acc\n"
              if ($log_r);
            $$in_r{hmm_id}  = '';
            $$in_r{hmm_acc} = '';
        }
    }

    #get type of HMM: global or fragment
    my $hmm_gf = '';
    if ( $$in_r{hmm} =~ /\S/ ) {
        &get_hmm_glo_frac( $$in_r{hmm}, \$hmm_gf, $$in_r{hmm_acc}, $log_header,
            $log_r );
    }
    elsif ( $$in_r{hmm_acc} ne '' ) {
        my $ref = &do_query(
            $db_proc,
"select hmm from $db..hmm2 where hmm_acc=\"$$in_r{hmm_acc}\" and is_current=1",
            0,
            $log_header,
            $log_r
        );
        &get_hmm_glo_frac( $$ref[0]{hmm}, \$hmm_gf, $$in_r{hmm_acc},
            $log_header, $log_r );
    }

    # insert or update hmm2
    if ( $$in_r{hmm_id} !~ /^\d+$/ ) {

        #check hmm_name to be unique
        #	goto XX;
        my $ref = &do_query( $db_proc,
            "select hmm_name, hmm_acc from $db..hmm2 where is_current=1" );
        for my $n (@$ref) {
            if ( $$n{hmm_name} eq $$in_r{hmm_name} ) {
                $$log_r .=
"***Error: $log_header--hmm name \"$$in_r{hmm_name}\" has been assigned to $$n{hmm_acc}, please rename your HMM\n"
                  if ($log_r);
                $check = 'error';
            }
        }
      XX:

        #check existing HMM acc for this iso_id:
        my $ref = &do_query( $db_proc,
            "select hmm_acc from $db..hmm2 where iso_id = $$in_r{iso_id}" );
        if ( @$ref > 0 ) {
            $$log_r .=
"***Error: $log_header--this iso_id has been assigned to one or more HMM acc before\n"
              if ($log_r);
            return "error";
        }
        my $ref = &do_query( $db_proc,
            "select iso_acc from $db..isology where iso_id = $$in_r{iso_id}" );
        if ( $$ref[0]{iso_acc} =~ /^(TIGR|PF)\d{5}$/ ) {
            $$log_r .=
"***Error: $log_header--this iso_id has been assigned to a HMM acc: $$ref[0]{iso_acc} before\n"
              if ($log_r);
            return "error";
        }
        return "error"
          if ( $check eq 'error' );

        #get new hmm_acc
        &convert_null(
            '',                     'null',
            \$$in_r{hmm_len},       \$$in_r{iso_id},
            \$$in_r{is_current},    \$$in_r{num_seed},
            \$$in_r{avg_score},     \$$in_r{std_dev},
            \$$in_r{min_score},     \$$in_r{max_score},
            \$$in_r{noise_cutoff},  \$$in_r{trusted_cutoff},
            \$$in_r{noise_cutoff2}, \$$in_r{trusted_cutoff2},
            \$$in_r{hmm_method_id}, \$$in_r{hmm_frag_method_id}
        );
        my $ref = &do_query( $db_proc,
"select max(hmm_acc) hmm_acc from $db..hmm2 where hmm_acc like \"TIGR%\""
        );
        $$in_r{hmm_acc} = (
            @$ref == 1 && $$ref[0]{hmm_acc} =~ /^TIGR(\d+)$/
            ? "TIGR" . sprintf( "%5.5d", $1 + 1 )
            : "TIGR00001"
        );

        # load hmm2 table
        my $hmm_comment   = $$in_r{hmm_comment};
        my $reference     = $$in_r{reference};
        my $hmm_com_name  = $$in_r{hmm_com_name};
        my $expanded_name = $$in_r{expanded_name};
        my $hmm           = $$in_r{hmm};
        my $hmm_frag      = $$in_r{hmm_frag};
        my $hmm_build     = $$in_r{hmm_build};
        my $private       = $$in_r{private};
        $hmm_comment   =~ s/\"/\"\"/g;
        $reference     =~ s/\"/\"\"/g;
        $hmm_com_name  =~ s/\"/\"\"/g;
        $expanded_name =~ s/\"/\"\"/g;
        $hmm           =~ s/\"/\"\"/g;
        $hmm_frag      =~ s/\"/\"\"/g;
        $hmm_build     =~ s/\"/\"\"/g;
        $private       =~ s/\"/\"\"/g;
        my $query_hmm2 =
"insert $db..hmm2 (hmm_acc, hmm_type, hmm_name, hmm_com_name, hmm, hmm_len, hmm_seed, hmm_frag, iso_id, search_prog, search_options, hmm_comment, is_current, author, prior, num_seed, avg_score, std_dev, min_score, max_score, ec_num, tc_num, iso_type, reference, noise_cutoff, trusted_cutoff, noise_cutoff2, trusted_cutoff2, hmm_build, private, gene_sym, expanded_name, hmm_method_id, hmm_frag_method_id, method_seed, entry_date, mod_date, hmm_mod_date, nomen_check) values (\"$$in_r{hmm_acc}\", \"$$in_r{hmm_type}\", \"$$in_r{hmm_name}\", \"$hmm_com_name\", \"$hmm\", $$in_r{hmm_len}, \"$$in_r{hmm_seed}\", \"$hmm_frag\", $$in_r{iso_id}, \"$$in_r{search_prog}\", \"$$in_r{search_options}\", \"$hmm_comment\", $$in_r{is_current}, \"$$in_r{author}\", \"$$in_r{prior}\", $$in_r{num_seed}, $$in_r{avg_score}, $$in_r{std_dev}, $$in_r{min_score}, $$in_r{max_score}, \"$$in_r{ec_num}\", \"$$in_r{tc_num}\", \"$$in_r{iso_type}\", \"$reference\", $$in_r{noise_cutoff}, $$in_r{trusted_cutoff}, $$in_r{noise_cutoff2}, $$in_r{trusted_cutoff2}, \"$hmm_build\", \"$private\", \"$$in_r{gene_sym}\", \"$expanded_name\", $$in_r{hmm_method_id}, $$in_r{hmm_frag_method_id}, \"$$in_r{method_seed}\", getdate(), getdate(), getdate(), $$in_r{nomen_check})";
        $query_hmm2 =~ s/, \"\", /, null, /g;
        &do_query( $db_proc, $query_hmm2, 1, $log_header, $log_r );
        &convert_null(
            'null',                 '',
            \$$in_r{hmm_len},       \$$in_r{iso_id},
            \$$in_r{is_current},    \$$in_r{num_seed},
            \$$in_r{avg_score},     \$$in_r{std_dev},
            \$$in_r{min_score},     \$$in_r{max_score},
            \$$in_r{noise_cutoff},  \$$in_r{trusted_cutoff},
            \$$in_r{noise_cutoff2}, \$$in_r{trusted_cutoff2},
            \$$in_r{hmm_method_id}, \$$in_r{hmm_frag_method_id}
        );
        my $ref = &do_query( $db_proc,
"select id, hmm_name, hmm_com_name from $db..hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}"
        );

        if (   @$ref == 0
            or $$ref[0]{hmm_name} ne $$in_r{hmm_name}
            or $$ref[0]{hmm_com_name} ne $$in_r{hmm_com_name} )
        {
            $$log_r .=
"***Error: $log_header--failed to load hmm $$in_r{hmm_acc} to HMM2 table\n"
              if ($log_r);
            return "error";
        }
        else {
            $$log_r .=
              "Message: $log_header--new HMM is inserted into hmm2 table\n"
              if ($log_r);
            $$in_r{hmm_id} = $$ref[0]{id};
        }
    }
    else {
        return "error"
          if ( $check eq 'error' );
        &convert_null(
            '',                     'null',
            \$$in_r{hmm_len},       \$$in_r{iso_id},
            \$$in_r{is_current},    \$$in_r{num_seed},
            \$$in_r{avg_score},     \$$in_r{std_dev},
            \$$in_r{min_score},     \$$in_r{max_score},
            \$$in_r{noise_cutoff},  \$$in_r{trusted_cutoff},
            \$$in_r{noise_cutoff2}, \$$in_r{trusted_cutoff2},
            \$$in_r{hmm_method_id}, \$$in_r{hmm_frag_method_id}
        );
        my $hmm_comment   = $$in_r{hmm_comment};
        my $reference     = $$in_r{reference};
        my $hmm_com_name  = $$in_r{hmm_com_name};
        my $expanded_name = $$in_r{expanded_name};
        my $hmm           = $$in_r{hmm};
        my $hmm_frag      = $$in_r{hmm_frag};
        my $hmm_build     = $$in_r{hmm_build};
        my $private       = $$in_r{private};
        $hmm_comment   =~ s/\"/\"\"/g;
        $reference     =~ s/\"/\"\"/g;
        $hmm_com_name  =~ s/\"/\"\"/g;
        $expanded_name =~ s/\"/\"\"/g;
        $hmm           =~ s/\"/\"\"/g;
        $hmm_frag      =~ s/\"/\"\"/g;
        $hmm_build     =~ s/\"/\"\"/g;
        $private       =~ s/\"/\"\"/g;

        #get old info from db
        my %data;
        &get_INFO_from_db( $db_proc, $db, $$in_r{hmm_acc}, "hmm_acc", \%data );
        if ( $data{hmm_id} != $$in_r{hmm_id} ) {
            print
"***Error: mismatch between HMM acc ($$in_r{hmm_acc}) and hmm_id ($data{hmm_id}, $$in_r{hmm_id})\n";
            return;
        }
        my $query_hmm2 = "update $db..hmm2 set ";
        $query_hmm2 .= "hmm_type=\"$$in_r{hmm_type}\", "
          if ( $$in_r{hmm_type} ne $data{hmm_type} );
        $query_hmm2 .= "hmm_name=\"$$in_r{hmm_name}\", "
          if ( $$in_r{hmm_name} ne $data{hmm_name} );
        $query_hmm2 .= "hmm_com_name=\"$hmm_com_name\", "
          if ( $hmm_com_name ne $data{hmm_com_name} );
        $query_hmm2 .= "hmm_len=$$in_r{hmm_len}, "
          if ( $$in_r{hmm_len} != $data{hmm_len} );
        $query_hmm2 .= "iso_id=$$in_r{iso_id}, "
          if ( $$in_r{iso_id} != $data{iso_id} );
        $query_hmm2 .= "search_prog=\"$$in_r{search_prog}\", "
          if ( $$in_r{search_prog} ne $data{search_prog} );
        $query_hmm2 .= "search_options=\"$$in_r{search_options}\", "
          if ( $$in_r{search_options} ne $data{search_options} );
        $query_hmm2 .= "hmm_comment=\"$hmm_comment\", "
          if ( $hmm_comment ne $data{hmm_comment} );
        $query_hmm2 .= "is_current=$$in_r{is_current}, "
          if ( $$in_r{is_current} != $data{is_current} );
        $query_hmm2 .= "author=\"$$in_r{author}\", "
          if ( $$in_r{author} ne $data{author} );
        $query_hmm2 .= "prior=\"$$in_r{prior}\", "
          if ( $$in_r{prior} ne $data{prior} );
        $query_hmm2 .= "num_seed=$$in_r{num_seed}, "
          if ( $$in_r{num_seed} != $data{num_seed} );
        $query_hmm2 .= "avg_score=$$in_r{avg_score}, "
          if ( $$in_r{avg_score} != $data{avg_score} );
        $query_hmm2 .= "std_dev=$$in_r{std_dev}, "
          if ( $$in_r{std_dev} != $data{std_dev} );
        $query_hmm2 .= "min_score=$$in_r{min_score}, "
          if ( $$in_r{min_score} != $data{min_score} );
        $query_hmm2 .= "max_score=$$in_r{max_score}, "
          if ( $$in_r{max_score} != $data{max_score} );
        $query_hmm2 .= "ec_num=\"$$in_r{ec_num}\", "
          if ( $$in_r{ec_num} ne $data{ec_num} );
        $query_hmm2 .= "tc_num=\"$$in_r{tc_num}\", "
          if ( $$in_r{tc_num} ne $data{tc_num} );

#	$query_hmm2 .= "struct_acc=\"$$in_r{struct_acc}\", " if($$in_r{struct_acc} ne $data{struct_acc});
        $query_hmm2 .= "iso_type=\"$$in_r{iso_type}\", "
          if ( $$in_r{iso_type} ne $data{iso_type} );
        $query_hmm2 .= "reference=\"$reference\", "
          if ( $reference ne $data{reference} );
        $query_hmm2 .= "noise_cutoff=$$in_r{noise_cutoff}, "
          if ( $$in_r{noise_cutoff} != $data{noise_cutoff} );
        $query_hmm2 .= "trusted_cutoff=$$in_r{trusted_cutoff}, "
          if ( $$in_r{trusted_cutoff} != $data{trusted_cutoff} );
        $query_hmm2 .= "noise_cutoff2=$$in_r{noise_cutoff2}, "
          if ( $$in_r{noise_cutoff2} != $data{noise_cutoff2} );
        $query_hmm2 .= "trusted_cutoff2=$$in_r{trusted_cutoff2}, "
          if ( $$in_r{trusted_cutoff2} != $data{trusted_cutoff2} );
        $query_hmm2 .= "private=\"$private\", "
          if ( $private ne $data{private} );
        $query_hmm2 .= "gene_sym=\"$$in_r{gene_sym}\", "
          if ( $$in_r{gene_sym} ne $data{gene_sym} );
        $query_hmm2 .= "expanded_name=\"$expanded_name\", "
          if ( $expanded_name ne $data{expanded_name} );
        $query_hmm2 .= "hmm_method_id=$$in_r{hmm_method_id}, "
          if ( $$in_r{hmm_method_id} != $data{hmm_method_id} );
        $query_hmm2 .= "hmm_frag_method_id=$$in_r{hmm_frag_method_id}, "
          if ( $$in_r{hmm_frag_method_id} != $data{hmm_frag_method_id} );
        $query_hmm2 .= "method_seed=\"$$in_r{method_seed}\", "
          if ( $$in_r{method_seed} ne $data{method_seed} );
        $query_hmm2 .= "nomen_check=$$in_r{nomen_check}, "
          if ( $$in_r{nomen_check} != $data{nomen_check} );
        $query_hmm2 .= "hmm=\"$hmm\", "
          if (  $$in_r{hmm} =~ /\S/
            and $$in_r{hmm} ne $data{hmm} );
        $query_hmm2 .= "hmm_frag=\"$hmm_frag\", "
          if (  $$in_r{hmm_frag} =~ /\S/
            and $$in_r{hmm_frag} ne $data{hmm_frag} );
        $query_hmm2 .= "hmm_build=\"$hmm_build\", "
          if (  $$in_r{hmm_build} =~ /\S/
            and $$in_r{hmm_build} ne $data{hmm_build} );
        $query_hmm2 .= "hmm_seed=\"$$in_r{hmm_seed}\", "
          if (  $$in_r{hmm_seed} =~ /\S/
            and $$in_r{hmm_seed} ne $data{hmm_seed} );
        $query_hmm2 =~ s/, $//;
        $query_hmm2 .= " where id=$$in_r{hmm_id}";
        $query_hmm2 =~ s/=\"\", /=null, /g;
        &do_query( $db_proc, $query_hmm2, 1, $log_header, $log_r )
          if ( $query_hmm2 !~ /set\s+where\s+/ );
        &convert_null(
            'null',                 '',
            \$$in_r{hmm_len},       \$$in_r{iso_id},
            \$$in_r{is_current},    \$$in_r{num_seed},
            \$$in_r{avg_score},     \$$in_r{std_dev},
            \$$in_r{min_score},     \$$in_r{max_score},
            \$$in_r{noise_cutoff},  \$$in_r{trusted_cutoff},
            \$$in_r{noise_cutoff2}, \$$in_r{trusted_cutoff2},
            \$$in_r{hmm_method_id}, \$$in_r{hmm_frag_method_id}
        );
    }

    # assign HMM model property: global, fragment, single, multiple
    if ( $hmm_gf ne '' ) {
        my $ref = &do_query( $db_proc,
"select id, status_vocab_id from $db..status where link_id=$$in_r{hmm_id} and status_vocab_id in ($status{hmm2}{hmm}{multi_frag}, $status{hmm2}{hmm}{multi_global}, $status{hmm2}{hmm}{single_frag}, $status{hmm2}{hmm}{single_global})"
        );
        if ( @$ref == 0 ) {
            &do_query(
                $db_proc,
"insert status (link_id, status_vocab_id, person) values ($$in_r{hmm_id}, $status{hmm2}{hmm}{$hmm_gf}, \"$$in_r{assigned_by}\")",
                1,
                $log_header,
                $log_r
            );
        }
        elsif ( @$ref == 1 ) {
            &do_query(
                $db_proc,
"update status set status_vocab_id=$status{hmm2}{hmm}{$hmm_gf} where link_id=$$in_r{hmm_id} and status_vocab_id in ($status{hmm2}{hmm}{multi_frag}, $status{hmm2}{hmm}{multi_global}, $status{hmm2}{hmm}{single_frag}, $status{hmm2}{hmm}{single_global})",
                1,
                $log_header,
                $log_r
            );
        }
        else {
            $$log_r .=
"***Error: $log_header--more than one status entry for the global/fragment/single/multiple property of HMM $$in_r{hmm_acc}, delete them all and insert new one\n"
              if ($log_r);
            &do_query(
                $db_proc,
"delete status where link_id=$$in_r{hmm_id} and status_vocab_id in ($status{hmm2}{hmm}{multi_frag}, $status{hmm2}{hmm}{multi_global}, $status{hmm2}{hmm}{single_frag}, $status{hmm2}{hmm}{single_global})",
                1,
                $log_header,
                $log_r
            );
            &do_query(
                $db_proc,
"insert status (link_id, status_vocab_id, person) values ($$in_r{hmm_id}, $status{hmm2}{hmm}{$hmm_gf}, \"$$in_r{assigned_by}\")",
                1,
                $log_header,
                $log_r
            );
        }
    }

    #update isology table
    # no ref_id is loaded
    $$in_r{iso_acc} = $$in_r{hmm_acc}
      if ( $$in_r{iso_acc} !~ /\S/ );
    $$in_r{user_class} = 0
      if ( $$in_r{user_class} !~ /^\d+$/ );
    my $class_name  = $$in_r{class_name};
    my $iso_comment = $$in_r{iso_comment};
    $class_name  =~ s/\"/\"\"/g;
    $iso_comment =~ s/\"/\"\"/g;
    my $query_isology =
"update isology set class_name=\"$class_name\", seq_type=\"$$in_r{seq_type}\", iso_type=\"$$in_r{iso_iso_type}\", method=\"$$in_r{iso_method}\", comment=\"$iso_comment\", assigned_by=\"$$in_r{assigned_by}\", user_class=$$in_r{user_class}, iso_acc=\"$$in_r{iso_acc}\", db_type=\"$$in_r{db_type}\" where iso_id=$$in_r{iso_id}";
    $query_isology =~ s/=\"\", /=null, /g;
    &do_query( $db_proc, $query_isology, 1, $log_header, $log_r );
    my $ref = &do_query( $db_proc,
        "select iso_acc from $db..isology where iso_id = $$in_r{iso_id}" );

    if (   @$ref != 1
        or $$ref[0]{iso_acc} ne $$in_r{iso_acc} )
    {
        $$log_r .=
"***Error: $log_header--failed to load iso_id $$in_r{iso_id} to isology table\n"
          if ($log_r);
        &do_query(
            $db_proc,
"delete hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}",
            1
        );
        $$log_r .=
          "Message: $log_header--HMM $$in_r{hmm_acc} is deleted from hmm2 table"
          if ($log_r);
        return "error";
    }

    #load SEED to alignment table
    my $ref = &do_query(
        $db_proc,
"select align_id from alignment where iso_id=$$in_r{iso_id} and align_type=\"seed\"",
        0,
        $log_header,
        $log_r
    );
    &convert_null( '', 'null', \$$in_r{seed_ref_id}, \$$in_r{seed_user_class} );
    if ( @$ref == 0 ) {
        my $seed_align_name = $$in_r{seed_align_name};
        my $seed_comment    = $$in_r{seed_comment};
        $seed_align_name =~ s/\"/\"\"/g;
        $seed_comment    =~ s/\"/\"\"/g;
        my $query_alignment =
"insert $db..alignment (align_name, alignment, method, comment, ref_id, assigned_by, user_class, iso_id, align_type) values (\"$seed_align_name\", \"$$in_r{seed_align}\", \"$$in_r{method_seed}\", \"$seed_comment\", $$in_r{seed_ref_id}, \"$$in_r{seed_align_assigned_by}\", $$in_r{seed_user_class}, $$in_r{iso_id}, \"$$in_r{seed_align_type}\")";
        my $ref =
          &do_query( $db_proc, $query_alignment, 0, "$log_header--seed--",
            $log_r );
        if ( @$ref != 1 ) {
            $$log_r .=
"***Error: $log_header--failed to load SEED alignment of iso_id $$in_r{iso_id} to alignment table\n"
              if ($log_r);
            &do_query(
                $db_proc,
"delete hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}",
                1,
                $log_header,
                $log_r
            );
            $$log_r .=
"Message: $log_header--HMM $$in_r{hmm_acc} is deleted from hmm2 table"
              if ($log_r);
            return "error";
        }
        else {
            $$in_r{seed_align_id} = $$ref[0]{'COL(1)'};
        }
    }
    elsif ( @$ref == 1 ) {
        if (   $$in_r{seed_align_id} =~ /^\d+$/
            && $$ref[0]{align_id} != $$in_r{seed_align_id} )
        {
            $$log_r .=
"***Error: $log_header--the seed align_id $$ref[0]{align_id} for iso_id $$in_r{iso_id} in DB is different from the seed align_id $$in_r{seed_align_id} in INFO file\n"
              if ($log_r);
            &do_query(
                $db_proc,
"delete hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}",
                1,
                $log_header,
                $log_r
            );
            $$log_r .=
"Message: $log_header--HMM $$in_r{hmm_acc} is deleted from hmm2 table"
              if ($log_r);
            return "error";
        }
        my $seed_align_name = $$in_r{seed_align_name};
        my $seed_comment    = $$in_r{seed_comment};
        $seed_align_name =~ s/\"/\"\"/g;
        $seed_comment    =~ s/\"/\"\"/g;
        my $query_alignment =
"update alignment set align_name=\"$seed_align_name\", method=\"$$in_r{method_seed}\", comment=\"$seed_comment\", ref_id=$$in_r{seed_ref_id}, assigned_by=\"$$in_r{seed_align_assigned_by}\", user_class=$$in_r{seed_user_class}, iso_id=$$in_r{iso_id}, align_type=\"$$in_r{seed_align_type}\"";
        $query_alignment .= ", alignment=\"$$in_r{seed_align}\""
          if ( $$in_r{seed_align} =~ /\S/ );
        $query_alignment .= " where align_id = $$ref[0]{align_id}";
        &do_query( $db_proc, $query_alignment, 1, "$log_header--seed--",
            $log_r );
    }
    else {
        $$log_r .=
"***Error: $log_header--this iso_id $$in_r{iso_id} has too many seed alignments linked to it\n"
          if ($log_r);
        &do_query(
            $db_proc,
"delete hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}",
            1
        );
        $$log_r .=
          "Message: $log_header--HMM $$in_r{hmm_acc} is deleted from hmm2 table"
          if ($log_r);
        return "error";
    }
    &convert_null( 'null', '', \$$in_r{seed_ref_id}, \$$in_r{seed_user_class} );

    #load FULL to alignment table
    my $ref = &do_query(
        $db_proc,
"select align_id from alignment where iso_id=$$in_r{iso_id} and align_type=\"all\"",
        0,
        $log_header,
        $log_r
    );
    &convert_null( '', 'null', \$$in_r{full_ref_id} );
    if ( @$ref == 0 ) {
        my $full_align_name = $$in_r{full_align_name};
        my $full_comment    = $$in_r{full_comment};
        $full_align_name =~ s/\"/\"\"/g;
        $full_comment    =~ s/\"/\"\"/g;
        my $query_alignment =
"insert $db..alignment (align_name, alignment, method, comment, ref_id, assigned_by, user_class, iso_id, align_type) values (\"$full_align_name\", \"$$in_r{full_align}\", \"$$in_r{method_full}\", \"$full_comment\", $$in_r{full_ref_id}, \"$$in_r{full_align_assigned_by}\", $$in_r{full_user_class}, $$in_r{iso_id}, \"$$in_r{full_align_type}\")";
        my $ref =
          &do_query( $db_proc, $query_alignment, 0, "$log_header--full--",
            $log_r );
        if ( @$ref != 1 ) {
            $$log_r .=
"***Error: $log_header--failed to load FULL alignment of iso_id $$in_r{iso_id} to alignment table\n"
              if ($log_r);
            &do_query(
                $db_proc,
"delete hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}",
                1,
                $log_header,
                $log_r
            );
            $$log_r .=
"Message: $log_header--HMM $$in_r{hmm_acc} is deleted from hmm2 table"
              if ($log_r);
            return "error";
        }
        else {
            $$in_r{full_align_id} = $$ref[0]{'COL(1)'};
        }
    }
    elsif ( @$ref == 1 ) {
        if (   $$in_r{full_align_id} =~ /^\d+$/
            && $$ref[0]{align_id} != $$in_r{full_align_id} )
        {
            $$log_r .=
"***Error: $log_header--the full align_id $$ref[0]{align_id} for iso_id $$in_r{iso_id} in DB is different from the full align_id $$in_r{full_align_id} in INFO file\n"
              if ($log_r);
            &do_query(
                $db_proc,
"delete hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}",
                1,
                $log_header,
                $log_r
            );
            $$log_r .=
"Message: $log_header--HMM $$in_r{hmm_acc} is deleted from hmm2 table"
              if ($log_r);
            return "error";
        }
        my $full_align_name = $$in_r{full_align_name};
        my $full_comment    = $$in_r{full_comment};
        $full_align_name =~ s/\"/\"\"/g;
        $full_comment    =~ s/\"/\"\"/g;
        my $query_alignment =
"update alignment set align_name=\"$full_align_name\", method=\"$$in_r{method_full}\", comment=\"$full_comment\", ref_id=$$in_r{full_ref_id}, assigned_by=\"$$in_r{full_align_assigned_by}\", user_class=$$in_r{full_user_class}, iso_id=$$in_r{iso_id}, align_type=\"$$in_r{full_align_type}\"";
        $query_alignment .= ", alignment=\"$$in_r{full_align}\""
          if ( $$in_r{full_align} =~ /\S/ );
        $query_alignment .= " where align_id = $$ref[0]{align_id}";
        &do_query( $db_proc, $query_alignment, 1, "$log_header--full--",
            $log_r );
    }
    else {
        $$log_r .=
"***Error: $log_header--this iso_id $$in_r{iso_id} has too many full alignments linked to it\n"
          if ($log_r);
        &do_query(
            $db_proc,
"delete hmm2 where hmm_acc = \"$$in_r{hmm_acc}\" and is_current = $$in_r{is_current}",
            1
        );
        $$log_r .=
          "Message: $log_header--HMM $$in_r{hmm_acc} is deleted from hmm2 table"
          if ($log_r);
        return "error";
    }
    &convert_null( 'null', '', \$$in_r{full_ref_id} );

    #load iso_link
    if ( $$in_r{full}{number} > 0 ) {
        my @iso_link;
        &construct_iso_link( $$in_r{full}, '', \@iso_link );
        &load_iso_link( 1, $$in_r{iso_id}, \@iso_link, $db, $db_proc, '', '',
            '' );
    }
}

sub get_hmm_glo_frac {
    use strict;
    my ( $hmm, $hmm_gf_r, $hmm_acc, $log_header, $log_r ) = @_;
    $log_header .= "--sub get_hmm_glo_frac--hmm_acc=$hmm_acc";
    for my $line ( split /\n/, $hmm ) {
        if ( $line =~ /hmmbuild\s+(.+?\s+-f|-f)\s+/ ) {
            $$hmm_gf_r = "multi_frag";
            last;
        }
        elsif ( $line =~ /hmmbuild\s+(.+?\s+-g|-g)\s+/ ) {
            $$hmm_gf_r = "single_global";
            last;
        }
        elsif ( $line =~ /hmmbuild\s+(.+?\s+-s|-s)\s+/ ) {
            $$hmm_gf_r = "single_frag";
            last;
        }
    }
    $$hmm_gf_r = "multi_global"
      if ( $$hmm_gf_r eq '' );
    return;
}

sub convert_null {
    use strict;
    my ( $from, $to, @arr ) = @_;
    for my $n (@arr) {
        $$n = $to
          if ( $$n eq $from );
    }
}

sub get_status_vocabulary {
    use strict;
    my ( $db_proc, $db, $status_r ) = @_;
    my $ref = &do_query(
        $db_proc,
"select id, status_type, table_name, field_name, color from $db..status_vocabulary where
            (table_name=\"iso_project\"      and field_name=\"project_id\" ) or
            (table_name=\"isology\"          and field_name=\"iso_id\"     ) or
            (table_name=\"alignment\"        and field_name=\"align_id\"   ) or
            (table_name=\"project_iso_link\" and field_name=\"id\"         ) or
            (table_name=\"characterized\"    and field_name=\"id\"         ) or
            (table_name=\"groups\"           and field_name=\"group_type\" ) or
            (table_name=\"hmm2\"             and field_name=\"iso_type\"   ) or
            (table_name=\"hmm2\"             and field_name=\"hmm\"        ) or
            (table_name=\"hmm2\"             and field_name=\"hmm_acc\"    ) or
            (table_name=\"hmm2\"             and field_name=\"hmm_type\"   ) or
            (table_name=\"hmm2\"             and field_name=\"method_seed\") or
            (table_name=\"hmm2\"             and field_name=\"hmm_mod_date\") or
            (table_name=\"hmm2\"             and field_name=\"related_hmm\") or
            (table_name=\"tree\"             and field_name=\"tree_id\")"
    );
    for (@$ref) {
        $$status_r{ $$_{table_name} }{ $$_{field_name} }{ $$_{status_type} } =
          $$_{id};
    }
}

sub create_mini_db_info {

    #create mini_db for a hmm in directory $dir/
    use strict;
    my ( $info_r, $hmm_acc, $dir, $format, $status_r, $nraa_size ) = @_;
    $$info_r{db_proc}->disconnect
      if ( defined $$info_r{db_proc} );
    $$info_r{db_proc} = &connect_db(
        $$info_r{db},   $$info_r{db_type}, $$info_r{server},
        $$info_r{user}, $$info_r{password}
    );
    my $ref = &doSQL(
        $info_r,
"select hmm_seed from $$info_r{db}..hmm2 where hmm_acc = \"$hmm_acc\" and is_current = 1",
        0,
        60
    );
    my %seq;
    $seq{ori} = $$ref[0]{hmm_seed};
    &parse_sequence( \%seq, 'mul', '', '', '', '', '' );
    my $tmp = &tmp_file( $dir, "tmp" );
    open( FH, ">$dir/${tmp}_$hmm_acc.seed" );
    my $out = &format_sequence( \%seq, 'fasta', '', '', 1, '', '', '' );
    print FH $out;
    close FH;

#    my $blast_result = system("/usr/local/projects/OGC/egad_test/working/blast_wrapper.pl -S $dir/${tmp}_$hmm_acc.seed -e new -d $dir -t ${tmp}_mini_db_$hmm_acc");
    my $blast_result = system(
"$ENV{HMM_SCRIPTS}/blast_wrapper.pl -S $dir/${tmp}_$hmm_acc.seed -e new -d $dir -t ${tmp}_mini_db_$hmm_acc"
    );
    system "/bin/mv $dir/${tmp}_mini_db_$hmm_acc $dir/mini_db_$hmm_acc";
    unlink(
        "$dir/${tmp}_$hmm_acc.seed",
        "$dir/${tmp}_$hmm_acc.seed.blast_wrapper_log",
        "$dir/$tmp"
    );

    if ( $blast_result != 25856 ) {    # means blast_wrapper failed
        unlink("$dir/mini_db_$hmm_acc");
        $$info_r{db_proc}->disconnect
          if ( defined $$info_r{db_proc} );
        return;
    }
    my %mini_db;
    &read_sequence( \%mini_db, "$dir/mini_db_$hmm_acc", '', '', '' );
    my $out;
    &merge_mini_db( '', \%mini_db, $format, \$out, '', '' );
    open( FH, ">$dir/mini_db_$hmm_acc" );
    print FH $out;
    close FH;
    &doSQL(
        $info_r,
"update $$info_r{db}..status set status_int = $nraa_size, link_id=h.id, status_char=h.hmm_mod_date from $$info_r{db}..hmm2 h, $$info_r{db}..status s where s.method=\"$hmm_acc\" and s.status_vocab_id=$$status_r{hmm2}{hmm_acc}{MINI_DB} and s.method=h.hmm_acc and h.is_current=1",
        1,
        60
    );
    $$info_r{db_proc}->disconnect
      if ( defined $$info_r{db_proc} );

#    print LOG "${\(scalar localtime)}: the child process building mini_db of $hmm_acc is done and should exit now\n";
    return;
}

sub merge_mini_db {

#merge mini_db1 to mini_db2. For change format of a mini_db file, set it as $db2 and set $db1 to null
# $db1_r and $db2_r has the same data structure as STOCKHOLM mul format for multiple sequences
# $format is "long" or "short"
# Note: since one prot may have different accession, it is important to make sure all prot using the
# same priority for assigning accession, i.e. YANK's EGAD->SP->GP->PIR etc.
    use strict;
    my ( $db1_r, $db2_r, $format, $out_r, $log_header, $log_r ) = @_;
    $log_header .= "--sub merge_mini_db";
    if ( $db1_r ne '' ) {
        my $db2_num = $$db2_r{stockholm}{GF}{mini_db}{number};
        for my $n ( 0 .. $$db1_r{stockholm}{GF}{mini_db}{number} - 1 ) {
            if ( $$db1_r{stockholm}{GF}{mini_db}{$n}{acc} =~
                /^([^\|]+\|[^\|\/]+)/ )
            {
                my $db1_acc = $1;
                my $found   = 0;
                for my $m ( 0 .. $db2_num - 1 ) {
                    if ( $$db2_r{stockholm}{GF}{mini_db}{$m}{acc} =~
                        /^\Q$db1_acc\E/ )
                    {
                        $found = 1;
                        last;
                    }
                }
                if ( $found == 0 ) {
                    $$db2_r{stockholm}{GF}{mini_db}{number} = 0
                      if ( !defined $$db2_r{stockholm}{GF}{mini_db}{number}
                        or $$db2_r{stockholm}{GF}{mini_db}{number} eq '' );
                    %{ $$db2_r{stockholm}{GF}{mini_db}
                          { $$db2_r{stockholm}{GF}{mini_db}{number} } } =
                      %{ $$db1_r{stockholm}{GF}{mini_db}{$n} };
                    ++$$db2_r{stockholm}{GF}{mini_db}{number};
                }
            }
        }
        my $db2_num = $$db2_r{number};
        for my $n ( 0 .. $$db1_r{number} - 1 ) {
            if ( $$db1_r{$n}{header} =~ /^([^\|]+\|[^\|\/]+)/ ) {
                my $db1_acc = $1;
                my $found   = 0;
                for my $m ( 0 .. $db2_num - 1 ) {
                    if ( $$db2_r{$m}{header} =~ /^\Q$db1_acc\E/ ) {
                        $found = 1;
                        last;
                    }
                }
                if ( $found == 0 ) {
                    %{ $$db2_r{ $$db2_r{number} } } = %{ $$db1_r{$n} };
                    ++$$db2_r{number};
                }
            }
        }
    }

    #    print "$$db1_r{number}, $$db2_r{number}\n";
    if ( $format eq 'short' ) {
        for my $n ( 0 .. $$db2_r{number} - 1 ) {
            $$db2_r{$n}{seq}     = '';
            $$db2_r{$n}{seq_gap} = '';
        }
        $$out_r = &format_sequence( $db2_r, 'mul', '', '', '', '', '', '' );
        $$out_r .= "\n";
        for my $n ( 0 .. $$db2_r{number} - 1 ) {
            if ( $$db2_r{$n}{header} =~ /^([^\|]+\|[^\|\/]+)/ ) {
                my $new_acc = $1;
                $$out_r .= "$1\n";
            }
        }
    }
    elsif ( $format eq 'long' ) {
        my %tmp_hash;
        for my $n ( 0 .. $$db2_r{number} - 1 ) {
            if ( $$db2_r{$n}{header} =~ /^([^\|]+\|[^\|\/]+)/ ) {
                %{ $tmp_hash{$1} } = ();
            }
            else {
                print "***Error1<br>\n";
            }
        }

        #problem here: sequences not in database are lost from merged mini_db
        &download_nraa( \%tmp_hash, '', '', 1, '', '' );

        #	&yank_prot(\%tmp_hash);
        for my $n ( 0 .. $$db2_r{number} - 1 ) {
            if ( $$db2_r{$n}{header} =~ /^([^\|]+\|[^\|\/]+)/ ) {
                my $acc = $1;
                $$db2_r{$n}{seq} = $tmp_hash{$acc}{seq}
                  if ( $tmp_hash{$acc}{seq} =~ /\S/ );
                $$db2_r{$n}{seq_gap} = '';
            }
            else {
                print "***Error2<br>\n";
            }
        }
        $$out_r = &format_sequence( $db2_r, 'mul', '', '', '', '', '', '' );
    }
    return 0;
}

sub shrink_seq {

    #re-numbering sequences if some sequences are deleted
    use strict;
    my ($seq_r) = @_;
    my $index = 0;
    for my $n ( 0 .. $$seq_r{number} - 1 ) {
        if ( defined $$seq_r{$n} ) {
            $$seq_r{$index} = $$seq_r{$n};

            #	    %{$$seq_r{$index}} = %{$$seq_r{$n}};
            ++$index;
        }
    }
    $$seq_r{number} = $index;
}

sub shrink_alignment {

    #shrink the length of sequences if all sequences have a gap at a position
    use strict;
    my ($seq_r) = @_;
    my $len     = length( $$seq_r{0}{seq_gap} );
    my $m       = 0;
    while ( $m < $len ) {
        my $is_gap = 1;
        for my $n ( 0 .. $$seq_r{number} - 1 ) {
            if ( substr( $$seq_r{$n}{seq_gap}, $m, 1 ) !~ /[.\-]/ ) {
                $is_gap = 0;
                last;
            }
        }
        if ( $is_gap == 1 ) {
            for my $n ( 0 .. $$seq_r{number} - 1 ) {
                $$seq_r{$n}{seq_gap} =
                    substr( $$seq_r{$n}{seq_gap}, 0, $m )
                  . substr( $$seq_r{$n}{seq_gap}, $m + 1 );
            }
            --$len;
        }
        else {
            ++$m;
        }
    }
}

sub add_prot_frag {

# $data_r holds coordinates of non-overlapping fragments of a protein
# $data_r is a array reference. Each array element is a hash reference with keys: "lend" and "rend"
# a new fragment ($lend, $rend) is added to $data_r.  Coordinates of fragments are re-calculated
# after adding the new fragment
    use strict;
    my ( $data_r, $lend, $rend, $add ) = @_;
    my ( $overlapped, @new_data ) = ( 0, () );
    $new_data[0] = {
        "lend" => $lend,
        "rend" => $rend
    };
    for my $n (@$data_r) {
        my $found = 0;
        for my $m (@new_data) {
            my ( $new_lend, $new_rend ) =
              &overlap_fragment( $$n{lend}, $$n{rend}, $$m{lend}, $$m{rend} );
            if ( $new_lend ne '' ) {
                ( $$m{lend}, $$m{rend} ) = ( $new_lend, $new_rend );
                $found      = 1;
                $overlapped = 1;
                last;
            }
        }
        push @new_data,
          {
            "lend" => $$n{lend},
            "rend" => $$n{rend}
          }
          if ( $found == 0 );
    }
    @$data_r = @new_data
      if ($add);
    return $overlapped;
}

sub overlap_fragment {

    # return null if the two fragments are not overlapping.
    # return expanded coordinates if the two fragments are overlapping
    use strict;
    my ( $lend1, $rend1, $lend2, $rend2 ) = @_;
    my ( $new_lend, $new_rend ) = ( '', '' );
    if (   $lend1 <= $lend2
        && $rend1 >= $rend2 )
    {
        ( $new_lend, $new_rend ) = ( $lend1, $rend1 );
    }
    elsif ($lend2 <= $lend1
        && $rend2 >= $rend1 )
    {
        ( $new_lend, $new_rend ) = ( $lend2, $rend2 );
    }
    elsif ($lend2 <= $lend1
        && $rend2 >= $lend1
        && $rend2 <= $rend1 )
    {
        ( $new_lend, $new_rend ) = ( $lend2, $rend1 );
    }
    elsif ($lend2 >= $lend1
        && $lend2 <= $rend1
        && $rend2 >= $rend1 )
    {
        ( $new_lend, $new_rend ) = ( $lend1, $rend2 );
    }
    return ( $new_lend, $new_rend );
}

sub yank_prot {

# get prot info using yank.
# $prot_r: hash reference. keys are the minimal acc, ie. DB|acc.
# new keys under $$prot_r{$acc} are: name, genus, species, length, seq, other_acc, first_acc
    use strict;
    my ( $prot_r, $log_header, $log_r, $prog, $cdbyank_index_file, $dir,
        $file_root )
      = @_;
    $log_header .= "--sub yank_prot";
    my ( $max_yank, $c1, $a1, @yank_list ) = ( 40, 0, 0, () );
    my $file = $$;

#    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB1 yank\n");
# create list for yank
    for my $acc (
        sort
        keys %$prot_r
      )
    {    ############## sort is not necessary
        $$prot_r{$acc} = {} if ( !defined $$prot_r{$acc} );
        ++$a1;
        if ( $a1 >= $max_yank ) {
            $a1 = 0;
            ++$c1;
            $yank_list[$c1] .= "\"$acc\" ";
        }
        else {
            $yank_list[$c1] .= "\"$acc\" ";
        }
    }

#    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB2 yank\n");
    if ( &yank_on("YANK_HANDLE") or $prog eq "cdbyank" ) {

#	&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB3 yank\n");
        my %yanked;
        for my $list (@yank_list) {
            my @file_cont;
            if ( $dir ne '' and $file_root ne '' ) {
                my $tmp_file = &tmp_file( $dir, $file_root );

                #create tmp_file for yanked sequences
                if ( $prog eq "yank" or $prog eq '' ) {
                    system "echo $list | yank > $dir/$tmp_file";
                }
                elsif ( $prog eq "cdbyank" ) {
                    system
"echo $list | cdbyank $cdbyank_index_file > $dir/$tmp_file";
                }

                #read tmp_file
                open( FH, "$dir/$tmp_file" );
                @file_cont = <FH>;
                close FH;

#get total number of acc and compare it to the entries in tmp_file.  How about some is missing?? not implemented
                $list =~ s/^\s+|\s+$//g;
                my $num_list = scalar( split /\s+/, $list );
                my $count = 0;
                for (@file_cont) {
                    ++$count
                      if ( $_ =~ /^\>/ );
                }
                print
"+++Warning: the number of prot acc in list $num_list is different from that in yank results $count\n"
                  if ( $num_list != $count );

                #                unlink("$dir/$tmp_file");

#	    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB4 yank\n");
            }
            else {
                if (   $prog eq "yank"
                    or $prog eq '' )
                {
                    open( FH,
                        "echo $list | /usr/local/packages/EMBOSS/bin/yank |" );
                }
                elsif ( $prog eq "cdbyank" ) {
                    open( FH,
"echo $list | /usr/local/bin/cdbyank $cdbyank_index_file |"
                    );
                }
                @file_cont = <FH>;
                close FH;
            }
            my $old_acc = '';
            for my $line (@file_cont) {
                chomp $line;
                if ( $line =~ /^>\S+$/ ) {
                    $line =~ s/$/ omni/;
                }    #bug fix for OMNI tag with no whitespace - put it in!
                if ( $line =~ /^>\S+\s+/ ) {

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC1 yank\n");
                    $line =~ s/^>//;
                    my @tmp_acc = split /\cA/, $line;
                    my $curr_acc = '';
                    for my $n (@tmp_acc) {
                        if ( $n =~ /^([^\|]+\|[^\|\/\s]+)/ ) {

#			    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC2 yank\n");
                            my $new_acc = $1;
                            if (   $new_acc ne ''
                                && defined $$prot_r{$new_acc}
                                && $yanked{$new_acc} != 1 )
                            {
                                $yanked{$new_acc} = 1;
                                $curr_acc = $new_acc;

#				&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC3 yank\n");
                                last;
                            }
                        }
                        else {
                            print "***Error: $log_header--wrong acc: $n\n";
                            $$log_r .= "***Error: $log_header--wrong acc: $n\n"
                              if ($log_r);
                        }
                    }

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC4 yank\n");
                    if ( $curr_acc =~ /\S/ ) {

#			&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC5 yank\n");
                        $$prot_r{$curr_acc}{first_acc} = $1
                          if ( $tmp_acc[0] =~ /^(\S+)\s/ );
                        for my $n (@tmp_acc) {
                            $$prot_r{$curr_acc}{other_acc}{$1} = 1
                              if ( $n =~ /^([^\|]+\|[^\|\/\s]+)/
                                && $1 ne $curr_acc );
                        }
                    }
                    else {
                        print
"***Error: $log_header--accession returned by YANK does not match any input accessions which should be in the format of DB|ID\n";
                        $log_r .=
"***Error: $log_header--accession returned by YANK does not match any input accessions which should be in the format of DB|ID\n"
                          if ($log_r);
                    }

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC6 yank\n");
                    if ( $line =~ s/\{([^\}]+)\}// ) {

#			&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC7 yank\n");
                        (
                            $$prot_r{$curr_acc}{genus},
                            $$prot_r{$curr_acc}{species}
                        ) = split( /\s+/, $1, 2 );
                    }
                    $tmp_acc[0] =~ s/\{([^\}]+)\}//;
                    $$prot_r{$curr_acc}{name} = $1
                      if ( $tmp_acc[0] =~ /^\S+\s+(.+?)\s*$/ );
                    if ( $old_acc ne '' ) {

#			&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC8 yank\n");
                        $$prot_r{$old_acc}{seq} =~ s/[^A-Za-z]//g
                          if ( $$prot_r{$old_acc}{seq} =~ /\S/ );
                        $$prot_r{$old_acc}{length} =
                          length( $$prot_r{$old_acc}{seq} )
                          if ( $old_acc ne '' );
                    }
                    $old_acc = $curr_acc;

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC9 yank\n");
                }
                else {
                    $$prot_r{$old_acc}{seq} .= $line
                      if ( $old_acc ne '' );
                }
            }

#	    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB5 yank\n");
            $$prot_r{$old_acc}{seq} =~ s/[^A-Za-z]//g
              if ( $$prot_r{$old_acc}{seq} =~ /\S/ );
            $$prot_r{$old_acc}{length} = length( $$prot_r{$old_acc}{seq} )
              if ( $old_acc ne '' );
        }
    }
    else {
        print "***Error: $log_header--YANK server is probably down\n";
        $$log_r .= "***Error: $log_header--YANK server is probably down\n"
          if ($log_r);
    }

#    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB6 yank\n");
}

sub yank_prot2 {

# get prot info using yank.
# $prot_r: hash reference. keys are the minimal acc, ie. DB|acc.
# new keys under $$prot_r{$acc} are: name, genus, species, length, seq, other_acc, first_acc
    use strict;
    use Cwd;
    use Cwd 'chdir';
    my ( $prot_r, $log_header, $log_r, $prog, $cdbyank_index_file, $dir,
        $file_root )
      = @_;
    $log_header .= "--sub yank_prot2";
    my ( $max_yank, $c1, $a1, @yank_list ) = ( 40, 0, 0, () );
    my $file = $$;

#    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB1 yank\n");
# create list for yank
    my $acc_file;
    my $tmp_dir;
    if ( -d "$ENV{WEBSERVER_TMP}" ) {
        $tmp_dir  = $ENV{WEBSERVER_TMP};
        $acc_file = $ENV{WEBSERVER_TMP} . "/.$$" . ".acc";
    }
    else {
        $tmp_dir  = "/tmp";
        $acc_file = "/tmp/" . ".$$" . ".acc";
    }
    my $pwd       = cwd();
    my $acc_count = 0;
    open( ACC, ">$acc_file" )
      || die "cannot open $acc_file to write: $!\n";
    for my $acc (
        sort
        keys %$prot_r
      )
    {    ############## sort is not necessary
        $$prot_r{$acc} = {}
          if ( !defined $$prot_r{$acc} );
        if ( $acc !~ /^gi/ ) {
            print ACC "$acc\n";
        }
        $acc_count++;
        ++$a1;
        if ( $a1 >= $max_yank ) {
            $a1 = 0;
            ++$c1;
            $yank_list[$c1] .= "\"$acc\" ";
        }
        else {
            $yank_list[$c1] .= "\"$acc\" ";
        }
    }
    close ACC;
    system "chmod 777 $acc_file";

#    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB2 yank\n");
#    if(&yank_on("YANK_HANDLE") or $prog eq "cdbyank") {
#	&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB3 yank\n");
    my %yanked;
    my @final_file;
    chdir $tmp_dir
      || die "cannot change to dir $tmp_dir\n";
    if (    $dir ne ''
        and $file_root ne '' )
    {
        my $tmp_file = &tmp_file( $dir, $file_root );

        #create tmp_file for yanked sequences
        if (   $prog eq "yank_panda"
            or $prog eq '' )
        {
            system "yank_panda --group=pep --f=$acc_file > $dir/$tmp_file";
            warn "$dir/$tmp_file\n";
            open( FH, "$dir/$tmp_file" );
            @final_file = <FH>;
            close FH;
            my $count = 0;
            for (@final_file) {
                ++$count
                  if ( $_ =~ /^\>/ );
            }
            print
"+++Warning: the number of prot acc in list $acc_count is different from that in yank results $count\n"
              if ( $acc_count != $count );
            unlink("$dir/$tmp_file");
        }
        elsif ( $prog eq "cdbyank" ) {
            for my $list (@yank_list) {
                my @file_cont;
                system
                  "echo $list | cdbyank $cdbyank_index_file > $dir/$tmp_file";
                open( FH, "$dir/$tmp_file" );
                @file_cont = <FH>;
                close FH;
                $list =~ s/^\s+|\s+$//g;
                my $num_list = scalar( split /\s+/, $list );
                my $count = 0;

                for (@file_cont) {
                    ++$count
                      if ( $_ =~ /^\>/ );
                }
                print
"+++Warning: the number of prot acc in list $num_list is different from that in yank results $count\n"
                  if ( $num_list != $count );
                unlink("$dir/$tmp_file");
                push( @final_file, @file_cont );
            }
        }
    }
    else {
        if (   $prog eq "yank_panda"
            or $prog eq '' )
        {
            open( FH, "/usr/local/common/yank_panda --f=$acc_file |" );
            @final_file = <FH>;
            close FH;
        }
        elsif ( $prog eq "cdbyank" ) {
            for my $list (@yank_list) {
                my @file_cont;
                open( FH,
                    "echo $list | /usr/local/bin/cdbyank $cdbyank_index_file |"
                );
                @file_cont = <FH>;
                close FH;
                push( @final_file, @file_cont );
            }
        }
    }
    chdir $pwd
      || die "cannot change to dir $pwd\n";
    my $old_acc = '';
    for my $line (@final_file) {
        chomp $line;

        #	    print "$line\n";
        if ( $line =~ /^>\S+$/ ) {
            $line =~ s/$/ omni/;
        }    #bug fix for OMNI tag with no whitespace - put it in!
        if ( $line =~ /^>\S+\s+/ ) {

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC1 yank\n");
            $line =~ s/^>//;
            my @tmp_acc;
            if ( $prog eq "cdbyank" ) {
                @tmp_acc = split /\cA/, $line;
            }
            elsif ($prog eq "yank_panda"
                or $prog eq "" )
            {
                @tmp_acc = split( /\s\^\|\^/, $line );
            }
            my $curr_acc = '';
            for my $n (@tmp_acc) {
                if ( $n =~ /^([^\|]+\|[^\|\/\s]+)/ ) {

#			    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC2 yank\n");
                    my $new_acc = $1;
                    if (   $new_acc ne ''
                        && defined $$prot_r{$new_acc}
                        && $yanked{$new_acc} != 1 )
                    {
                        $yanked{$new_acc} = 1;
                        $curr_acc = $new_acc;

#				&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC3 yank\n");
                        last;
                    }
                }
                else {
                    print "***Error: $log_header--wrong acc: $n\n";
                    $$log_r .= "***Error: $log_header--wrong acc: $n\n"
                      if ($log_r);
                }
            }

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC4 yank\n");
            if ( $curr_acc =~ /\S/ ) {

#			&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC5 yank\n");
                $$prot_r{$curr_acc}{first_acc} = $1
                  if ( $tmp_acc[0] =~ /^(\S+)\s/ );
                for my $n (@tmp_acc) {
                    $$prot_r{$curr_acc}{other_acc}{$1} = 1
                      if ( $n =~ /^([^\|]+\|[^\|\/\s]+)/
                        && $1 ne $curr_acc );
                }
            }
            else {
                print
"***Error: $log_header--accession returned by YANK does not match any input accessions which should be in the format of DB|ID\n";
                $log_r .=
"***Error: $log_header--accession returned by YANK does not match any input accessions which should be in the format of DB|ID\n"
                  if ($log_r);
            }

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC6 yank\n");
            if ( $line =~ s/\{([^\}]+)\}// ) {

#			&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC7 yank\n");
                ( $$prot_r{$curr_acc}{genus}, $$prot_r{$curr_acc}{species} ) =
                  split( /\s+/, $1, 2 );
            }
            $tmp_acc[0] =~ s/\{([^\}]+)\}//;
            $$prot_r{$curr_acc}{name} = $1
              if ( $tmp_acc[0] =~ /^\S+\s+(.+?)\s*$/ );
            if ( $old_acc ne '' ) {

#			&print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC8 yank\n");
                $$prot_r{$old_acc}{seq} =~ s/[^A-Za-z]//g
                  if ( $$prot_r{$old_acc}{seq} =~ /\S/ );
                $$prot_r{$old_acc}{length} = length( $$prot_r{$old_acc}{seq} )
                  if ( $old_acc ne '' );
            }
            $old_acc = $curr_acc;

#		    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "CC9 yank\n");
        }
        else {
            $$prot_r{$old_acc}{seq} .= $line
              if ( $old_acc ne '' );
        }
    }

#	    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB5 yank\n");
    $$prot_r{$old_acc}{seq} =~ s/[^A-Za-z]//g
      if ( $$prot_r{$old_acc}{seq} =~ /\S/ );
    $$prot_r{$old_acc}{length} = length( $$prot_r{$old_acc}{seq} )
      if ( $old_acc ne '' );

#    } else {
#	print "***Error: $log_header--YANK server is probably down\n";
#	$$log_r .= "***Error: $log_header--YANK server is probably down\n" if($log_r);
#    }
    unlink("$acc_file");

#    &print_cyc("/usr/local/projects/OGC/egad_test/ogc_cyc_log/$file", "BB6 yank\n");
}

sub yank_on {

    # test whether yank server is on or off
    use strict;
    my ($handle) = @_;
    $handle = (
          $handle ne ''
        ? $handle
        : "YANK_HANDLE"
    );
    my $result = '';
    {
        no strict "refs";
        open( $handle,
            "/usr/local/packages/EMBOSS/bin/yank -t -a \"EGAD\|11\" | " );
        while (<$handle>) {
            $result .= $_;
        }
        close $handle;
    }
    return (
        $result =~ /^>SP/
        ? 1
        : 0
    );
}

sub get_INFO_from_file {

    # read file and save data in hash reference $data_r
    use strict;
    my ( $data_r, $file, $log_r ) = @_;
    if ( -s $file ) {
        open( FH, "$file" );
        while ( my $line = <FH> ) {
            chomp $line;
            $line =~ s/^\s*|\s*$//g;

            #single line entry
            if ( $line =~ /^AC\s+(.+)/ ) {
                $log_r .= "***Error: extra AC field found\n"
                  if ( $$data_r{'hmm_acc'} ne '' );
                $$data_r{'hmm_acc'} = $1;
            }
            elsif ( $line =~ /^ID\s+(.+)/ ) {
                $log_r .= "***Error: extra ID field found\n"
                  if ( $$data_r{'hmm_name'} ne '' );
                $$data_r{'hmm_name'} = $1;
            }
            elsif ( $line =~ /^DE\s+(.+)/ ) {
                $log_r .= "***Error: extra DE field found\n"
                  if ( $$data_r{'hmm_com_name'} ne '' );
                $$data_r{'hmm_com_name'} = $1;
            }
            elsif ( $line =~ /^AU\s+(.+)/ ) {
                $log_r .= "***Error: extra AU field found\n"
                  if ( $$data_r{'author'} ne '' );
                $$data_r{'author'} = $1;
            }
            elsif ( $line =~ /^RI\s+(.+)/ ) {
                $log_r .= "***Error: extra RI field found\n"
                  if ( $$data_r{'role_id'} ne '' );
                $$data_r{'role_id'} = $1;
            }
            elsif ( $line =~ /^EC\s+(.+)/ ) {
                $log_r .= "***Error: extra EC field found\n"
                  if ( $$data_r{'ec_num'} ne '' );
                $$data_r{'ec_num'} = $1;
            }
            elsif ( $line =~ /^TR\s+(.+)/ ) {
                $log_r .= "***Error: extra TR field found\n"
                  if ( $$data_r{'tc_num'} ne '' );
                $$data_r{'tc_num'} = $1;
            }
            elsif ( $line =~ /^AL\s+(.+)/ ) {
                $log_r .= "***Error: extra AL field found\n"
                  if ( $$data_r{'method_seed'} ne '' );
                $$data_r{'method_seed'} = $1;
            }
            elsif ( $line =~ /^IT\s+(.+)/ ) {
                $log_r .= "***Error: extra IT field found\n"
                  if ( $$data_r{'iso_type'} ne '' );
                $$data_r{'iso_type'} = $1;
            }
            elsif ( $line =~ /^TP\s+(.+)/ ) {
                $log_r .= "***Error: extra TP field found\n"
                  if ( $$data_r{'hmm_type'} ne '' );
                $$data_r{'hmm_type'} = $1;
            }
            elsif ( $line =~ /^GA\s+(.+)/ ) {
                $log_r .= "***Error: extra GA field found\n"
                  if ( $$data_r{'search_prog'} ne '' );
                $$data_r{'search_prog'} = $1;
            }

#	    elsif($line =~ /^SA\s+(.+)/) { $log_r .= "***Error: extra SA field found\n" if($$data_r{'struct_acc'}     ne ''); $$data_r{'struct_acc'}     = $1;}
            elsif ( $line =~ /^GS\s+(.+)/ ) {
                $log_r .= "***Error: extra GS field found\n"
                  if ( $$data_r{'gene_sym'} ne '' );
                $$data_r{'gene_sym'} = $1;
            }
            elsif ( $line =~ /^EN\s+(.+)/ ) {
                $log_r .= "***Error: extra EN field found\n"
                  if ( $$data_r{'expanded_name'} ne '' );
                $$data_r{'expanded_name'} = $1;
            }
            elsif ( $line =~ /^SO\s+(.+)/ ) {
                $log_r .= "***Error: extra SO field found\n"
                  if ( $$data_r{'search_options'} ne '' );
                $$data_r{'search_options'} = $1;
            }
            elsif ( $line =~ /^RH\s+(.+)/ ) {
                $log_r .= "***Error: extra RH field found\n"
                  if ( $$data_r{'related_hmm'} ne '' );
                $$data_r{'related_hmm'} = $1;
            }
            elsif ( $line =~
                /^NC\s+(-?(\d+\.?\d*|\d*\.\d+))(\s+(-?(\d+\.?\d*|\d*\.\d+)))?/ )
            {
                $log_r .= "***Error: extra NC field found\n"
                  if ( $$data_r{'noise_cutoff'} ne ''
                    or $$data_r{'noise_cutoff2'} ne '' );
                $$data_r{'noise_cutoff'}  = $1;
                $$data_r{'noise_cutoff2'} = (
                      $4 ne ''
                    ? $4
                    : $1
                );
            }
            elsif ( $line =~
                /^TC\s+(-?(\d+\.?\d*|\d*\.\d+))(\s+(-?(\d+\.?\d*|\d*\.\d+)))?/ )
            {
                $log_r .= "***Error: extra TC field found\n"
                  if ( $$data_r{'trusted_cutoff'} ne ''
                    or $$data_r{'trusted_cutoff2'} ne '' );
                $$data_r{'trusted_cutoff'}  = $1;
                $$data_r{'trusted_cutoff2'} = (
                      $4 ne ''
                    ? $4
                    : $1
                );
            }
            elsif ( $line =~ /^TIGR\s+hmm_id\s+(\d+)/ ) {
                $log_r .= "***Error: extra hmm_id field found\n"
                  if ( $$data_r{hmm_id} ne '' );
                $$data_r{hmm_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_align_id\s+(\d+)/ ) {
                $log_r .= "***Error: extra seed_align_id field found\n"
                  if ( $$data_r{seed_align_id} ne '' );
                $$data_r{seed_align_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_align_id\s+(\d+)/ ) {
                $log_r .= "***Error: extra full_align_id field found\n"
                  if ( $$data_r{full_align_id} ne '' );
                $$data_r{full_align_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+hmm_file\s+(.+)\s*$/ ) {
                $log_r .= "***Error: extra hmm_file field found\n"
                  if ( $$data_r{hmm_file} ne '' );
                $$data_r{hmm_file} = $1;
            }
            elsif ( $line =~ /^TIGR\s+frag_file\s+(.+)\s*$/ ) {
                $log_r .= "***Error: extra frag_file field found\n"
                  if ( $$data_r{frag_file} ne '' );
                $$data_r{frag_file} = $1;
            }
            elsif ( $line =~ /^TIGR\s+build_file\s+(.+)\s*$/ ) {
                $log_r .= "***Error: extra build_file field found\n"
                  if ( $$data_r{build_file} ne '' );
                $$data_r{build_file} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_file\s+(.+)\s*$/ ) {
                $log_r .= "***Error: extra seed_file field found\n"
                  if ( $$data_r{seed_file} ne '' );
                $$data_r{seed_file} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_file\s+(.+)\s*$/ ) {
                $log_r .= "***Error: extra full_file field found\n"
                  if ( $$data_r{full_file} ne '' );
                $$data_r{full_file} = $1;
            }
            elsif ( $line =~ /^TIGR\s+iso_id\s+(\d+)/ ) {
                $log_r .= "***Error: extra iso_id field found\n"
                  if ( $$data_r{iso_id} ne '' );
                $$data_r{iso_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+is_current\s+(-?\d+)/ ) {
                $log_r .= "***Error: extra is_current field found\n"
                  if ( $$data_r{is_current} ne '' );
                $$data_r{is_current} = $1;
            }
            elsif ( $line =~ /^TIGR\s+num_seed\s+(\d+)/ ) {
                $log_r .= "***Error: extra num_seed field found\n"
                  if ( $$data_r{num_seed} ne '' );
                $$data_r{num_seed} = $1;
            }
            elsif ( $line =~ /^TIGR\s+prior\s+(.+)/ ) {
                $log_r .= "***Error: extra prior field found\n"
                  if ( $$data_r{prior} ne '' );
                $$data_r{prior} = $1;
            }
            elsif ( $line =~ /^TIGR\s+avg_score\s+(.+)/ ) {
                $log_r .= "***Error: extra avg_score field found\n"
                  if ( $$data_r{avg_score} ne '' );
                $$data_r{avg_score} = $1;
            }
            elsif ( $line =~ /^TIGR\s+std_dev\s+(.+)/ ) {
                $log_r .= "***Error: extra std_dev field found\n"
                  if ( $$data_r{std_dev} ne '' );
                $$data_r{std_dev} = $1;
            }
            elsif ( $line =~ /^TIGR\s+min_score\s+(.+)/ ) {
                $log_r .= "***Error: extra min_score field found\n"
                  if ( $$data_r{min_score} ne '' );
                $$data_r{min_score} = $1;
            }
            elsif ( $line =~ /^TIGR\s+max_score\s+(.+)/ ) {
                $log_r .= "***Error: extra max_score field found\n"
                  if ( $$data_r{max_score} ne '' );
                $$data_r{max_score} = $1;
            }
            elsif ( $line =~ /^TIGR\s+hmm_len\s+(\d+)/ ) {
                $log_r .= "***Error: extra hmm_len field found\n"
                  if ( $$data_r{hmm_len} ne '' );
                $$data_r{hmm_len} = $1;
            }
            elsif ( $line =~ /^TIGR\s+iso_acc\s+(.+)/ ) {
                $log_r .= "***Error: extra iso_acc field found\n"
                  if ( $$data_r{iso_acc} ne '' );
                $$data_r{iso_acc} = $1;
            }
            elsif ( $line =~ /^TIGR\s+iso_iso_type\s+(.+)/ ) {
                $log_r .= "***Error: extra iso_iso_type field found\n"
                  if ( $$data_r{iso_iso_type} ne '' );
                $$data_r{iso_iso_type} = $1;
            }
            elsif ( $line =~ /^TIGR\s+iso_method\s+(.+)/ ) {
                $log_r .= "***Error: extra iso_method field found\n"
                  if ( $$data_r{iso_method} ne '' );
                $$data_r{iso_method} = $1;
            }
            elsif ( $line =~ /^TIGR\s+assigned_by\s+(.+)/ ) {
                $log_r .= "***Error: extra assigned_by field found\n"
                  if ( $$data_r{assigned_by} ne '' );
                $$data_r{assigned_by} = $1;
            }
            elsif ( $line =~ /^TIGR\s+class_name\s+(.+)/ ) {
                $log_r .= "***Error: extra class_name field found\n"
                  if ( $$data_r{class_name} ne '' );
                $$data_r{class_name} = $1;
            }
            elsif ( $line =~ /^TIGR\s+db_type\s+(.+)/ ) {
                $log_r .= "***Error: extra db_type field found\n"
                  if ( $$data_r{db_type} ne '' );
                $$data_r{db_type} = $1;
            }
            elsif ( $line =~ /^TIGR\s+iso_ref_id\s+(\d+)/ ) {
                $log_r .= "***Error: extra iso_ref_id field found\n"
                  if ( $$data_r{iso_ref_id} ne '' );
                $$data_r{iso_ref_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+user_class\s+(.+)/ ) {
                $log_r .= "***Error: extra user_class field found\n"
                  if ( $$data_r{user_class} ne '' );
                $$data_r{user_class} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seq_type\s+(.+)/ ) {
                $log_r .= "***Error: extra seq_type field found\n"
                  if ( $$data_r{seq_type} ne '' );
                $$data_r{seq_type} = $1;
            }
            elsif ( $line =~ /^TIGR\s+iso_comment\s+(.+)/ ) {
                $log_r .= "***Error: extra iso_comment field found\n"
                  if ( $$data_r{iso_comment} ne '' );
                $$data_r{iso_comment} = $1;
            }
            elsif ( $line =~ /^TIGR\s+alignment_iso_id\s+(.+)/ ) {
                $log_r .= "***Error: extra alignment_iso_id field found\n"
                  if ( $$data_r{alignment_iso_id} ne '' );
                $$data_r{alignment_iso_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_align_name\s+(.+)/ ) {
                $log_r .= "***Error: extra seed_align_name field found\n"
                  if ( $$data_r{seed_align_name} ne '' );
                $$data_r{seed_align_name} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_user_class\s+(.+)/ ) {
                $log_r .= "***Error: extra seed_user_class field found\n"
                  if ( $$data_r{seed_user_class} ne '' );
                $$data_r{seed_user_class} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_comment\s+(.+)/ ) {
                $log_r .= "***Error: extra seed_comment field found\n"
                  if ( $$data_r{seed_comment} ne '' );
                $$data_r{seed_comment} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_ref_id\s+(.+)/ ) {
                $log_r .= "***Error: extra seed_ref_id field found\n"
                  if ( $$data_r{seed_ref_id} ne '' );
                $$data_r{seed_ref_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_align_type\s+(.+)/ ) {
                $log_r .= "***Error: extra seed_align_type field found\n"
                  if ( $$data_r{seed_align_type} ne '' );
                $$data_r{seed_align_type} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_align_name\s+(.+)/ ) {
                $log_r .= "***Error: extra full_align_name field found\n"
                  if ( $$data_r{full_align_name} ne '' );
                $$data_r{full_align_name} = $1;
            }
            elsif ( $line =~ /^TIGR\s+method_full\s+(.+)/ ) {
                $log_r .= "***Error: extra method_full field found\n"
                  if ( $$data_r{method_full} ne '' );
                $$data_r{method_full} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_user_class\s+(.+)/ ) {
                $log_r .= "***Error: extra full_user_class field found\n"
                  if ( $$data_r{full_user_class} ne '' );
                $$data_r{full_user_class} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_comment\s+(.+)/ ) {
                $log_r .= "***Error: extra full_comment field found\n"
                  if ( $$data_r{full_comment} ne '' );
                $$data_r{full_comment} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_ref_id\s+(.+)/ ) {
                $log_r .= "***Error: extra full_ref_id field found\n"
                  if ( $$data_r{full_ref_id} ne '' );
                $$data_r{full_ref_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_align_type\s+(.+)/ ) {
                $log_r .= "***Error: extra full_align_type field found\n"
                  if ( $$data_r{full_align_type} ne '' );
                $$data_r{full_align_type} = $1;
            }
            elsif ( $line =~ /^TIGR\s+hmm_method_id\s+(.+)/ ) {
                $log_r .= "***Error: extra hmm_method_id field found\n"
                  if ( $$data_r{hmm_method_id} ne '' );
                $$data_r{hmm_method_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+hmm_frag_method_id\s+(.+)/ ) {
                $log_r .= "***Error: extra hmm_frag_method_id field found\n"
                  if ( $$data_r{hmm_frag_method_id} ne '' );
                $$data_r{hmm_frag_method_id} = $1;
            }
            elsif ( $line =~ /^TIGR\s+seed_align_assigned_by\s+(.+)/ ) {
                $log_r .= "***Error: extra seed_align_assigned_by field found\n"
                  if ( $$data_r{seed_align_assigned_by} ne '' );
                $$data_r{seed_align_assigned_by} = $1;
            }
            elsif ( $line =~ /^TIGR\s+full_align_assigned_by\s+(.+)/ ) {
                $log_r .= "***Error: extra full_align_assigned_by field found\n"
                  if ( $$data_r{full_align_assigned_by} ne '' );
                $$data_r{full_align_assigned_by} = $1;
            }
            elsif ( $line =~ /^TIGR\s+nomen_check\s+(\d)\s*$/ ) {
                $log_r .= "***Error: extra nomen_check field found\n"
                  if ( $$data_r{nomen_check} ne '' );
                $$data_r{nomen_check} = $1;
            }

#	    elsif($line =~ /^TIGR\s+X\s+(.+)\s*$/) { $log_r .= "***Error: extra X field found\n" if($$data_r{X} ne ''); $$data_r{X} = $1;}
#multiple line entry
            elsif ( $line =~ /^CC\s+(.+)$/ ) {
                $$data_r{hmm_comment} .= "$1\n";
            }
            elsif ( $line =~ /^PC\s+(.+)$/ ) {
                $$data_r{private} .= "$1\n";
            }

            #reference collection
            elsif ( $line =~ /^TIGR\s+reference\s+(.+)$/ ) {
                $$data_r{reference} .= "$1\n";
            }
            elsif ( $line =~ /^DC\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^DR\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^RC\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^RN\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^RT\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^RM\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^RL\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^RA\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^KW\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^PI\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^BD\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^SE\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^SQ\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^A2\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^NS\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^PR\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^BM\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            elsif ( $line =~ /^XD\s+/ ) {
                $$data_r{reference} .= "$line\n";
            }
            else {
                $$data_r{reference} .= "$line\n";
                $log_r .=
"***Error: unknown field $line found in file $file, save it to reference field\n"
                  if ($log_r);
            }
        }
    }

    #fixing reference
    chomp $$data_r{reference};
    chomp $$data_r{hmm_comment};
    chomp $$data_r{private};
    return 0;
}

sub get_BUILD_from_file {
    use strict;
    my ( $data_r, $build_file, $log_r ) = @_;
    if ( -s $build_file
        && open( FH, "$build_file" ) )
    {
        while (<FH>) {
            $$data_r{hmm_build} .= $_;
            $$data_r{num_seed} = $1
              if ( $_ =~ /^Number\s+of\s+sequences:\s+(\d+)\s*$/ );
            $$data_r{avg_score} = $1
              if ( $_ =~ /^Average\s+score:\s+(-?\d*\.?\d*)\s+bits$/ );
            $$data_r{std_dev} = $1
              if ( $_ =~ /^Std\.\s+deviation:\s+(-?\d*\.?\d*)\s+bits$/ );
            $$data_r{min_score} = $1
              if ( $_ =~ /^Minimum\s+score:\s+(-?\d*\.?\d*)\s+bits$/ );
            $$data_r{max_score} = $1
              if ( $_ =~ /^Maximum\s+score:\s+(-?\d*\.?\d*)\s+bits$/ );
            $$data_r{hmm_len_check} = $1
              if ( $_ =~
                /^Constructed\s+a\s+profile\s+HMM\s+\(length\s+(\d+)\)\s*$/ );
            $$data_r{prior} = $1
              if ( $_ =~ /^Prior strategy:\s+(.+)$/ );
        }
        close FH;
    }
}

sub get_HMM_from_file {
    use strict;
    my ( $data_r, $hmm_file, $log_r ) = @_;
    if ( -s $hmm_file
        && open( FH, "$hmm_file" ) )
    {
        while (<FH>) {
            $$data_r{hmm} .= $_;
            $$data_r{hmm_len} = $1
              if ( $_ =~ /^LENG\s+(\d+)\s*$/ );
        }
        close FH;
    }
    chomp $$data_r{hmm};
}

sub get_FRAG_from_file {
    use strict;
    my ( $data_r, $frag_file, $log_r ) = @_;
    if ( -s $frag_file
        && open( FH, "$frag_file" ) )
    {
        while (<FH>) {
            $$data_r{hmm_frag} .= $_;
        }
        close FH;
    }
    chomp $$data_r{hmm_frag};
}

sub get_SEED_from_file {
    use strict;
    my ( $data_r, $seed_file, $log_r ) = @_;
    if (   -s $seed_file
        && -r $seed_file )
    {
        my %seed;
        &read_sequence( \%seed, $seed_file, '', '', '' );
        $$data_r{seed} = \%seed;
    }
    else {
        $log_r .= "***Error: can not read SEED file\n"
          if ($log_r);
    }
}

sub get_FULL_from_file {
    use strict;
    my ( $data_r, $full_file, $log_r ) = @_;
    if (   -s $full_file
        && -r $full_file )
    {
        my %full;
        &read_sequence( \%full, $full_file, '', '', '' );
        $$data_r{full} = \%full;
    }
}

#    $name = $1 if($name =~ /^([^\.]+)/);
#    my($hmm_file, $hmm_frag_file, $build_file, $seed_file, $full_file);
#    $hmm_file      = "$dir/$name.HMM"   if(-s "$dir/$name.HMM");
#    $hmm_file      = "$dir/$name.hmm"   if(-s "$dir/$name.hmm"   && $hmm_file eq '');
#    $hmm_frag_file = "$dir/$name.FRAG"  if(-s "$dir/$name.FRAG");
#    $hmm_frag_file = "$dir/$name.frag"  if(-s "$dir/$name.frag"  && $hmm_frag_file eq '');
#    $build_file    = "$dir/$name.BUILD" if(-s "$dir/$name.BUILD");
#    $build_file    = "$dir/$name.build" if(-s "$dir/$name.build" && $build_file eq '');
#    $seed_file     = "$dir/$name.SEED"  if(-s "$dir/$name.SEED");
#    $seed_file     = "$dir/$name.seed"  if(-s "$dir/$name.seed"  && $seed_file eq '');
#    $full_file     = "$dir/$name.FULL"  if(-s "$dir/$name.FULL");
#    $full_file     = "$dir/$name.full"  if(-s "$dir/$name.full"  && $full_file eq '');
sub get_INFO_from_db {

    #create INFO information using data in DB
    use strict;
    my ( $db_proc, $db, $query_id, $query_type, $data_r ) = @_;
    my ( $query, $ref1, $ref2, $ref3, $ref4, $ref5, $ref6 ) =
      ( '', '', '', '', '', '', '' );
    if (   $query_type eq 'hmm_acc'
        && $query_id ne '' )
    {
        $ref1 = &do_query( $db_proc,
"select * from hmm2 where hmm_acc = \"$query_id\" and is_current = 1"
        );
        return
          if ( $$ref1[0]{iso_id} eq '' );
        $ref2 = &do_query( $db_proc,
            "select * from isology where iso_id = $$ref1[0]{iso_id}" );
        $ref3 = &do_query( $db_proc,
"select * from alignment where iso_id = $$ref1[0]{iso_id} and align_type=\"seed\""
        );
        $ref4 = &do_query( $db_proc,
"select * from alignment where iso_id = $$ref1[0]{iso_id} and align_type=\"all\""
        );
        $ref5 = &do_query( $db_proc,
            "select * from hmm_role_link where hmm_acc=\"$query_id\"" );
        $ref6 = &do_query( $db_proc,
            "select * from hmm_go_link where hmm_acc=\"$query_id\"" );
    }
    elsif ($query_type eq 'iso_id'
        && $query_id ne '' )
    {
        $ref1 = &do_query( $db_proc,
            "select * from hmm2 where iso_id = $query_id and is_current = 1" );
        my $hmm_acc = $$ref1[0]{hmm_acc};
        $ref2 = &do_query( $db_proc,
            "select * from isology where iso_id = $query_id" );
        $ref3 = &do_query( $db_proc,
"select * from alignment where iso_id = $query_id and align_type=\"seed\""
        );
        $ref4 = &do_query( $db_proc,
"select * from alignment where iso_id = $query_id and align_type=\"all\""
        );
        $ref5 =
          &do_query( $db_proc,
            "select * from hmm_role_link where hmm_acc=\"$hmm_acc\"" )
          if ( $hmm_acc =~ /\S/ );
        $ref6 =
          &do_query( $db_proc,
            "select * from hmm_go_link where hmm_acc=\"$hmm_acc\"" )
          if ( $hmm_acc =~ /\S/ );
    }
    else {
        return;
    }
    $$data_r{hmm_id} = $$ref1[0]{id}
      if ( $$ref1[0]{id} ne '' );
    $$data_r{hmm_acc} = $$ref1[0]{hmm_acc}
      if ( $$ref1[0]{hmm_acc} ne '' );
    $$data_r{hmm_type} = $$ref1[0]{hmm_type}
      if ( $$ref1[0]{hmm_type} ne '' );
    $$data_r{hmm_name} = $$ref1[0]{hmm_name}
      if ( $$ref1[0]{hmm_name} ne '' );
    $$data_r{hmm_com_name} = $$ref1[0]{hmm_com_name}
      if ( $$ref1[0]{hmm_com_name} ne '' );
    $$data_r{hmm} = $$ref1[0]{hmm}
      if ( $$ref1[0]{hmm} ne '' );
    $$data_r{hmm_len} = $$ref1[0]{hmm_len}
      if ( $$ref1[0]{hmm_len} ne '' );
    $$data_r{hmm_seed} = $$ref1[0]{hmm_seed}
      if ( $$ref1[0]{hmm_seed} ne '' );
    $$data_r{hmm_frag} = $$ref1[0]{hmm_frag}
      if ( $$ref1[0]{hmm_frag} ne '' );
    $$data_r{iso_id} = $$ref1[0]{iso_id}
      if ( $$ref1[0]{iso_id} ne '' );
    $$data_r{search_prog} = $$ref1[0]{search_prog}
      if ( $$ref1[0]{search_prog} ne '' );
    $$data_r{search_options} = $$ref1[0]{search_options}
      if ( $$ref1[0]{search_options} ne '' );
    $$data_r{hmm_comment} = $$ref1[0]{hmm_comment}
      if ( $$ref1[0]{hmm_comment} ne '' );
    $$data_r{is_current} = $$ref1[0]{is_current}
      if ( $$ref1[0]{is_current} ne '' );
    $$data_r{author} = $$ref1[0]{author}
      if ( $$ref1[0]{author} ne '' );
    $$data_r{entry_date} = $$ref1[0]{entry_date}
      if ( $$ref1[0]{entry_date} ne '' );
    $$data_r{prior} = $$ref1[0]{prior}
      if ( $$ref1[0]{prior} ne '' );
    $$data_r{num_seed} = $$ref1[0]{num_seed}
      if ( $$ref1[0]{num_seed} ne '' );
    $$data_r{avg_score} = $$ref1[0]{avg_score}
      if ( $$ref1[0]{avg_score} ne '' );
    $$data_r{std_dev} = $$ref1[0]{std_dev}
      if ( $$ref1[0]{std_dev} ne '' );
    $$data_r{min_score} = $$ref1[0]{min_score}
      if ( $$ref1[0]{min_score} ne '' );
    $$data_r{max_score} = $$ref1[0]{max_score}
      if ( $$ref1[0]{max_score} ne '' );
    $$data_r{role_id} = $$ref1[0]{role_id}
      if ( $$ref1[0]{role_id} ne '' );
    $$data_r{ec_num} = $$ref1[0]{ec_num}
      if ( $$ref1[0]{ec_num} ne '' );
    $$data_r{tc_num} = $$ref1[0]{tc_num}
      if ( $$ref1[0]{tc_num} ne '' );

#    $$data_r{struct_acc} =      $$ref1[0]{struct_acc}      if($$ref1[0]{struct_acc} ne '');
    $$data_r{mod_date} = $$ref1[0]{mod_date}
      if ( $$ref1[0]{mod_date} ne '' );
    $$data_r{iso_type} = $$ref1[0]{iso_type}
      if ( $$ref1[0]{iso_type} ne '' );
    $$data_r{reference} = $$ref1[0]{reference}
      if ( $$ref1[0]{reference} ne '' );
    $$data_r{noise_cutoff} = $$ref1[0]{noise_cutoff}
      if ( $$ref1[0]{noise_cutoff} ne '' );
    $$data_r{trusted_cutoff} = $$ref1[0]{trusted_cutoff}
      if ( $$ref1[0]{trusted_cutoff} ne '' );
    $$data_r{noise_cutoff2} = $$ref1[0]{noise_cutoff2}
      if ( $$ref1[0]{noise_cutoff2} ne '' );
    $$data_r{trusted_cutoff2} = $$ref1[0]{trusted_cutoff2}
      if ( $$ref1[0]{trusted_cutoff2} ne '' );
    $$data_r{hmm_mod_date} = $$ref1[0]{hmm_mod_date}
      if ( $$ref1[0]{hmm_mod_date} ne '' );
    $$data_r{hmm_build} = $$ref1[0]{hmm_build}
      if ( $$ref1[0]{hmm_build} ne '' );
    $$data_r{private} = $$ref1[0]{private}
      if ( $$ref1[0]{private} ne '' );
    $$data_r{gene_sym} = $$ref1[0]{gene_sym}
      if ( $$ref1[0]{gene_sym} ne '' );
    $$data_r{expanded_name} = $$ref1[0]{expanded_name}
      if ( $$ref1[0]{expanded_name} ne '' );
    $$data_r{hmm_method_id} = $$ref1[0]{hmm_method_id}
      if ( $$ref1[0]{hmm_method_id} ne '' );
    $$data_r{hmm_frag_method_id} = $$ref1[0]{hmm_frag_method_id}
      if ( $$ref1[0]{hmm_frag_method_id} ne '' );
    $$data_r{method_seed} = $$ref1[0]{method_seed}
      if ( $$ref1[0]{method_seed} ne '' );
    $$data_r{nomen_check} = $$ref1[0]{nomen_check}
      if ( $$ref1[0]{nomen_check} ne '' );
    $$data_r{SAC} = $1
      if (
        $$data_r{hmm_build} =~ /\nSearch algorithm configuration:\s*(.+)\n/ );

    # data in isology table
    $$data_r{iso_iso_id} = $$ref2[0]{iso_id}
      if ( $$ref2[0]{iso_id} ne '' );
    $$data_r{class_name} = $$ref2[0]{class_name}
      if ( $$ref2[0]{class_name} ne '' );
    $$data_r{seq_type} = $$ref2[0]{seq_type}
      if ( $$ref2[0]{seq_type} ne '' );
    $$data_r{iso_iso_type} = $$ref2[0]{iso_type}
      if ( $$ref2[0]{iso_type} ne '' );
    $$data_r{iso_method} = $$ref2[0]{method}
      if ( $$ref2[0]{method} ne '' );
    $$data_r{iso_comment} = $$ref2[0]{comment}
      if ( $$ref2[0]{comment} ne '' );
    $$data_r{iso_ref_id} = $$ref2[0]{iso_ref_id}
      if ( $$ref2[0]{iso_ref_id} ne '' );
    $$data_r{assigned_by} = $$ref2[0]{assigned_by}
      if ( $$ref2[0]{assigned_by} ne '' );
    $$data_r{user_class} = $$ref2[0]{user_class}
      if ( $$ref2[0]{user_class} ne '' );
    $$data_r{iso_acc} = $$ref2[0]{iso_acc}
      if ( $$ref2[0]{iso_acc} ne '' );
    $$data_r{db_type} = $$ref2[0]{db_type}
      if ( $$ref2[0]{db_type} ne '' );

    # data in alignment
    $$data_r{seed_align_id} = $$ref3[0]{align_id}
      if ( $$ref3[0]{align_id} ne '' );
    $$data_r{seed_align_name} = $$ref3[0]{align_name}
      if ( $$ref3[0]{align_name} ne '' );
    $$data_r{seed_user_class} = $$ref3[0]{user_class}
      if ( $$ref3[0]{user_class} ne '' );
    $$data_r{seed_comment} = $$ref3[0]{comment}
      if ( $$ref3[0]{comment} ne '' );
    $$data_r{seed_ref_id} = $$ref3[0]{ref_id}
      if ( $$ref3[0]{ref_id} ne '' );
    $$data_r{seed_align_assigned_by} = $$ref3[0]{assigned_by}
      if ( $$ref3[0]{assigned_by} ne '' );
    $$data_r{seed_align_type} = $$ref3[0]{align_type}
      if ( $$ref3[0]{align_type} ne '' );
    $$data_r{seed_align} = $$ref3[0]{alignment}
      if ( $$ref3[0]{alignment} ne '' );
    $$data_r{full_align_id} = $$ref4[0]{align_id}
      if ( $$ref4[0]{align_id} ne '' );
    $$data_r{full_align_name} = $$ref4[0]{align_name}
      if ( $$ref4[0]{align_name} ne '' );
    $$data_r{method_full} = $$ref4[0]{method}
      if ( $$ref4[0]{method} ne '' );
    $$data_r{full_user_class} = $$ref4[0]{user_class}
      if ( $$ref4[0]{user_class} ne '' );
    $$data_r{full_comment} = $$ref4[0]{comment}
      if ( $$ref4[0]{comment} ne '' );
    $$data_r{full_ref_id} = $$ref4[0]{ref_id}
      if ( $$ref4[0]{ref_id} ne '' );
    $$data_r{full_align_assigned_by} = $$ref4[0]{assigned_by}
      if ( $$ref4[0]{assigned_by} ne '' );
    $$data_r{full_align_type} = $$ref4[0]{align_type}
      if ( $$ref4[0]{align_type} ne '' );
    $$data_r{full_align} = $$ref4[0]{alignment}
      if ( $$ref4[0]{alignment} ne '' );

    #data in hmm_role_link
    if ( ref $ref5 eq "ARRAY" ) {
        for my $n (@$ref5) {
            $$data_r{role}{ $$n{hrl_id} }{role_id} = $$n{role_id}
              if ( $$n{hrl_id} ne '' );
            $$data_r{role}{ $$n{hrl_id} }{curated} = $$n{curated}
              if ( $$n{hrl_id} ne '' );
            $$data_r{role}{ $$n{hrl_id} }{owner} = $$n{owner}
              if ( $$n{hrl_id} ne '' );
            $$data_r{role}{ $$n{hrl_id} }{mod_date} = $$n{mod_date}
              if ( $$n{hrl_id} ne '' );
        }
    }
    if ( ref $ref6 eq "ARRAY" ) {
        for my $n (@$ref6) {
            $$data_r{go}{ $$n{id} }{go_id} = $$n{go_term}
              if ( $$n{id} ne '' );
            $$data_r{go}{ $$n{id} }{curated} = $$n{curated}
              if ( $$n{id} ne '' );
            $$data_r{go}{ $$n{id} }{owner} = $$n{owner}
              if ( $$n{id} ne '' );
            $$data_r{go}{ $$n{id} }{mod_date} = $$n{mod_date}
              if ( $$n{id} ne '' );
        }
    }
    chomp $$data_r{hmm};
    chomp $$data_r{hmm_frag};
    chomp $$data_r{hmm_seed};
    chomp $$data_r{seed_align};
    chomp $$data_r{full_align};

    #process seed and full
    if ( $$data_r{hmm_seed} ne '' ) {
        $$data_r{seed}{ori} = $$data_r{hmm_seed};
        &parse_sequence( $$data_r{seed}, 'mul', '', '', '', '', '' );
    }
    if ( $$data_r{full_align} ne '' ) {
        $$data_r{full}{ori} = $$data_r{full_align};
        &parse_sequence( $$data_r{full}, 'fasta', '', '', '', '', '' );
    }
}

sub create_INFO_file {

    #create INFO information using data in DB
    #level 1: special INFO file for daemon etc.
    #level 2: INFO file for internal use
    #level 3: INFO file for external use
    use strict;
    my ( $hash_r, $level ) = @_;
    my $out;
    $out .= "ID  $$hash_r{hmm_name}\n"
      if ( $$hash_r{hmm_name} =~ /\S+/ );
    $out .= "AC  $$hash_r{hmm_acc}\n"
      if ( $$hash_r{hmm_acc} =~ /\S+/ );
    $out .= "DE  $$hash_r{hmm_com_name}\n"
      if ( $$hash_r{hmm_com_name} =~ /\S+/ );
    $out .= "AU  $$hash_r{author}\n"
      if ( $$hash_r{author} =~ /\S+/ );
    $out .= "TC  $$hash_r{trusted_cutoff} $$hash_r{trusted_cutoff2}\n"
      if ( $$hash_r{trusted_cutoff} =~ /\S+/ );
    $out .= "NC  $$hash_r{noise_cutoff} $$hash_r{noise_cutoff2}\n"
      if ( $$hash_r{noise_cutoff} =~ /\S+/ );
    $out .= "AL  $1\n"
      if ( $$hash_r{method_seed} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "GA  $1\n"
      if ( $$hash_r{search_prog} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "IT  $1\n"
      if ( $$hash_r{iso_type} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "RI  $1\n"
      if ( $$hash_r{role_id} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "SO  $1\n"
      if ( $$hash_r{search_options} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "EN  $1\n"
      if ( $$hash_r{expanded_name} =~ /^\s*(\S+.+)\s*$/ );

    #    $out .= "SA  $1\n" if($$hash_r{struct_acc}     =~ /^\s*(\S+.+)\s*$/);
    $out .= "GS  $1\n"
      if ( $$hash_r{gene_sym} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "EC  $1\n"
      if ( $$hash_r{ec_num} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "TR  $1\n"
      if ( $$hash_r{tc_num} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "TP  $1\n"
      if ( $$hash_r{hmm_type} =~ /^\s*(\S+.+)\s*$/ );
    $out .= "RH  $1\n"
      if ( $$hash_r{related_hmm} =~ /^\s*(\S+.+)\s*$/ );
    if ( $$hash_r{hmm_comment} =~ /\S/ ) {
        my $tmp = $$hash_r{hmm_comment};
        $tmp =~ s/\A/CC  /;
        $tmp =~ s/\n/\nCC  /g;
        $out .= $tmp;
        chomp $out;
        $out .= "\n";
    }
    if (   $level =~ /^(1|2)$/
        && $$hash_r{private} =~ /\S/ )
    {
        my $tmp = $$hash_r{private};
        $tmp =~ s/\A/PC  /;
        $tmp =~ s/\n/\nPC  /g;
        $out .= $tmp;
        chomp $out;
        $out .= "\n";
    }
    my $reference = $$hash_r{reference};
    $reference =~ s/\nAL\s+\S.+/\n/g
      if ( $$hash_r{method_seed} =~ /\S/ );
    $out .= $reference
      if ( $reference =~ /\S/ );
    chomp $out;
    $out .= "\n";
    if ( $level == 1 ) {
        $out .= "TIGR  hmm_id $$hash_r{hmm_id}\n"
          if ( $$hash_r{hmm_id} =~ /\S+/ );
        $out .= "TIGR  iso_id $$hash_r{iso_id}\n"
          if ( $$hash_r{iso_id} =~ /\S+/ );
        $out .= "TIGR  is_current $$hash_r{is_current}\n"
          if ( $$hash_r{is_current} =~ /\S+/ );
        $out .= "TIGR  hmm_method_id $$hash_r{hmm_method_id}\n"
          if ( $$hash_r{hmm_method_id} =~ /\S+/ );
        $out .= "TIGR  hmm_frag_method_id $$hash_r{hmm_frag_method_id}\n"
          if ( $$hash_r{hmm_frag_method_id} =~ /\S+/ );
        $out .= "TIGR  num_seed $$hash_r{num_seed}\n"
          if ( $$hash_r{num_seed} =~ /\S+/ );
        $out .= "TIGR  prior $$hash_r{prior}\n"
          if ( $$hash_r{prior} =~ /\S+/ );
        $out .= "TIGR  avg_score $$hash_r{avg_score}\n"
          if ( $$hash_r{avg_score} =~ /\S+/ );
        $out .= "TIGR  std_dev $$hash_r{std_dev}\n"
          if ( $$hash_r{std_dev} =~ /\S+/ );
        $out .= "TIGR  min_score $$hash_r{min_score}\n"
          if ( $$hash_r{min_score} =~ /\S+/ );
        $out .= "TIGR  max_score $$hash_r{max_score}\n"
          if ( $$hash_r{max_score} =~ /\S+/ );
        $out .= "TIGR  hmm_len $$hash_r{hmm_len}\n"
          if ( $$hash_r{hmm_len} =~ /\S+/ );
        $out .= "TIGR  iso_acc $$hash_r{iso_acc}\n"
          if ( $$hash_r{iso_acc} =~ /\S+/ );
        $out .= "TIGR  iso_iso_type $$hash_r{iso_iso_type}\n"
          if ( $$hash_r{iso_iso_type} =~ /\S+/ );
        $out .= "TIGR  iso_method $$hash_r{iso_method}\n"
          if ( $$hash_r{iso_method} =~ /\S+/ );
        $out .= "TIGR  assigned_by $$hash_r{assigned_by}\n"
          if ( $$hash_r{assigned_by} =~ /\S+/ );
        $out .= "TIGR  class_name $$hash_r{class_name}\n"
          if ( $$hash_r{class_name} =~ /\S+/ );
        $out .= "TIGR  db_type $$hash_r{db_type}\n"
          if ( $$hash_r{db_type} =~ /\S+/ );
        $out .= "TIGR  iso_ref_id $$hash_r{iso_ref_id}\n"
          if ( $$hash_r{iso_ref_id} =~ /\S+/ );
        $out .= "TIGR  user_class $$hash_r{user_class}\n"
          if ( $$hash_r{user_class} =~ /\S+/ );
        $out .= "TIGR  seq_type $$hash_r{seq_type}\n"
          if ( $$hash_r{seq_type} =~ /\S+/ );
        $out .= "TIGR  iso_comment $$hash_r{iso_comment}\n"
          if ( $$hash_r{iso_comment} =~ /\S+/ );
        $out .= "TIGR  alignment_iso_id $$hash_r{alignment_iso_id}\n"
          if ( $$hash_r{alignment_iso_id} =~ /\S+/ );
        $out .= "TIGR  seed_align_id $$hash_r{seed_align_id}\n"
          if ( $$hash_r{seed_align_id} =~ /\S+/ );
        $out .= "TIGR  seed_align_name $$hash_r{seed_align_name}\n"
          if ( $$hash_r{seed_align_name} =~ /\S+/ );
        $out .= "TIGR  seed_user_class $$hash_r{seed_user_class}\n"
          if ( $$hash_r{seed_user_class} =~ /\S+/ );
        $out .= "TIGR  seed_comment $$hash_r{seed_comment}\n"
          if ( $$hash_r{seed_comment} =~ /\S+/ );
        $out .= "TIGR  seed_ref_id $$hash_r{seed_ref_id}\n"
          if ( $$hash_r{seed_ref_id} =~ /\S+/ );
        $out .=
          "TIGR  seed_align_assigned_by $$hash_r{seed_align_assigned_by}\n"
          if ( $$hash_r{seed_align_assigned_by} =~ /\S+/ );
        $out .= "TIGR  seed_align_type $$hash_r{seed_align_type}\n"
          if ( $$hash_r{seed_align_type} =~ /\S+/ );
        $out .= "TIGR  full_align_id $$hash_r{full_align_id}\n"
          if ( $$hash_r{full_align_id} =~ /\S+/ );
        $out .= "TIGR  full_align_name $$hash_r{full_align_name}\n"
          if ( $$hash_r{full_align_name} =~ /\S+/ );
        $out .= "TIGR  method_full $$hash_r{method_full}\n"
          if ( $$hash_r{method_full} =~ /\S+/ );
        $out .= "TIGR  full_user_class $$hash_r{full_user_class}\n"
          if ( $$hash_r{full_user_class} =~ /\S+/ );
        $out .= "TIGR  full_comment $$hash_r{full_comment}\n"
          if ( $$hash_r{full_comment} =~ /\S+/ );
        $out .= "TIGR  full_ref_id $$hash_r{full_ref_id}\n"
          if ( $$hash_r{full_ref_id} =~ /\S+/ );
        $out .=
          "TIGR  full_align_assigned_by $$hash_r{full_align_assigned_by}\n"
          if ( $$hash_r{full_align_assigned_by} =~ /\S+/ );
        $out .= "TIGR  full_align_type $$hash_r{full_align_type}\n"
          if ( $$hash_r{full_align_type} =~ /\S+/ );
        $out .= "TIGR  hmm_file $$hash_r{hmm_file}\n"
          if ( $$hash_r{hmm_file} =~ /\S+/ );
        $out .= "TIGR  frag_file $$hash_r{frag_file}\n"
          if ( $$hash_r{frag_file} =~ /\S+/ );
        $out .= "TIGR  build_file $$hash_r{build_file}\n"
          if ( $$hash_r{build_file} =~ /\S+/ );
        $out .= "TIGR  seed_file $$hash_r{seed_file}\n"
          if ( $$hash_r{seed_file} =~ /\S+/ );
        $out .= "TIGR  full_file $$hash_r{full_file}\n"
          if ( $$hash_r{full_file} =~ /\S+/ );
        $out .= "TIGR  nomen_check $$hash_r{nomen_check}\n"
          if ( $$hash_r{nomen_check} =~ /\S+/ );

        #       $out .= "TIGR   $$hash_r{}\n" if($$hash_r{} =~ /\S+/);
    }
    if ( $level == 3 ) {
        if ( $out =~ /\nRH\s+(\S.+)\n/ ) {
            my $ori_rh = $1;
            $ori_rh =~ s/\Q$$hash_r{hmm_acc}\E//;
            $ori_rh =~ s/^\s+|\s+$//g;
            if ( $ori_rh =~ /\S/ ) {
                $out =~ s/\nRH\s+(\S.+)\n/\nRH  $ori_rh\n/;
            }
            else {
                $out =~ s/\nRH\s+(\S.+)\n/\n/;
            }
        }
        $out =~ s/\nTP\s+TOGA.*/\nTP  TIGRFAMs/;
        $out =~ s/\nA2.*/\n/g;
        $out =~ s/\nBD.*/\n/g;
        $out =~ s/\nGA.*/\n/g;
        $out =~ s/\nSE.*/\n/g;
        $out =~ s/\n[A-Z0-9]{2}\s*$/\n/mg;
    }
    $out =~ s/\n\s*\n/\n/g;
    return $out;
}

sub zero_launch {

    #other system calls also triggers SIG{CHLD}, but no $pid is returned here.
    # so, it can be used to distinguish all system calls from launch_child_proc.
    #any global variable used in this subroutine must be prefixed with main::
    use strict;
    my $child_id = waitpid( -1, &WNOHANG )
      ;    # &WNOHANG is defined by 'use POSIX ":sys_wait_h";'
    if ( $main::child{$child_id} == 1 ) {
        delete $main::child{$child_id};
        --$main::num_child;
        print
"${\(scalar localtime)}: child process $child_id exits and current number of child process is $main::num_child\n";

    #	print LOG "${\(scalar localtime)}: child process $child_id is finished\n";
    }
}

sub launch_child_proc {
    use strict;
    use Carp;
    my ( $sub_r, $para_r, $max_child ) = @_;
    my ( $sub_log, $sub_status ) = ( '', 'error' );
  FORK:
    my $count = 0;
    if ( $main::num_child < $max_child ) {
        ++$count;
        my $pid;
        if ( $pid = fork ) {

 #fork returns child pid to parent, and this section runs in the PARENT process.
            ++$main::num_child;
            $main::child{$pid} = 1;
            $sub_log .=
"successful fork, pid: $pid, total \# of child process is $main::num_child\n";
            $sub_status = 'success';

#	    print LOG "${\(scalar localtime)}: child process $pid is launched with parameters: @$para_r\n";
        }
        elsif ( defined $pid
            && $pid ne '' )
        {

#defined, not string-equal to the null string, CAN be numeric zero.
#fork, if successful, returns 0 to child, and this section runs in the CHILD process.
            &$sub_r(@$para_r);
            if (   $$para_r[1] =~ /^(PF|TIGR)\d+/
                && $$para_r[2] =~ /ogc_mini_db$/ )
            {

#		print  LOG "${\(scalar localtime)}: returned from sub create_mini_db_info in child process for $$para_r[1]\n" if($LOG);
#		print  "${\(scalar localtime)}: returned from sub create_mini_db_info in child process for $$para_r[1]\n";
            }
            print "${\(scalar localtime)}: return from child process\n";

            #	    exit 0;
            die;

            #	    croak "returned from child process\n";
        }
        elsif ($count < 6
            && $! =~ /No more process/ )
        {
            print "fork failed, retry after 10 seconds\n";
            sleep 10;
            redo FORK;
        }
        else {
            $sub_log .= "can't fork: $!\n";
            $sub_status = 'error';
        }
    }
    return ( $sub_status, $sub_log );
}

sub update_iso_link_coor {
    use strict;
    my ( $db_proc, $db, $revise, $acc, $LOG ) = @_;
    my $ref = &do_query( $db_proc,
"select iso_id from $db..hmm2 where hmm_acc = \"$acc\" and is_current = 1"
    );
    return
      if ( @$ref == 0 );
    my $iso_id = $$ref[0]{iso_id};
    print $LOG
      "acc: $acc, iso_id: $iso_id------------------------------------------\n"
      if ($LOG);
    print
      "acc: $acc, iso_id: $iso_id------------------------------------------\n"
      if ($LOG);
    my $ref = &do_query( $db_proc,
"select alignment from $db..alignment where iso_id = $iso_id and align_type = \"all\""
    );
    my %seq;
    $seq{ori} = $$ref[0]{alignment};
    &parse_sequence( \%seq, 'fasta', '', '', '', $LOG, '' );

#get full length protein sequences from protein table for a list of protein EGAD ids
    my %prot;
    for my $n ( 0 .. $seq{number} - 1 ) {
        next
          if ( $seq{$n}{prot_id} !~ /^\d+$/ );
        my $ref = &do_query( $db_proc,
"select prot_id, prot_seq from $db..protein where prot_id = $seq{$n}{prot_id}"
        );
        $prot{ $$ref[0]{prot_id} } = $$ref[0]{prot_seq}
          if ( @$ref == 1 );
    }

    #get old info from iso_link
    my $ref = &do_query( $db_proc,
        "select prot_id, lend, rend from $db..iso_link where iso_id = $iso_id"
    );
    my %old_iso_link;
    for (@$ref) {
        $old_iso_link{ $$_{prot_id} }{lend} = $$_{lend};
        $old_iso_link{ $$_{prot_id} }{rend} = $$_{rend};
    }

    #match two sequence
    my $count = 0;
    for my $n ( 0 .. $seq{number} - 1 ) {
        next
          if ( $seq{$n}{prot_id} !~ /^\d+$/ );
        ++$count;
        my $prot_id = $seq{$n}{prot_id};
        my ( $success, $a_left, $a_right, $b_left, $b_right, $lend, $rend ) =
          &compare_two_seq( uc $seq{$n}{seq}, uc $prot{$prot_id} )
          if ( $seq{$n}{seq} ne ''
            && $prot{$prot_id} ne '' );
        print $LOG
"+++Warning: not a substring: aligned sequence is truncated at ($a_left, $a_right) at ends while its length is ${\(length $seq{$n}{seq})}\n$n, $prot_id\n$seq{$n}{seq}\n$prot{$prot_id}\n\n"
          if ( $a_left != 1
            or $a_right < length( $seq{$n}{seq} ) );
        if ($success) {
            $seq{$n}{lend} = $lend;
            $seq{$n}{rend} = $rend;
        }
        if (   $lend != $old_iso_link{$prot_id}{lend}
            or $rend != $old_iso_link{$prot_id}{rend} )
        {
            print $LOG
"+++Warning: mismatch between old and calc coor for seq\# $n in iso_link: old($old_iso_link{$prot_id}{lend}, $old_iso_link{$prot_id}{rend}), calc($lend, $rend)\n$n, $prot_id, length of them are ${\(length $seq{$n}{seq})}, ${\(length $prot{$prot_id})}\n$seq{$n}{seq}\n$prot{$prot_id}\n\n"
              if ( $old_iso_link{$prot_id}{lend} ne ''
                or $old_iso_link{$prot_id}{rend} ne '' );
        }
    }
    print $LOG
      "CC $count, ${\(scalar(keys %prot))}, ${\(scalar(keys %old_iso_link))}\n";
    my @new_iso_link;
    &construct_iso_link( \%seq, $iso_id, \@new_iso_link );
    &load_iso_link( 1, $iso_id, \@new_iso_link, $db, $db_proc, $LOG, '', '' )
      if ($revise);
}

sub compare_two_seq {

    #seq $b (probably full length protein sequence) should be longer than $a
    use strict;
    my ( $a, $b ) = @_;
    my ( $success, $a_left, $a_right, $b_left, $b_right ) =
      &match_two_string( $a, $b );
    $success = 0
      if (
        $success
        && (   ( $a_right - $a_left + 1 ) / length($a) < 0.80
            && ( $b_right - $b_left + 1 ) / length($b) < 0.80 )
      );
    my ( $lend, $rend ) = ( 0, 0 );
    if (   $success
        && $b_left - $a_left + 1 > 0
        && $b_left - $a_left + length($a) <= length($b) )
    {
        $lend = $b_left - $a_left + 1;
        $rend = $b_left - $a_left + length($a);
    }
    return ( $success, $a_left, $a_right, $b_left, $b_right, $lend, $rend );
}

sub get_coor_from_nraa {
    use strict;
    my ( $db_file, $seq_r ) = @_;
    for my $n ( 0 .. $$seq_r{number} - 1 ) {
        my ($h) = split( /\s+/, $$seq_r{$n}{header} );
        if ( $h =~ /^GB/ ) {
            $h =~ s/GB\|(.*)\|(.*)\|(.*)/GB\|$1\|$2/gi;
        }
        my $result = '';
        open( FH, "yank_panda -a \"$h\" | " );

        #	open(FH, "yank -t -a \"$h\" | ");
        while (<FH>) {
            $result .= $_;
        }
        close FH;
        if (   $result ne ''
            && $result =~ /^>/ )
        {
            my %tmp_seq;
            $tmp_seq{ori} = $result;
            &parse_sequence( \%tmp_seq, 'fasta', '', '', '', '', '' );
            my $s1 = $$seq_r{$n}{seq};
            my $s2 = $tmp_seq{0}{seq};
            my ( $success, $s1_left, $s1_right, $s2_left, $s2_right ) =
              &match_two_string( $s1, $s2 );
            $success = 0
              if (
                $success
                && (   ( $s1_right - $s1_left + 1 ) / length($s1) < 0.80
                    && ( $s2_right - $s2_left + 1 ) / length($s2) < 0.80 )
              );
            if (   $success
                && $s2_left - $s1_left + 1 > 0
                && $s2_left - $s1_left + length($s1) <= length($s2) )
            {
                $$seq_r{$n}{lend} = $s2_left - $s1_left + 1;
                $$seq_r{$n}{rend} = $s2_left - $s1_left + length($s1);
            }
        }
    }
}

sub tmp_dir {

#create/delete directories under $base_dir depending on $mode. For $mode="create", if $base_dir="/use/local/project/tmp_dir"
#and $new_dir="myDir", the directory /usr/local/project/tmp_dir/myDir1/ will be created. If myDir1 already exists, myDir2 will
#be created.  If myDir1 is deleted after myDir2 is created, another myDir1 will be created.
#If $mode="delete" and $new_dir="myDir1", the $new_dir will be removed completely if permission is set properly.  This means,
#the calling program should call it twice for creating and deleting it.
#If $clean_time is set, any directory named $new_dir with a date older than that time from the current time will be deleted.
    use strict;

    # clean not implemented yet.
    my ( $mode, $base_dir, $new_dir, $created_dir_r, $clean_time, $log_r ) = @_;
    my $curr_dir = cwd();

    #check that it exists and is a directory
    unless ( -d $base_dir ) {
        $$log_r .=
          "the base directory $base_dir does not exist or is not a direcotry\n"
          if ( $log_r ne '' );
        chdir $curr_dir;
        return 1;
    }
    if ( $mode eq 'create' ) {

   #check that it has enough permission and try to modify permission if possible
        unless ( -r $base_dir
            && -w $base_dir
            && -x $base_dir )
        {
            chmod 0755, $base_dir;
            if (   !-r $base_dir
                or !-w $base_dir
                or !-x $base_dir )
            {
                $$log_r .=
"the program does not have enough permission in base directory $base_dir\n"
                  if ( $log_r ne '' );
                chdir $curr_dir;
                return 1;
            }
        }
        if ( chdir $base_dir ) {
            my @list = <$new_dir*>;
            my %tmp;
            for my $name (@list) {
                $tmp{$1} = 1
                  if ( $name =~ /^$new_dir(\d+)$/ );
            }
            my $num = 0;
            while ( $tmp{$num} == 1 ) {
                ++$num;
            }
            system "rm -rf $new_dir$num";
            mkdir "$new_dir$num", 0755;
            $$created_dir_r = "$new_dir$num";
        }
    }
    elsif ( $mode eq 'delete' ) {

   #check that it has enough permission and try to modify permission if possible
        unless ( -w $base_dir
            && -x $base_dir )
        {
            chmod 0755, $base_dir;
            if (   !-w $base_dir
                or !-x $base_dir )
            {
                $$log_r .=
"the program does not have enough permission in base directory $base_dir\n"
                  if ( $log_r ne '' );
                chdir $curr_dir;
                return 1;
            }
        }
        chmod 0777, "$base_dir/$new_dir";
        system "chmod 777 $base_dir/$new_dir";
        system "rm -rf $base_dir/$new_dir";
        if ( -d "$base_dir/$new_dir" ) {
            $$log_r = "+++Warning: can not delete $base_dir/$new_dir\n"
              if ( $log_r ne '' );
            chdir $curr_dir;
            return 1;
        }
        else {
            $$log_r = "$base_dir/$new_dir is deleted successfully\n"
              if ( $log_r ne '' );
        }
    }
    elsif ( $mode eq 'clean' ) {
    }
    else {
        $$log_r .= "unknown mode\n"
          if ( $log_r ne '' );
        chdir $curr_dir;
        return 1;
    }
    chdir $curr_dir;
    return 0;
}

sub tmp_file {

# create and return a file name for temporary usage.  The program that calls it should delete the temporary file.
# the file is in the directory passed to the subroutine. If the file already exists in the directory, the file
# name will be filename1, filename2, etc.
    use strict;
    my ( $dir, $file ) = @_;
    my @list = <$dir/$file*>;
    for (@list) {
        $_ = $1
          if ( $_ =~ /([^\/]+)$/ );
    }
    my %tmp;
    for my $name (@list) {
        $tmp{$1} = 1
          if ( $name =~ /^\Q$file\E(\d+)$/ );
    }
    my $num = 0;
    while ( $tmp{$num} == 1 ) {
        ++$num;
    }
    system("touch $dir/$file$num");
    return "$file$num";
}

sub load_tree {
    use strict;
    my (
        $db,           $db_proc,  $tree_file, $tree,
        $iso_id,       $align_id, $tree_name, $tree_method,
        $tree_comment, $tree_assigned_by
    ) = @_;
    if ( $tree eq ''
        && open( FH, "$tree_file" ) )
    {
        while (<FH>) {
            $tree .= $_;
        }
        close FH;
    }
    return
      if ( $tree eq '' );

#    print "insert $db..tree (tree_name, tree, method, assigned_by, comment) values (\"$tree_name\", \"$tree\", \"$tree_method\", \"$tree_assigned_by\", \"$tree_comment\")\n";
    my $ref = &do_query( $db_proc,
"insert $db..tree (tree_name, tree, method, assigned_by, comment) values (\"$tree_name\", \"$tree\", \"$tree_method\", \"$tree_assigned_by\", \"$tree_comment\")"
    );
    my $tree_id = $$ref[0]{'COL(1)'};
    &do_query(
        $db_proc,
"insert $db..tree_link (iso_id, tree_id, align_id) values ($iso_id, $tree_id, $align_id)",
        1
      )
      if ( $iso_id ne ''
        && $align_id ne '' );
}

sub calc_coor_for_seq {
    use strict;
    my (
        $db,      $db_proc, $seq_r, $blastp, $flag_r,
        $fatal_r, $header,  $LOG1,  $LOG2
    ) = @_;
    my ( $prot_id_list, %prot );
    for my $k ( 0 .. $$seq_r{number} - 1 ) {
        $prot_id_list .= "$$seq_r{$k}{prot_id}, "
          if ( $$seq_r{$k}{prot_id} ne '' );
    }
    $prot_id_list =~ s/, $//;
    my $ref = &do_query( $db_proc,
"select prot_id, prot_seq from $db..protein where prot_id in ($prot_id_list) order by prot_id"
    ) if ( $prot_id_list ne '' );
    for (@$ref) {
        $prot{ $$_{prot_id} }{prot_seq} = $$_{prot_seq}
          if ( $$_{prot_seq} ne ''
            && $$_{prot_id} ne '' );
    }
    &check_lend_rend_in_MA(
        $seq_r,  '',       \%prot,  $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    );
}

sub delete_proj {
    use strict;

#the order of deletion might be important, delete status, alignemnt, isology, iso_link, iso_project, project_iso_link
#note that hmm2 table is not involved yet
    my ( $db, $db_proc, $proj_id, $ogc_dir ) = @_;

    #set up some parameters
    my ( %status_vocab_ids, $iso_id_list );
    my $ref = &do_query( $db_proc,
        "select id, table_name from $db..status_vocabulary" );
    for (@$ref) {
        $status_vocab_ids{ $$_{table_name} } .= "$$_{id}, ";
    }
    for my $k ( keys %status_vocab_ids ) {
        $status_vocab_ids{$k} =~ s/, $//;
    }
    $ref = &do_query( $db_proc,
        "select iso_id from $db..project_iso_link where project_id = $proj_id"
    );
    for (@$ref) {
        $iso_id_list .= "$$_{iso_id}, ";
    }
    $iso_id_list =~ s/, $//;

    #delete data from DB
    &do_query(
        $db_proc,
"delete $db..status from $db..status s, $db..alignment a where link_id = iso_id and iso_id in ($iso_id_list) and status_vocab_id in ($status_vocab_ids{alignment})",
        1
    );
    &do_query(
        $db_proc,
"delete $db..status where link_id in ($iso_id_list) and status_vocab_id in ($status_vocab_ids{isology})",
        1
    );
    &do_query(
        $db_proc,
"delete $db..status where link_id in (select id from $db..project_iso_link where project_id = $proj_id) and status_vocab_id in ($status_vocab_ids{project_iso_link})",
        1
    );
    &do_query(
        $db_proc,
"delete $db..status where link_id = $proj_id and status_vocab_id in ($status_vocab_ids{iso_project})",
        1
    );
    &do_query( $db_proc, "delete $db..status where project_id = $proj_id", 1 );
    &do_query( $db_proc, "delete $db..alignment where iso_id in ($iso_id_list)",
        1 );
    &do_query( $db_proc, "delete $db..isology where iso_id in ($iso_id_list)",
        1 );
    &do_query( $db_proc, "delete $db..iso_link where iso_id in ($iso_id_list)",
        1 );
    &do_query( $db_proc,
        "delete $db..project_iso_link where project_id = $proj_id", 1 );
    &do_query( $db_proc, "delete $db..iso_project where project_id = $proj_id",
        1 );
    system "rm -rf $ogc_dir/ogc_results/$proj_id";
}

sub check_info_complete {
    use strict;
    my ( $file, $hmm_name_r ) = @_;
    open( FH, "$file" );
    my ( $data, %data );
    while (<FH>) {
        $data .= $_;
    }
    close FH;
    &get_INFO_from_file( \%data, $file );

    #check not null
    for (
        'hmm_name', 'hmm_com_name',  'iso_type',
        'author',   'expanded_name', 'method_seed'
      )
    {
        return 0
          if ( !defined $data{$_}
            or $data{$_} eq '' );
    }

    #check not null and digital
    for ( 'noise_cutoff', 'trusted_cutoff' ) {
        return 0
          if ( !defined $data{$_}
            or $data{$_} eq ''
            or $data{$_} =~ /[^0-9.+\-]/ );
    }

    #check if IT(iso_type)  is equivalog,  RI is required
    return 0
      if (
        $data{iso_type} eq 'equivalog'
        && ( !defined $data{role_id}
            or $data{role_id} eq '' )
      );

    #check other condition
    $$hmm_name_r = $data{hmm_name};
    return 1;
}

sub load_db {
    my ( $db, $db_proc, $inup, $table, $return_r, $where, %field ) = @_;
    my $noreturn = 0;
    $noreturn = 1
      if ( $return_r eq '' );
}

sub read_hmmbuild_para {
    use strict;
    my ( $para_r, $file ) = @_;
    open( FH, "$file" );
    <FH>;
    while (<FH>) {
        chomp $_;
        my ( $name, $tmp ) = split( /<>/, $_, 2 );
        (
            $$para_r{$name}{com_name}, $$para_r{$name}{default},
            $$para_r{$name}{command},  $$para_r{$name}{comment},
            $$para_r{$name}{author},   $$para_r{$name}{alias}
        ) = split( /<>/, $tmp );
    }
    close FH;
}

sub get_projID_from_dir {
    my ($proj_r) = @_;
    my @proj_id = <*>;
    for my $proj_id (@proj_id) {
        if ( $proj_id =~ /^\d+$/
            && -d $proj_id )
        {
            chdir "$proj_id";
            &get_isoID_from_dir( \%{ $$proj_r{$proj_id} } );
            chdir "..";
        }
    }
}

sub process_check_hmm_data {
    use strict;
    my ( $db, $db_proc, $class, $data_r, $name, $fatal_r, $LOG1, $LOG2 ) = @_;
    my ( @HMM, @FRAG );
    @HMM = split( /\n/, $$data_r{hmm}{value} )
      if ( $$data_r{hmm}{value} );
    @FRAG = split( /\n/, $$data_r{frag}{value} )
      if ( $$data_r{frag}{value} );
    print $LOG1 "***Error: unknown meaning of file name $name\n"
      if ( $class !~ /^(TOGA|PFAM|OGC)$/ );
    print $LOG2 "***Error: unknown meaning of file name $name\n"
      if ( $class !~ /^(TOGA|PFAM|OGC)$/ );

# check hmm_acc.  hmm_acc might be found in INFO stored in $data{hmm_acc}{value}, HMM, FRAG, or as the file name for PFAM
    my ( $HMM_acc, $FRAG_acc, $file_acc );
    for (@HMM) {
        if ( $_ =~ /^ACC\s+(\S*)\s*/ ) {
            $HMM_acc = $1;
            last;
        }
    }
    for (@FRAG) {
        if ( $_ =~ /^ACC\s+(\S*)\s*/ ) {
            $FRAG_acc = $1;
            last;
        }
    }
    $file_acc = $name
      if ( $class eq 'PFAM' );
    for my $s1 ( $$data_r{hmm_acc}{value}, $HMM_acc, $FRAG_acc, $file_acc ) {
        for my $s2 ( $$data_r{hmm_acc}{value}, $HMM_acc, $FRAG_acc, $file_acc )
        {
            if (   $s1 =~ /\S/
                && $s2 =~ /\S/
                && $s1 ne $s2 )
            {
                print $LOG1
                  "***Error: mismatch of hmm_acc among data, $s1, $s2\n";
                print $LOG2
                  "***Error: mismatch of hmm_acc among data, $s1, $s2\n";
                $$fatal_r .= "mismatch of hmm_acc among data, $s1, $s2\n"
                  if ($fatal_r);
            }
        }
    }
    for ( $HMM_acc, $FRAG_acc, $file_acc ) {
        $$data_r{hmm_acc}{value} = $_
          if ( $_ ne '' );
    }
    $$data_r{hmm_acc}{value} = &get_avail_hmm_acc( $db, $db_proc, 1 )
      if ( $$data_r{hmm_acc}{value} eq '' );
    print "$name, hmm_acc: $$data_r{hmm_acc}{value}\n";

    #update HMM FRAG if necessary.
    # check hmm_type that might be found in INFO file
    print $LOG1 "***Error: unknown HMM type: $$data_r{hmm_type}{value}\n"
      if ( $$data_r{hmm_type}{value} !~ /(TOGA|PFAM)/
        && $$data_r{hmm_type}{value} =~ /\S/ );
    print $LOG2 "***Error: unknown HMM type: $$data_r{hmm_type}{value}\n"
      if ( $$data_r{hmm_type}{value} !~ /(TOGA|PFAM)/
        && $$data_r{hmm_type}{value} =~ /\S/ );
    print $LOG1
      "***Error: wrong HMM type for $class: $$data_r{hmm_type}{value}\n"
      if (
        (
               $class eq 'PFAM'
            && index( $$data_r{hmm_type}{value}, 'PFAM' ) < 0
            && $$data_r{hmm_type}{value} =~ /\S/
        )
        or ( $class eq 'TOGA'
            && index( $$data_r{hmm_type}{value}, 'TOGA' ) < 0 )
        && $$data_r{hmm_type}{value} =~ /\S/
      );
    print $LOG2
      "***Error: wrong HMM type for $class: $$data_r{hmm_type}{value}\n"
      if (
        (
               $class eq 'PFAM'
            && index( $$data_r{hmm_type}{value}, 'PFAM' ) < 0
            && $$data_r{hmm_type}{value} =~ /\S/
        )
        or ( $class eq 'TOGA'
            && index( $$data_r{hmm_type}{value}, 'TOGA' ) < 0 )
        && $$data_r{hmm_type}{value} =~ /\S/
      );

#check hmm_name that might be found in field ID in INFO stored in $data{hmm_name}{value}, field NAME in PFAM HMM, field NAME in PFAM FRAG, and file name for TOGA
    my ( $HMM_name, $FRAG_name, $file_name );
    for (@HMM) {
        if (   $class eq 'PFAM'
            && $_ =~ /^NAME\s+(\S*)\s*/
            && $_ !~ /^NAME\s+PF\d{5,}\s*$/ )
        {
            $HMM_name = $1;
            last;
        }
        elsif ($class eq 'TOGA'
            && $_ =~ /^NAME\s+(\S*)\s*/
            && $_ !~ /^NAME\s+TIGR\d{5,}\s*$/ )
        {
            $HMM_name = $1;
            last;
        }
    }
    for (@FRAG) {
        if (   $class eq 'PFAM'
            && $_ =~ /^NAME\s+(\S*)\s*/
            && $_ !~ /^NAME\s+PF\d{5,}\s*$/ )
        {
            $FRAG_name = $1;
            last;
        }
        if (   $class eq 'TOGA'
            && $_ =~ /^NAME\s+(\S*)\s*/
            && $_ !~ /^NAME\s+TIGR\d{5,}\s*$/ )
        {
            $FRAG_name = $1;
            last;
        }
    }
    $file_name = $name
      if ( $class eq 'TOGA' );
    for my $s1 ( $$data_r{hmm_name}{value}, $HMM_name, $FRAG_name, $file_name )
    {
        for my $s2 ( $$data_r{hmm_name}{value}, $HMM_name, $FRAG_name,
            $file_name )
        {
            if (   $s1 =~ /\S/
                && $s2 =~ /\S/
                && $s1 ne $s2 )
            {
                print $LOG1
"+++Warning: mismatch of hmm_name among data, $s1, $s2, using file name or ID in INFO file if they match each other\n";
                print $LOG2
"+++Warning: mismatch of hmm_name among data, $s1, $s2, using file name or ID in INFO file if they match each other\n";
            }
        }
    }
    if (   $$data_r{hmm_name}{value} =~ /\S/
        && $file_name =~ /\S/
        && $$data_r{hmm_name}{value} ne $file_name )
    {
        print $LOG1
"***Error: mismatch of hmm_name between file name and ID in INFO file, $$data_r{hmm_name}{value}, $file_name\n";
        print $LOG2
"***Error: mismatch of hmm_name between file name and ID in INFO file, $$data_r{hmm_name}{value}, $file_name\n";
        $$fatal_r .=
"mismatch of hmm_name between file name and ID in INFO file, $$data_r{hmm_name}{value}, $file_name\n"
          if ($fatal_r);
    }
    $$data_r{hmm_name}{value} = $file_name
      if ( $file_name =~ /\S/ );
    print "$name, hmm_name: $$data_r{hmm_name}{value}\n";

#update HMM FRAG if necessary.
#check hmm_com_name that might be found in field DE in INFO, field DESC in HMM and FRAG.
#check hmm_len that mighr be found in BUILD stored in $data{hmm_len}{value}, HMM, and FRAG.
    my ( $HMM_LENG, $HMM_length, $FRAG_LENG, $FRAG_length );
    for (@HMM) {
        if ( $_ =~ /^LENG\s+(\d+)\s*$/ ) {
            $HMM_LENG = $1;
            last;
        }
    }
    for (@FRAG) {
        if ( $_ =~ /^LENG\s+(\d+)\s*$/ ) {
            $FRAG_LENG = $1;
            last;
        }
    }
    my $ind = @HMM;
    while ( $ind > -1 ) {
        --$ind;
        if ( $HMM[$ind] =~ /^\s+(\d+)\s+/ ) {
            $HMM_length = $1;
            last;
        }
    }
    $ind = @FRAG;
    while ( $ind > -1 ) {
        --$ind;
        if ( $FRAG[$ind] =~ /^\s+(\d+)\s+/ ) {
            $FRAG_length = $1;
            last;
        }
    }
    for my $s1 ( $$data_r{hmm_len}{value},
        $HMM_LENG, $HMM_length, $FRAG_LENG, $FRAG_length )
    {
        for my $s2 ( $$data_r{hmm_len}{value},
            $HMM_LENG, $HMM_length, $FRAG_LENG, $FRAG_length )
        {
            if (   $s1 =~ /\S/
                && $s2 =~ /\S/
                && $s1 != $s2 )
            {
                print $LOG1
"+++Warning: wrong hmm_len among data, $s1, $s2, using actual hmm_len\n";
                print $LOG2
"+++Warning: wrong hmm_len among data, $s1, $s2, using actual hmm_len\n";
            }
        }
    }
    if (   $HMM_length =~ /\d/
        && $FRAG_length =~ /\d/
        && $HMM_length != $FRAG_length )
    {
        print $LOG1 "***Error: HMM and FRAG models have different length\n";
        print $LOG2 "***Error: HMM and FRAG models have different length\n";
        $$fatal_r .= "HMM and FRAG models have different length\n"
          if ($fatal_r);
    }
    $$data_r{hmm_len}{value} = $HMM_length;

# update HMM, FRAG
# check number of SEED that contained in SQ field in INFO for PFAM, BUILD file stored in $data_r{num_seed}{value}, and calculated from actual SEED file
    my ( $calc_num_seed, $info_num_seed );
    print $LOG1
"***Error: mismatch of number of SEED sequences, use calculated sequence number\n"
      if (
        (
               $$data_r{seed}{number} != $$data_r{num_seed}{value}
            && $$data_r{seed}{number}    ne ''
            && $$data_r{num_seed}{value} ne ''
        )
        || (   $$data_r{seed}{number} != $$data_r{SQ}
            && $$data_r{seed}{number} ne ''
            && $$data_r{SQ} ne '' )
        || (   $$data_r{num_seed}{value} != $$data_r{SQ}
            && $$data_r{num_seed}{value} ne ''
            && $$data_r{SQ} ne '' )
      );
    print $LOG2
"***Error: mismatch of number of SEED sequences, use calculated sequence number\n"
      if (
        (
               $$data_r{seed}{number} != $$data_r{num_seed}{value}
            && $$data_r{seed}{number}    ne ''
            && $$data_r{num_seed}{value} ne ''
        )
        || (   $$data_r{seed}{number} != $$data_r{SQ}
            && $$data_r{seed}{number} ne ''
            && $$data_r{SQ} ne '' )
        || (   $$data_r{num_seed}{value} != $$data_r{SQ}
            && $$data_r{num_seed}{value} ne ''
            && $$data_r{SQ} ne '' )
      );
    $$data_r{num_seed}{value} = $$data_r{seed}{number}
      if ( $$data_r{seed}{number} > 0 );
    $$data_r{num_seed}{value} = $$data_r{SQ}
      if ( $$data_r{num_seed}{value} eq ''
        && $$data_r{SQ} > 0 );

    #check other
}

sub init_align {
    use strict;
    my ( $iso_id, $para_r ) = @_;
    my $cmd =
"/usr/local/bin/clustalw $$para_r{option} -infile=$$para_r{input} -outfile=$$para_r{output} > $$para_r{output}.log";
    system($cmd);
}

sub trim_seq {
    use strict;
    my ( $seq_r, $left_r, $right_r, $percent_cutoff, $seq_name, $log_header,
        $log_r )
      = @_;
    $log_header .= "--sub trim_seq";
    for my $n ( 1 .. $$seq_r{number} - 1 ) {
        $$log_r .=
"***Error: $log_header--different lengths of aligned sequences $seq_name: ${\length($$seq_r{$n}{seq_gap})}, ${\length($$seq_r{0}{seq_gap})}\n"
          if (
            length( $$seq_r{$n}{seq_gap} ) != length( $$seq_r{0}{seq_gap} ) );
    }
    for my $m ( 0 .. length( $$seq_r{0}{seq_gap} ) - 1 ) {
        my $count = 0;
        for my $n ( 0 .. $$seq_r{number} - 1 ) {
            ++$count
              if ( substr( $$seq_r{$n}{'seq_gap'}, $m, 1 ) =~ /^[a-zA-Z]$/ );
        }
        if ( $count >= $$seq_r{number} * $percent_cutoff / 100 ) {
            $$left_r = $m;
            last;
        }
    }
    for ( my $m = length( $$seq_r{0}{seq_gap} ) - 1 ; $m >= 0 ; --$m ) {
        my $count = 0;
        for my $n ( 0 .. $$seq_r{number} - 1 ) {
            ++$count
              if ( substr( $$seq_r{$n}{'seq_gap'}, $m, 1 ) =~ /^[a-zA-Z]$/ );
        }
        if ( $count >= $$seq_r{number} * $percent_cutoff / 100 ) {
            $$right_r = $m;
            last;
        }
    }
}

sub trim_align {
    use strict;
    my ( $iso_id, $ver, $para_r ) = @_;
    my @file = <*>;
    my %seq;
    &read_sequence( \%seq, "${iso_id}_${ver}_clus.msf", 'msf', '', '' );
    open( F, ">${iso_id}_${ver}_clus_trim.fa" );
    for my $n ( 1 .. $seq{number} - 1 ) {
        print
"***Error: different lengths of aligned sequences ${\length($seq{$n}{seq_gap})}, ${\length($seq{0}{seq_gap})}\n"
          if ( length( $seq{$n}{seq_gap} ) != length( $seq{0}{seq_gap} ) );
    }
    my ( $left, $right );
    for my $m ( 0 .. length( $seq{0}{seq_gap} ) - 1 ) {
        my $count = 0;
        for my $n ( 0 .. $seq{number} - 1 ) {
            ++$count
              if ( substr( $seq{$n}{'seq_gap'}, $m, 1 ) =~ /^[a-zA-Z]$/ );
        }
        if ( $count >= $seq{number} / 2 - 0.51 ) {
            $left = $m;
            last;
        }
    }
    for ( my $m = length( $seq{0}{seq_gap} ) - 1 ; $m >= 0 ; --$m ) {
        my $count = 0;
        for my $n ( 0 .. $seq{number} - 1 ) {
            ++$count
              if ( substr( $seq{$n}{'seq_gap'}, $m, 1 ) =~ /^[a-zA-Z]$/ );
        }
        if ( $count >= $seq{number} / 2 - 0.51 ) {
            $right = $m;
            last;
        }
    }
    if ( $left >= $right ) {
        $left  = 0;
        $right = length( $seq{0}{seq_gap} ) - 1;
    }
    for my $n ( 0 .. $seq{number} - 1 ) {
        $seq{$n}{seq_gap} =
          substr( $seq{$n}{seq_gap}, $left, $right - $left + 1 );
    }
    print F &format_sequence( \%seq, 'fasta', 'all', 60, 1, '', '', '' );
    close F;
    my $cmd =
"clustalw -output=gcg -infile=${iso_id}_${ver}_clus_trim.fa -outfile=${iso_id}_${ver}_clus_trim.msf > ${iso_id}_${ver}_clus_trim.msf.log";
    system($cmd);

    #    unlink "${iso_id}_${ver}_clus_trim.fa";
}

sub hmmbuild {
    use strict;
    my ( $iso_id, $ver, $type, $option ) = @_;
    system(
"hmmbuild $option ${iso_id}_${ver}_$type.HMM ${iso_id}_${ver}.msf > ${iso_id}_${ver}_$type.build"
    );
}

sub calibrate {
    use strict;
    my ( $iso_id, $para_r ) = @_;
    my @file = <*>;
    my $ver = &get_last_version( \@file, "${iso_id}_", "_2p.HMM", '' );
    my $cmd =
"hmmcalibrate --num 5000 --histfile ${iso_id}_${ver}_2p.hist ${iso_id}_${ver}_2p.HMM > ${iso_id}_${ver}_2p.cali";
    system($cmd);
    $cmd =
"hmmcalibrate --num 5000 --histfile ${iso_id}_${ver}_2pf.hist ${iso_id}_${ver}_2pf.HMM > ${iso_id}_${ver}_2pf.cali";
    system($cmd);
    $cmd =
"hmmcalibrate --num 5000 --histfile ${iso_id}_${ver}_2b.hist ${iso_id}_${ver}_2b.HMM > ${iso_id}_${ver}_2b.cali";
    system($cmd);
    $cmd =
"hmmcalibrate --num 5000 --histfile ${iso_id}_${ver}_2bw.hist ${iso_id}_${ver}_2bw.HMM > ${iso_id}_${ver}_2bw.cali";
    system($cmd);
}

sub mini_db {
    use strict;
    my ( $proj_name, $iso_id, $para_r ) = @_;
    my @file = <*>;
    my $ver = &get_last_version( \@file, "${iso_id}_", ".fa", '' );
    open( PARA, "../para_$proj_name" );
    my %para = &get_proj_para( \*PARA );
    close PARA;
    for ( @{ $para{mini_db}{$iso_id} } ) {
        return
          if ( $_ == $ver );
    }
    my $cmd =
"$ENV{HMM_SCRIPTS}/presearch_db.pl  -o \"-matrix BLOSUM62 -B700 -V700 P=4\" -A ${iso_id}_$ver.fa -s ../mini_db_$proj_name -d /usr/local/db/panda/nraa/nraa_clean -c 10";

    #    system($cmd);
    open( PARA, ">../para_$proj_name" );
    push @{ $para{mini_db}{$iso_id} }, $ver;

    #remove unused key-value pairs in current command.
    $para{mini_db}{option}{o} = '-matrix BLOSUM62 -B700 -V700 P=4';
    &write_proj_para( \%para, \*PARA );
    close PARA;
}

sub write_proj_para {
    use strict;
    my ( $para_r, $PARA ) = @_;
    for my $k ( keys %$para_r ) {
        if ( $k eq 'mini_db' ) {
            print $PARA "mini_db:\n";
            for my $k2 ( keys %{ $$para_r{mini_db} } ) {
                if ( $k2 !~ /\D/ ) {
                    print $PARA "$k2";
                    for ( @{ $$para_r{mini_db}{$k2} } ) {
                        print $PARA "|$_";
                    }
                    print $PARA "\n";
                }
                elsif ( $k2 eq 'option' ) {
                    print $PARA "$k2";
                    for ( keys %{ $$para_r{mini_db}{$k2} } ) {
                        print $PARA "|$_|$$para_r{mini_db}{$k2}{$_}";
                    }
                    print $PARA "\n";
                }
            }
            print $PARA "\/\/\n";
        }
    }
}

sub get_proj_para {
    use strict;
    my ($PARA) = @_;
    my %para;
    while ( my $line = <$PARA> ) {
        chomp $line;
        if ( $line eq 'mini_db:' ) {
            while () {
                my $s = <$PARA>;
                chomp $s;
                last
                  if ( $s eq '//' );
                my ( $k1, $k2 ) = split( /\|/, $s, 2 );
                if ( $k1 eq 'option' ) {
                    %{ $para{mini_db}{option} } = split( /\|/, $s );
                }
                elsif ( $k1 !~ /\D/ ) {
                    @{ $para{mini_db}{$k1} } = split( /\|/, $k2 );
                }
            }
        }
    }
    return %para;
}

sub hmmsearch {
    use strict;
    my ( $proj_id, $iso_id, $ver, $n, $type ) = @_;
    system(
"hmmsearch -Z $n ${iso_id}_${ver}_$type.HMM ../mini_db_$proj_id > hits_${iso_id}_${ver}_${type}_mini"
    );
}

sub full {
    use strict;
    my ( $proj_id, $iso_id, $ver, $method, $cutoff, $dir ) = @_;
    $dir .= $dir . '/'
      if ( $dir ne ''
        && $dir !~ /\/$/ );
    open( FH1, "${dir}hits_${iso_id}_${ver}_${method}_mini" );
    open( FH2, ">${dir}hits_${iso_id}_${ver}_${method}_mini.score" );
    while ( my $line = <FH1> ) {
        print FH2 "$1 $2 $3 $4\n"
          if ( $line =~
/(\S+): domain \d+ of \d+, from (\d+) to (\d+): score (-?\d+\.\d+), E =/
          );
    }
    close FH1;
    close FH2;
    my $status = system(
"$ENV{HMM_SCRIPTS}/create_hmm_alignment.pl -b $cutoff -c ${dir}hits_${iso_id}_${ver}_${method}_mini.score -d ${dir}../mini_db_$proj_id -H ${dir}${iso_id}_${ver}_${method}.HMM -M -S ${dir}hits_${iso_id}_${ver}_${method}.msf > ${dir}hits_${iso_id}_${ver}_${method}.msf.log"
    );
    unlink "${dir}hits_${iso_id}_${ver}_${method}_mini.score";
    unlink "${dir}hits_${iso_id}_${ver}_${method}.msf.log";
    return $status;
}

sub build_hmm {
    use strict;
    my ( $proj_name, $iso_r, $LOG1, $LOG2 ) = @_;
    my %data;
    &init_align( $$iso_r{iso_id}, $$iso_r{init_align} );
    &trim_align( $$iso_r{iso_id}, $$iso_r{trim_align} );
    &hmmbuild( $$iso_r{iso_id}, $$iso_r{hmmbuild} );

    #    &calibrate($$iso_r{iso_id}, $$iso_r{calibrate});
    &mini_db( $proj_name, $$iso_r{iso_id}, $$iso_r{mini_db} );
    &hmmsearch( $proj_name, $$iso_r{iso_id}, $$iso_r{hmmsearch} );
    &full( $proj_name, $$iso_r{iso_id}, $$iso_r{full} );
}

sub get_last_version {

    # $file_r: array reference for file names
    # $s1, $s2: strings before and after version number in file name
    # $ver_r: array reference for sorted version values for output
    # return the largest version value
    use strict;
    my ( $file_r, $s1, $s2, $ver_r ) = @_;
    my @num;
    for (@$file_r) {
        chomp $_;
        if ( $_ =~ /^$s1(\d+)$s2$/ ) {
            push @num, $1;
        }
    }
    @$ver_r = sort { $a <=> $b } @num
      if ($ver_r);
    if ( @num == 0 ) {
        return 0;
    }
    else {
        return "${[sort {$a<=>$b} @num]}[$#num]";
    }
}

sub get_avail_hmm_acc {

    # $mode: '1' to get smallest id available
    use strict;
    my ( $db, $db_proc, $mode ) = @_;
    my $query = "select hmm_acc from $db..hmm2 where hmm_acc like 'TIGR%'";
    my $ref = &do_query( $db_proc, $query );
    my ( @list, $value );
    for (@$ref) {
        push @list, substr( $$_{hmm_acc}, 4 );
    }
    @list = sort { $a <=> $b } @list;
    if ( $mode == 1 ) {
        for my $n ( 1 .. @list - 1 ) {
            if ( $list[$n] - $list[ $n - 1 ] > 1 ) {
                $value = $list[ $n - 1 ] + 1;
                last;
            }
        }
    }
    $value = $list[$#list] + 1
      if ( $value eq '' );
    return "TIGR" . '0' x ( 5 - length($value) ) . "$value";
}
#################################################################################################################
#login SYBYANK as yank, password sybyank
sub replace_egad_acc {
    use strict;
    my ( $seq_r, $omni_cmr, $user, $password, $group, $log_header, $log_r ) =
      @_;
    $log_header .= "--sub replace_egad_acc";
    my %order;
    $order{1} = {
        "OMNI" => 0,
        "SP"   => 1,
        "PIR"  => 2,
        "GP"   => 3
    };
    my $sybyank = &connect_db( 'nraa', 'Sybase', 'SYBYANK', $user, $password );
    for my $n ( 0 .. $$seq_r{number} - 1 ) {
        my ( $match, $db, $id, @tmp ) = ( '', '', '', () );
        if ( $$seq_r{$n}{header} =~ /^([^\|\/\:\s]+)(\||\:)([^\|\/\:\s]+)/ ) {
            $db = $1;
            $id = $3;
        }
        else {
            $$log_r .=
"***Error: $log_header--can not find db and id from sequence header $$seq_r{$n}{header}, changing accession failed\n"
              if ($log_r);
            next;
        }

        #	print "$db, $id\n";
        my $ref = &do_query(
            $sybyank,
"select a2.accession_db, a2.accession_id from accession a1, accession a2 where a1.accession_db=\"$db\" and a1.accession_id=\"$id\" and a1.seq_id=a2.seq_id",
            0,
            $log_header,
            $log_r
        );
        for my $m (@$ref) {
            if ( $$m{accession_db} ne "EGAD" ) {
                $tmp[ $order{$group}{ $$m{accession_db} } ]{db} =
                  $$m{accession_db};
                $tmp[ $order{$group}{ $$m{accession_db} } ]{id} =
                  $$m{accession_id};
            }
        }
        for my $m (@tmp) {
            if (    defined $m
                and $$m{db} ne ''
                and $$m{id} ne '' )
            {

#		print "replace $$seq_r{$n}{db}, $$seq_r{$n}{0}, $$seq_r{$n}{number}, $$seq_r{$n}{prot_id}, $$seq_r{$n}{locus} by ";
                $$seq_r{$n}{db} = (
                    $$m{db} eq "OMNI" && $omni_cmr eq "CMR"
                    ? "CMR"
                    : $$m{db}
                );
                $$seq_r{$n}{0} = $$m{id};
                $$seq_r{$n}{number} = 1
                  if ( $$seq_r{$n}{number} < 1 );
                $$seq_r{$n}{prot_id} = '';
                $$seq_r{$n}{locus}   = '';

#		print "$$seq_r{$n}{db}, $$seq_r{$n}{0}, $$seq_r{$n}{number}, $$seq_r{$n}{prot_id}, $$seq_r{$n}{locus}\n";
                last;
            }
        }
    }
    $sybyank->disconnect
      if ( defined $sybyank );
}

sub download_nraa {

    # $get_seq = 1 for get seq from nraa
    use strict;
    my ( $prot_r, $log_header, $log_r, $get_seq, $user, $password ) = @_;
    $user = "access"
      if ( $user !~ /\S/ );
    $password = "access"
      if ( $password !~ /\S/ );
    my $sybyank = &connect_db( 'nraa', 'Sybase', 'SYBYANK', $user, $password );
    die
      if ( !defined $sybyank );
    for my $acc ( keys %$prot_r ) {

        if ( $acc =~ /^([^\|\/\s]+)\|([^\|\/\s]+)/ ) {
            my $sth =
              $sybyank->prepare("exec sybyank_fasta_header \"$1\", \"$2\"");
            $sth->execute();
            my $ref = $sth->fetchall_arrayref( {} );
            $sth->finish;
            if ( @$ref > 0 ) {
                my $seq_id = $$ref[0]{seq_id};
                $$prot_r{$acc}{number}    = @$ref;
                $$prot_r{$acc}{first_acc} = $1
                  if (
                    $$ref[0]{accession_info} =~ /^([^\|\/\s]+\|[^\|\/\s]+)/ );
                for my $n ( 0 .. @$ref - 1 ) {
                    (
                        $$prot_r{$acc}{$n}{header},
                        $$prot_r{$acc}{$n}{prot_info}
                    ) = split( /\s+/, $$ref[$n]{accession_info}, 2 );
                    $$prot_r{$acc}{$n}{header}    =~ s/^\s+|\s+$//g;
                    $$prot_r{$acc}{$n}{prot_info} =~ s/^\s+|\s+$//g;
                    my @header = split /\|/, $$prot_r{$acc}{$n}{header};
                    $$prot_r{$acc}{$n}{db} = $header[0];
                    $$prot_r{$acc}{$n}{id} = $header[1];
                    if ( $$prot_r{$acc}{name} !~ /\S/ ) {
                        $$prot_r{$acc}{name} = $$prot_r{$acc}{$n}{prot_info};
                        ( $$prot_r{$acc}{genus}, $$prot_r{$acc}{species} ) =
                          split( /\s+/, $1, 2 )
                          if ( $$prot_r{$acc}{name} =~ s/\{([^\}]+)\}// );
                    }
                    ( $$prot_r{$acc}{genus}, $$prot_r{$acc}{species} ) =
                      split( /\s+/, $1, 2 )
                      if (  $$prot_r{$acc}{species} !~ /S/
                        and $$prot_r{$acc}{$n}{prot_info} =~ /\{([^\}]+)\}/ );
                    ( $$prot_r{$acc}{genus}, $$prot_r{$acc}{species} ) =
                      split( /\s+/, $1, 2 )
                      if (  $$prot_r{$acc}{species} !~ /S/
                        and $$prot_r{$acc}{$n}{db} eq "PIR"
                        and $$prot_r{$acc}{$n}{prot_info} =~
                        /\(([^\)]+)\)\s*$/ );
                    $$prot_r{$acc}{other_acc}{"$header[0]|$header[1]"} = 1
                      if ( $acc ne "$header[0]|$header[1]" );
                }
                if ( $get_seq == 1 ) {
                    my $sth =
                      $sybyank->prepare("sybyank_fasta_sequence $seq_id");
                    $sth->execute();
                    my $ref = $sth->fetchall_arrayref( {} );
                    $sth->finish;
                    for my $n (@$ref) {
                        $$prot_r{$acc}{seq} .= uc "$$n{seq_line}";
                    }
                    $$prot_r{$acc}{seq} =~ s/[^A-Za-z]//g;
                    $$prot_r{$acc}{length} = length( $$prot_r{$acc}{seq} );
                }
            }
        }
    }
    $sybyank->disconnect
      if ( defined $sybyank );
}

sub format_sequence {

# if stockholm header should not be printed, set $$seq_r{stockholm}='' before calling this subroutine
# empty sequence are deleted
# seq_r: hash reference for single or multiple sequences
# format: format of sequence for output
# seq_num: the number of a sequence for formatting (0 .. number-1). otherwise formatting all sequences.
# gap: 0 using ungapped sequence in seq_r hash, 1 otherwise
# LOG1: log file for short output
# LOG2: log file for long output
    use strict;
    my ( $seq_r, $format, $seq_num, $length, $gap, $header, $LOG1, $LOG2 ) = @_;
    my $index = 0;
    for my $n ( 0 .. $$seq_r{number} - 1 ) {
        if ( defined $$seq_r{$n}
            && $$seq_r{$n} ne '' )
        {
            $$seq_r{$index} = $$seq_r{$n}
              if ( $index != $n );
            ++$index;
        }
    }
    $$seq_r{number} = $index;
    my ( $r1, $r2, $result );
    if (   $seq_num ne ''
        && $seq_num !~ /\D/
        && $seq_num >= 0
        && $seq_num <= $$seq_r{number} - 1 )
    {
        $r1 = $seq_num;
        $r2 = $seq_num;
    }
    else {
        $r1 = 0;
        $r2 = $$seq_r{number} - 1;
    }

    #format to mul
    if ( $format eq "mul" ) {

        #format STOCKHOLM information
        if ( $$seq_r{stockholm} ne '' ) {
            $result .= "# STOCKHOLM\n";
            $result .=
              "#=GF hmmbuild_acc $$seq_r{stockholm}{GF}{hmmbuild_acc}\n"
              if ( defined $$seq_r{stockholm}{GF}{hmmbuild_acc}
                && $$seq_r{stockholm}{GF}{hmmbuild_acc} ne '' );
            $result .= "#=GF project_id $$seq_r{stockholm}{GF}{project_id}\n"
              if ( defined $$seq_r{stockholm}{GF}{project_id}
                && $$seq_r{stockholm}{GF}{project_id} ne '' );
            $result .= "#=GF iso_id $$seq_r{stockholm}{GF}{iso_id}\n"
              if ( defined $$seq_r{stockholm}{GF}{iso_id}
                && $$seq_r{stockholm}{GF}{iso_id} ne '' );
            for
              my $n ( 0 .. $$seq_r{stockholm}{GF}{hmmbuild_method}{number} - 1 )
            {
                $result .=
"#=GF hmmbuild_method $$seq_r{stockholm}{GF}{hmmbuild_method}{$n}\n";
            }
            for my $k ( 0 .. $$seq_r{stockholm}{GF}{mini_db}{number} - 1 ) {
                my $acc  = $$seq_r{stockholm}{GF}{mini_db}{$k}{acc};
                my $seq  = $$seq_r{stockholm}{GF}{mini_db}{$k}{seq};
                my $cmd  = $$seq_r{stockholm}{GF}{mini_db}{$k}{cmd};
                my $name = $$seq_r{stockholm}{GF}{mini_db}{$k}{db_name};
                my $size = $$seq_r{stockholm}{GF}{mini_db}{$k}{db_size};
                my $date = $$seq_r{stockholm}{GF}{mini_db}{$k}{db_date};
                $result .=
"#=GF MINI_DB ACC $acc SEQ $seq CMD \"$cmd\" DB $name $size \"$date\"\n";
            }
        }

        #format sequence information
        foreach my $n ( $r1 .. $r2 ) {
            my $sequence =
              ($gap)
              ? $$seq_r{$n}{seq_gap}
              : $$seq_r{$n}{seq};
            next
              unless ( $sequence =~ /\w/ );
            my $line = $$seq_r{$n}{db}
              if ( $$seq_r{$n}{db} ne '' );
            $line .= '|' . $$seq_r{$n}{prot_id}
              if ( defined $$seq_r{$n}{prot_id}
                && $$seq_r{$n}{prot_id} ne '' );
            $line .= '|' . $$seq_r{$n}{locus}
              if ( defined $$seq_r{$n}{locus}
                && $$seq_r{$n}{locus} ne '' );
            for my $m ( 0 .. $$seq_r{$n}{number} - 1 ) {
                last
                  if ( length($line) + 1 + length( $$seq_r{$n}{$m} ) > 30 );
                $line .= '|' . $$seq_r{$n}{$m};    # if($$seq_r{$n}{$m} ne '');
            }
            $line =~ s/^\|//;
            $line .= '/' . $$seq_r{$n}{lend} . "-" . $$seq_r{$n}{rend}
              if ( defined $$seq_r{$n}{lend}
                && defined $$seq_r{$n}{rend}
                && $$seq_r{$n}{lend} > 0
                && $$seq_r{$n}{rend} > 0 );
            my $tmp = (
                length($line) < 40
                ? 40 - length($line)
                : 1
            );
            $line   .= ( " " x $tmp ) . $sequence . "\n";
            $result .= $line;
#### note that it may fail if the header is too long, eg for GP|XXXXXX
        }
    }

    #format to fasta
    elsif ( $format eq "fasta" ) {
        for my $n ( $r1 .. $r2 ) {
            my $sequence =
              ($gap)
              ? $$seq_r{$n}{seq_gap}
              : $$seq_r{$n}{seq};
            next
              unless ( $sequence =~ /\w/ );
            $length = 60
              if ( $length eq '' );
            my $line = '>';
            $line .= $$seq_r{$n}{db}
              if ( $$seq_r{$n}{db} ne '' );
            $line .= '|' . $$seq_r{$n}{prot_id}
              if ( defined $$seq_r{$n}{prot_id}
                && $$seq_r{$n}{prot_id} ne '' );
            $line .= '|' . $$seq_r{$n}{locus}
              if ( defined $$seq_r{$n}{locus}
                && $$seq_r{$n}{locus} ne '' );
            if ( defined $$seq_r{$n}{number}
                && $$seq_r{$n}{number} =~ /^\d+$/ )
            {

                for ( 0 .. $$seq_r{$n}{number} - 1 ) {
                    $line .= '|' . $$seq_r{$n}{$_};
                }
            }
            $line =~ s/^>\|/>/;
            $line .= '/' . $$seq_r{$n}{lend} . "-" . $$seq_r{$n}{rend}
              if ( defined $$seq_r{$n}{lend}
                && defined $$seq_r{$n}{rend}
                && $$seq_r{$n}{lend} > 0
                && $$seq_r{$n}{rend} > 0 );
            $line .= ' ' . $$seq_r{$n}{comment}
              if ( $$seq_r{$n}{comment} );
            $line .= "\n";
            while () {
                if ( length($sequence) > $length ) {
                    $sequence =~ s/^([a-zA-Z\-.]{$length})//;
                    $line .= $1 . "\n";
                }
                else {
                    $line .= $sequence . "\n";
                    last;
                }
            }
            $result .= $line;
        }
    }
    elsif ( $format eq "msf" ) {
        my $head_len = 40;
        my ( %new_acc, %tmp_seq );
        my $len = length( $$seq_r{0}{seq_gap} );
        $result = "PileUp\n\n   MSF:  $len  Type: P  Check: 1111  ..\n\n"
          ;    #added DUMMY checksum for HMMER 2.2g hmmbuild. DHH.

        #	$result = "PileUp\n\n\n\n   MSF:  $len  Type: P  Check:   ..\n\n";
        for my $n ( $r1 .. $r2 ) {
            $tmp_seq{$n} = $$seq_r{$n}{seq_gap};
            $new_acc{$n} = $$seq_r{$n}{db}
              if ( $$seq_r{$n}{db} ne '' );
            $new_acc{$n} .= '|' . $$seq_r{$n}{prot_id}
              if ( $$seq_r{$n}{prot_id} ne '' );
            $new_acc{$n} .= '|' . $$seq_r{$n}{locus}
              if ( $$seq_r{$n}{locus} ne '' );
            for my $m ( 0 .. $$seq_r{$n}{number} - 1 ) {
                last
                  if (
                    length( $new_acc{$n} ) + 1 + length( $$seq_r{$n}{$m} ) >
                    30 );
                $new_acc{$n} .=
                  '|' . $$seq_r{$n}{$m};    # if($$seq_r{$n}{$m} ne '');
            }
            $new_acc{$n} =~ s/^\|//;
            $new_acc{$n} =~ s/\|+$//;
            $new_acc{$n} .= '/' . $$seq_r{$n}{lend} . "-" . $$seq_r{$n}{rend}
              if ( defined $$seq_r{$n}{lend}
                && defined $$seq_r{$n}{rend}
                && $$seq_r{$n}{lend} > 0
                && $$seq_r{$n}{rend} > 0 );
            $new_acc{$n} = substr( $new_acc{$n}, 0, $head_len );
            $new_acc{$n} .= ' ' x ( $head_len - length( $new_acc{$n} ) )
              if ( $head_len > length( $new_acc{$n} ) );
            $result .=
              " Name: $new_acc{$n} Len:  $len  Check:   0  Weight:  1.0\n";
        }
        $result .= "\n//\n\n\n";
        while ( length( $tmp_seq{$r1} ) > 0 ) {
            for my $n ( $r1 .. $r2 ) {
                $result .= $new_acc{$n};
                for my $m ( 0 .. 4 ) {
                    if ( length( $tmp_seq{$n} ) >= 10 ) {
                        $result .= ' ' . $1
                          if ( $tmp_seq{$n} =~ s/^([a-zA-Z\-.]{10})// );
                    }
                    else {
                        $result .= ' ' . $tmp_seq{$n};
                        $tmp_seq{$n} = '';
                    }
                }
                $result .= "\n";
            }
            $result .= "\n\n";
        }
        chomp $result;
        chomp $result;
    }
    else {
        print $LOG1
"SEQ_FORMAT_ERROR\t$header\tunknown format: $format in 'format_sequence'\n"
          if ($LOG1);
        print $LOG2
"SEQ_FORMAT_ERROR\t$header\tunknown format: $format in 'format_sequence'\n"
          if ($LOG2);
    }
    chomp $result;
    return $result;
}
#######################################################################################################################################
sub check_data {
    my ( $info_r, $data_r, $flag_r, $fatal_r, $LOG1, $LOG2 ) = @_;
    if ( $$data_r{hmm_seed}{ori} ) {
        print $LOG2 "$$data_r{header} check hmm_seed .........\n";
        &check_MA(
            $$data_r{hmm_seed},      'mul',
            \%{ $$data_r{protein} }, $$info_r{db_proc},
            $$info_r{blastp},        \$$flag_r{hmm_seed},
            $fatal_r,                $$data_r{header},
            $LOG1,                   $LOG2
        );
    }
    if ( $$data_r{align_seed}{ori} ) {
        print $LOG2 "$$data_r{header} check align_seed .........\n";
        &check_MA(
            $$data_r{align_seed},    'fasta',
            \%{ $$data_r{protein} }, $$info_r{db_proc},
            $$info_r{blastp},        \$$flag_r{align_seed},
            $fatal_r,                $$data_r{header},
            $LOG1,                   $LOG2
        );
    }
    if ( $$data_r{align_full}{ori} ) {
        print $LOG2 "$$data_r{header} check align_full .........\n";
        &check_MA(
            $$data_r{align_full},    'fasta',
            \%{ $$data_r{protein} }, $$info_r{db_proc},
            $$info_r{blastp},        \$$flag_r{align_full},
            $fatal_r,                $$data_r{header},
            $LOG1,                   $LOG2
        );
    }
    if (   $$data_r{hmm_seed}{ori}
        && $$data_r{align_seed}{ori} )
    {
        print $LOG2
"$$data_r{header} check hmm_seed(MA1) with align_seed(MA2) .........\n";
        &check_MA_MA( $$data_r{hmm_seed}, $$data_r{align_seed}, 'eq', $fatal_r,
            $$data_r{header}, $LOG1, $LOG2 );
    }
    if (   $$data_r{align_full}{ori}
        && $$data_r{align_seed}{ori} )
    {
        print $LOG2
"$$data_r{header} check align_seed(MA1) with align_full(MA2) .........\n";
        &check_MA_MA( $$data_r{align_seed}, $$data_r{align_full}, 'le',
            $fatal_r, $$data_r{header}, $LOG1, $LOG2 );
    }
    if (   $$data_r{hmm_seed}{ori}
        && $$data_r{align_full}{ori} )
    {
        print $LOG2
"$$data_r{header} check hmm_seed(MA1) with align_full(MA2) .........\n";
        &check_MA_MA( $$data_r{hmm_seed}, $$data_r{align_full}, 'le', $fatal_r,
            $$data_r{header}, $LOG1, $LOG2 );
    }
    if (   $$data_r{iso_link}
        && $$data_r{align_full}{ori} )
    {
        print $LOG2 "check lend and rend in iso_link table.............\n";
        &check_lend_rend_in_iso_link( $$data_r{align_full}, $$data_r{iso_link},
            $$data_r{iso_id}, \$$flag_r{iso_link}, $fatal_r, $$data_r{header},
            $LOG1, $LOG2 );
    }
}

sub check_MA {
    my (
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    ) = @_;
    &check_MA_content(
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    );
    &check_locus_in_mul_format(
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    ) if ( $type eq 'mul' );
    &check_MA_prot_id(
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    );
    &check_lend_rend_in_MA(
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    );
}

sub check_hmm2 {

    #suggestion: remove "like .." in query, replace * in query to items.
    my ( $info_r, $LOG1, $LOG2 ) = @_;
    my ( $ref, $query, $num );
    $query =
      "select * from egad..hmm2 where hmm_acc like 'TIGR%' and is_current = 1";
    $ref = &do_query( $$info_r{db_proc}, $query );
    $num = @$ref;
    for my $n1 ( 0 .. $num - 2 ) {
        for my $n2 ( $n1 + 1 .. $num - 1 ) {
            if ( $$ref[$n1]{hmm_name} eq $$ref[$n2]{hmm_name} ) {
                print $LOG2
"$$ref[$n1]{hmm_name}, $$ref[$n1]{hmm_acc}, $$ref[$n2]{hmm_acc}********************************************************************\n";
                for ( keys %{ $$ref[$n1] } ) {
                    print LOG2
                      "key:$_, value1:$$ref[$n1]{$_}, value2:$$ref[$n2]{$_}\n"
                      if ( $$ref[$n1]{$_} ne $$ref[$n2]{$_} );
                }
            }
            if ( $$ref[$n1]{hmm_acc} eq $$ref[$n2]{hmm_acc} ) {
                print $LOG2
"$$ref[$n1]{hmm_acc}, $$ref[$n2]{hmm_acc}********************************************************************\n";
                for ( keys %{ $$ref[$n1] } ) {
                    print LOG2
                      "key:$_, value1:$$ref[$n1]{$_}, value2:$$ref[$n2]{$_}\n"
                      if ( $$ref[$n1]{$_} ne $$ref[$n2]{$_} );
                }
            }
            if ( $$ref[$n1]{iso_id} eq $$ref[$n2]{iso_id} ) {
                print $LOG2
"$$ref[$n1]{iso_id}, $$ref[$n1]{hmm_acc}, $$ref[$n2]{hmm_acc}********************************************************************\n";
                for ( keys %{ $$ref[$n1] } ) {
                    print LOG2
                      "key:$_, value1:$$ref[$n1]{$_}, value2:$$ref[$n2]{$_}\n"
                      if ( $$ref[$n1]{$_} ne $$ref[$n2]{$_} );
                }
            }
        }
    }
    $query = "select * from egad..isology where iso_acc like 'TIGR%'";
    $ref   = &do_query( $$info_r{db_proc}, $query );
    $num   = @$ref;
    for my $n1 ( 0 .. $num - 2 ) {
        for my $n2 ( $n1 + 1 .. $num - 1 ) {
            if ( $$ref[$n1]{iso_acc} eq $$ref[$n2]{iso_acc} ) {
                print $LOG2
"$$ref[$n1]{iso_acc}, $$ref[$n1]{iso_id}, $$ref[$n2]{iso_id}********************************************************************\n";
                for ( keys %{ $$ref[$n1] } ) {
                    print LOG2
                      "key:$_, value1:$$ref[$n1]{$_}, value2:$$ref[$n2]{$_}\n"
                      if ( $$ref[$n1]{$_} ne $$ref[$n2]{$_} );
                }
            }
        }
    }
    $query = "select iso_id from egad..isology where iso_acc like 'TIGR%'";
    $ref   = &do_query( $$info_r{db_proc}, $query );
    $num   = @$ref;
    for my $n1 ( 0 .. $num - 1 ) {
        $query =
"select i.iso_acc i_acc, h.hmm_acc h_acc from isology i, hmm2 h where i.iso_id = h.iso_id and i.iso_id = $$ref[$n1]{iso_id}";
        my $ref1 = &do_query( $$info_r{db_proc}, $query );
        if ( $$ref1[0]{i_acc} ne $$ref1[0]{h_acc} ) {
            print $LOG2
"$$ref[$n1]{iso_id}, $$ref1[0]{i_acc}, $$ref1[0]{h_acc}********************************************************************\n";
        }
    }
    $query = "select id from egad..hmm2 where hmm_acc like 'TIGR%'";
    $ref   = &do_query( $$info_r{db_proc}, $query );
    $num   = @$ref;
    for my $n1 ( 0 .. $num - 1 ) {
        $query =
"select i.iso_acc i_acc, h.hmm_acc h_acc from isology i, hmm2 h where i.iso_id = h.iso_id and h.id = $$ref[$n1]{id}";
        my $ref1 = &do_query( $$info_r{db_proc}, $query );
        if ( $$ref1[0]{i_acc} ne $$ref1[0]{h_acc} ) {
            print $LOG2
"$$ref[$n1]{id}, $$ref1[0]{i_acc}, $$ref1[0]{h_acc}********************************************************************\n";
        }
    }
}

sub check_lend_rend_in_iso_link {
    use strict;
    my ( $ma_r, $iso_link_r, $iso_id, $flag_r, $fatal_r, $header, $LOG1, $LOG2 )
      = @_;
    my @new_iso_link;
    &construct_iso_link( $ma_r, $iso_id, \@new_iso_link );

    #    for(@new_iso_link) {
    #	print $LOG2 "$_->{prot_id}, $_->{iso_id}, $_->{lend}, $_->{rend}\n";
    #    }
    for my $k1 (@new_iso_link) {
        my $found = 0;
        for my $k2 (@$iso_link_r) {
            if (   $$k1{prot_id} == $$k2{prot_id}
                && $$k1{lend} == $$k2{lend}
                && $$k1{rend} == $$k2{rend} )
            {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            print $LOG1
"BAD_ISOLINK\t$header\tmissing prot_id or wrong lend/rend. rebuild isolink. there were ${\scalar(@$iso_link_r)} lines ";
            print $LOG2
"BAD_ISOLINK\t$header\tmissing prot_id or wrong lend/rend. rebuild isolink. there were ${\scalar(@$iso_link_r)} lines ";
            @$iso_link_r = @new_iso_link;
            $$flag_r .= "the iso_link table is re-constructed\n";
            print $LOG1 "there are ${\scalar(@$iso_link_r)} lines now\n";
            print $LOG2 "there are ${\scalar(@$iso_link_r)} lines now\n";
            last;
        }
    }
    for my $k1 (@$iso_link_r) {
        my $found = 0;
        for my $k2 (@new_iso_link) {
            if (   $$k1{prot_id} == $$k2{prot_id}
                && $$k1{lend} == $$k2{lend}
                && $$k1{rend} == $$k2{rend} )
            {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            print $LOG1
"BAD_ISOLINK\t$header\tmissing prot_id or wrong lend/rend. rebuild isolink. there were ${\scalar(@$iso_link_r)} lines ";
            print $LOG2
"BAD_ISOLINK\t$header\tmissing prot_id or wrong lend/rend. rebuild isolink. there were ${\scalar(@$iso_link_r)} lines ";
            @$iso_link_r = @new_iso_link;
            $$flag_r .= "the iso_link table is re-constructed\n";
            print $LOG1 "there are ${\scalar(@$iso_link_r)} lines now\n";
            print $LOG2 "there are ${\scalar(@$iso_link_r)} lines now\n";
            last;
        }
    }
}

sub construct_iso_link {

    # $ma_r is hash reference for sequences
    # $iso_id is not used here
    # iso_link_r is the array reference that hold the lnking information
    use strict;
    my ( $ma_r, $iso_id, $iso_link_r, $status_vocab_id ) = @_;
    my $m = -1;
    foreach my $k ( 0 .. $$ma_r{number} - 1 ) {
        if ( $$ma_r{$k}{header} =~ /^([^\/\|\s]+)\|([^\/\|\s]+)/ ) {
            ++$m;
            my $acc_db = $1;
            my $acc_id = $2;
            $$iso_link_r[$m]{acc_db}          = $acc_db;
            $$iso_link_r[$m]{acc_id}          = $acc_id;
            $$iso_link_r[$m]{lend}            = $$ma_r{$k}{lend};
            $$iso_link_r[$m]{rend}            = $$ma_r{$k}{rend};
            $$iso_link_r[$m]{seq_id}          = 0;
            $$iso_link_r[$m]{status_vocab_id} = $status_vocab_id;

            if ( $$ma_r{stockholm}{GF}{score}{"EGAD|$$ma_r{$k}{prot_id}"} ) {
                $$iso_link_r[$m]{score} =
                  $$ma_r{stockholm}{GF}{score}{"EGAD|$$ma_r{$k}{prot_id}"}{b};
                $$iso_link_r[$m]{e_val} =
                  $$ma_r{stockholm}{GF}{score}{"EGAD|$$ma_r{$k}{prot_id}"}{e};
            }
            elsif ($$ma_r{$k}{score} ne ''
                or $$ma_r{$k}{e_val} ne '' )
            {
                $$iso_link_r[$m]{score} = $$ma_r{$k}{score};
                $$iso_link_r[$m]{e_val} = $$ma_r{$k}{e_val};
            }
        }
    }
}

sub load_iso_link {

    # $revise: 1 to revise the database, 0 otherwise.
    # hash reference $iso_link_ref: holding linking info
    use strict;
    my ( $revise, $iso_id, $iso_link_ref, $db, $db_proc, $LOG, $DB_LOG,
        $header ) = @_;
    &do_query( $db_proc, "delete from $db..iso_link where iso_id = $iso_id", 1 )
      if ($revise);
    for (@$iso_link_ref) {
        &convert_null( '', 'null', \$$_{lend}, \$$_{rend}, \$$_{score},
            \$$_{e_val} );
        my $query =
"insert into $db..iso_link (seq_id, iso_id, lend, rend, score, e_val, acc_db, acc_id) values ($$_{seq_id}, $iso_id, $$_{lend}, $$_{rend}, $$_{score}, $$_{e_val}, \"$$_{acc_db}\", \"$$_{acc_id}\")";
        &convert_null( 'null', '', \$$_{lend}, \$$_{rend}, \$$_{score},
            \$$_{e_val} );
        &do_query( $db_proc, $query, 1 )
          if ($revise);
    }
}

sub load_align_link {

    # $revise: 1 to revise the database, 0 otherwise.
    # hash reference $align_link_ref: holding linking info
    use strict;
    my (
        $revise, $align_id, $align_link_ref,
        $db,     $db_proc,  $LOG,
        $DB_LOG, $header,   $status_vocab_id
    ) = @_;
    my $query =
"delete $db..align_link where align_id = $align_id and status_vocab_id = $status_vocab_id";
    print $LOG "$header delete query = $query\n"
      if ($LOG);
    if ($revise) {
        &do_query( $db_proc, $query, 1 );
        print $DB_LOG "$query\n"
          if ($DB_LOG);
    }
    for (@$align_link_ref) {
        $$_{prot_id} = "null"
          if ( $$_{prot_id} eq '' );
        $$_{lend} = "null"
          if ( $$_{lend} eq '' );
        $$_{rend} = "null"
          if ( $$_{rend} eq '' );
        $$_{score} = "null"
          if ( $$_{score} eq '' );
        $$_{e_val} = "null"
          if ( $$_{e_val} eq '' );
        $$_{status_vocab_id} = "null"
          if ( $$_{status_vocab_id} eq '' );
###inserted check --> any data with NULL prot_id is trash and should not be entered into align_link table ... 07/19/01
        if ( $$_{prot_id} ne "null" ) {
            $query =
"insert $db..align_link (align_id, prot_id, lend, rend, score, e_val, status_vocab_id) values ($align_id, $$_{prot_id}, $$_{lend}, $$_{rend}, $$_{score}, $$_{e_val}, $$_{status_vocab_id})";
            print $LOG "$header insert query = $query\n"
              if ($LOG);
        }
        else {
            print $LOG
"ERROR! $header insert query = $query\n prot_id cannot equal NULL.\n"
              if ($LOG);
        }
        if ($revise) {
            &do_query( $db_proc, $query, 1 );
            print $DB_LOG "$query\n"
              if ($DB_LOG);
        }
        $$_{lend} = ''
          if ( $$_{lend} eq 'null' );
        $$_{rend} = ''
          if ( $$_{rend} eq 'null' );
        $$_{score} = ""
          if ( $$_{score} eq 'null' );
        $$_{e_val} = ""
          if ( $$_{e_val} eq 'null' );
        $$_{status_vocab_id} = ""
          if ( $$_{status_vocab_id} eq 'null' );
    }
}

sub check_MA_MA {

    #done
    use strict;
    my ( $ma1_r, $ma2_r, $cond, $fatal_r, $header, $LOG1, $LOG2 ) = @_;
    if (   $$ma1_r{number} != $$ma2_r{number}
        && $cond eq 'eq' )
    {
        print $LOG1
"MA_MA_ERROR\t$header\tthe number of sequences in the two MAs are not equal\n";
        print $LOG2
"MA_MA_ERROR\t$header\tthe number of sequences in the two MAs are not equal\n";
    }
    elsif ($$ma1_r{number} > $$ma2_r{number}
        && $cond eq 'le' )
    {
        print $LOG1
"MA_MA_ERROR\t$header\tthe number of sequences in MA1 is larger than that in MA2\n";
        print $LOG2
"MA_MA_ERROR\t$header\tthe number of sequences in MA1 is larger than that in MA2\n";
    }
    elsif ($$ma1_r{number} < $$ma2_r{number}
        && $cond eq 'ge' )
    {
        print $LOG1
"MA_MA_ERROR\t$header\tthe number of sequences in MA1 is smaller than that in MA2\n";
        print $LOG2
"MA_MA_ERROR\t$header\tthe number of sequences in MA1 is smaller than that in MA2\n";
    }
    foreach my $k1 ( 0 .. $$ma1_r{number} - 1 ) {
        my $found;
        my $s1 = uc $$ma1_r{$k1}{seq};
        if ( $s1 eq '' ) {
            print $LOG1 "SEQ_ERROR\t$header\tempty sequence $k1 in MA1\n";
            print $LOG2 "SEQ_ERROR\t$header\tempty sequence $k1 in MA1\n";
            next;
        }
        foreach my $k2 ( 0 .. $$ma2_r{number} - 1 ) {
            $found = 1;
            my $s2 = uc $$ma2_r{$k2}{seq};
            if ( $s2 ne '' ) {
                for ( keys %{ $$ma1_r{$k1} }, keys %{ $$ma2_r{$k2} } ) {

        #		    print $LOG2 "$k1, $k2, $_, $$ma1_r{$k1}{$_}, $$ma2_r{$k2}{$_}\n";
                    if (   $$ma1_r{$k1}{$_} ne $$ma2_r{$k2}{$_}
                        && $_ =~ /^(db|prot_id|locus|\d+)$/ )
                    {
                        $found = 0;
                        last;
                    }
                }
                if ($found) {
                    my ( $success, $a_left, $a_right, $b_left, $b_right ) =
                      &match_two_string( $s1, $s2 );
                    if (
                        $success
                        && (   ( $a_right - $a_left + 1 ) / length($s1) < 0.80
                            && ( $b_right - $b_left + 1 ) / length($s2) < 0.80 )
                      )
                    {
                        print $LOG2
"$k1, $k2, $_, ${\(length($s1))}, ${\(length($s2))}, $a_left, $a_right, $b_left, $b_right\n";
                        $success = 0;
                    }
                    if ($success) {
                        print $LOG2
"SEQ_MATCH_WARNING\t$header\tthe seq $k1 in MA1 is not identical to seq $k2 in MA2. length are (${\(length($s1))}, ${\(length($s2))}), coordinates are ($a_left, $a_right, $b_left, $b_right)\n"
                          if ( $a_left != 1
                            or $a_right != length($s1)
                            or $b_left != 1
                            or $b_right != length($s2) );
                        last;
                    }
                }
            }
            else {
                print $LOG1 "SEQ_ERROR\t$header\tempty sequence $k2 in MA2\n";
                print $LOG2 "SEQ_ERROR\t$header\tempty sequence $k2 in MA2\n";
            }
        }
        unless ($found) {
            print $LOG1
"MA_MA_ERROR\t$header\tthe sequence MA1 $k1 is not found in MA2\n";
            print $LOG2
"MA_MA_ERROR\t$header\tthe sequence MA1 $k1 is not found in MA2\n";

#	    $$fatal_r .= "Mismatch of sequences bwtween hmm_seed/align_seed/align_full\n";
        }
    }
}

sub match_two_string {

 # take $a and $b, return the coordinates of matching section on both sequences.
 #done
    use strict;
    my ( $a, $b ) = @_;
    my ( $success, $s1_l, $s1_r, $s2_l, $s2_r );
    my $swap;

    #put $b as shorter sequence
    if ( length($a) < length($b) ) {
        ( $a, $b ) = ( $b, $a );
        $swap = 1;
    }

    #initial trial, truncate both ends evenly
    for ( 0 .. 8 ) {
        if ( length($b) >= $_ * 2 + 1
            && index( $a, substr( $b, $_, length($b) - $_ * 2 ) ) >= 0 )
        {
            $success = 1;
            ( $s1_l, $s1_r, $s2_l, $s2_r ) =
              &extend_ends( $a, $b, substr( $b, $_, length($b) - $_ * 2 ) );
            last;
        }
    }

    # recursive trial for matching
    unless ($success) {
        ( $success, my $c ) = &recursive_match( $a, $b, $b );
        ( $s1_l, $s1_r, $s2_l, $s2_r ) = &extend_ends( $a, $b, $c )
          if ($success);
    }
    if ($swap) {
        return ( $success, $s2_l, $s2_r, $s1_l, $s1_r );
    }
    else {
        return ( $success, $s1_l, $s1_r, $s2_l, $s2_r );
    }
}

sub extend_ends {

    #extend matching section $c to maximum length
    #done
    use strict;
    my ( $a, $b, $c ) = @_;
    my $s2_l = index( $b, $c ) + 1;
    my $s2_r = index( $b, $c ) + length($c);
    while ( index( $a, substr( $b, $s2_l - 1, $s2_r - $s2_l + 1 ) ) >= 0 ) {
        last
          if ( $s2_r > length($b) );
        ++$s2_r;
    }
    --$s2_r;
    while ( index( $a, substr( $b, $s2_l - 1, $s2_r - $s2_l + 1 ) ) >= 0 ) {
        last
          if ( $s2_l < 1 );
        --$s2_l;
    }
    ++$s2_l;
    my $s1_l = index( $a, substr( $b, $s2_l - 1, $s2_r - $s2_l + 1 ) ) + 1;
    my $s1_r = $s1_l + $s2_r - $s2_l;
    return ( $s1_l, $s1_r, $s2_l, $s2_r );
}

sub recursive_match {

    #recursive matching, split $c to two parts.
    #done
    use strict;
    my ( $a, $b, $c ) = @_;
    if ( index( $a, $c ) >= 0 ) {
        my ( $sa_l, $sa_r, $sb_l, $sb_r ) = &extend_ends( $a, $b, $c );
        return ( 1, substr( $b, $sb_l - 1, $sb_r - $sb_l + 1 ) );
    }
    elsif ( length($c) < 6 ) {
        return ( 0, '' );
    }
    else {
        my $half = int( length($c) / 2 + 0.5 );
        my ( $success1, $c1 ) =
          &recursive_match( $a, $b, substr( $c, 0, $half ) );
        my ( $success2, $c2 ) = &recursive_match( $a, $b, substr( $c, $half ) );
        if (   $success1
            && $success2
            && length($c1) <= length($c2) )
        {
            return ( 1, $c2 );
        }
        elsif ($success1
            && $success2
            && length($c1) > length($c2) )
        {
            return ( 1, $c1 );
        }
        elsif ( $success1
            && !$success2 )
        {
            return ( 1, $c1 );
        }
        elsif ( !$success1
            && $success2 )
        {
            return ( 1, $c2 );
        }
        else {
            return ( 0, '' );
        }
    }
}

sub check_lend_rend_in_MA {

    #done
    use strict;
    my (
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    ) = @_;
    foreach my $k ( 0 .. $$ma_r{number} - 1 ) {

# check sequences in MA with full length protein sequence, and calculate lend and rend
        if ( $$ma_r{$k}{prot_id} ne '' ) {
            my $prot_id = $$ma_r{$k}{prot_id};
            my $a       = uc $$prot_r{$prot_id}{prot_seq};
            my $b       = uc $$ma_r{$k}{seq};
            if ( $a eq '' ) {
                print $LOG1
"SEQ_ERROR\t$header\tempty sequence in protein for prot_id $prot_id\n"
                  if ($LOG1);
                print $LOG2
"SEQ_ERROR\t$header\tempty sequence in protein for prot_id $prot_id\n"
                  if ($LOG2);
                next;
            }
            if ( $b eq '' ) {
                print $LOG1
"SEQ_ERROR\t$header\tempty sequence in alignment for prot_id $prot_id\n"
                  if ($LOG1);
                print $LOG2
"SEQ_ERROR\t$header\tempty sequence in alignment for prot_id $prot_id\n"
                  if ($LOG2);
                next;
            }
            my ( $success, $a_left, $a_right, $b_left, $b_right ) =
              &match_two_string( $a, $b );
            if ( $success
                && ( ( $b_right - $b_left + 1 ) / length($b) < 0.85 ) )
            {
                print $LOG1
"SEQ_MATCH_ERROR\t$header\ttoo much truncation. length, lend, rend on seq $k in MA are ${\(length$b)}, $b_left, $b_right\n"
                  if ($LOG1);
                print $LOG2
"SEQ_MATCH_ERROR\t$header\ttoo much truncation. length, lend, rend on seq $k in MA are ${\(length$b)}, $b_left, $b_right\n"
                  if ($LOG2);
                $success = 0;
            }
            if ($success) {
                my $c_left  = $b_left - 1;
                my $c_right = length($b) - $b_right;
                if (   $c_left > 0
                    || $c_right > 0 )
                {
                    print $LOG1
"SEQ_MATCH_WARNING\t$header\tthe seq $k in MA is truncated ($c_left, $c_right) times to match that from protein table\n"
                      if ($LOG1);
                    print $LOG2
"SEQ_MATCH_WARNING\t$header\tthe seq $k in MA is truncated ($c_left, $c_right) times to match that from protein table\n"
                      if ($LOG2);
                }
                my $cal_left  = $a_left - $b_left + 1;
                my $cal_right = $a_right + length($b) - $b_right;
                if (
                       0 < $cal_left
                    && $cal_left < $cal_right
                    && (   $cal_left != $$ma_r{$k}{lend}
                        || $cal_right != $$ma_r{$k}{rend} )
                  )
                {
                    print $LOG1
"HEADER_ERROR\t$header\twrong coordinates in seq $k. Calculated: ($cal_left, $cal_right), in DB: ($$ma_r{$k}{lend}, $$ma_r{$k}{rend})\n"
                      if ($LOG1);
                    print $LOG2
"HEADER_ERROR\t$header\twrong coordinates in seq $k. Calculated: ($cal_left, $cal_right), in DB: ($$ma_r{$k}{lend}, $$ma_r{$k}{rend})\n"
                      if ($LOG2);
                    print $LOG1
"HEADER_UPDATE\t$header\tthe old lend and rend will be replaced by new ones in MA $k\n"
                      if ($LOG1);
                    print $LOG2
"HEADER_UPDATE\t$header\tthe old lend and rend will be replaced by new ones in MA $k\n"
                      if ($LOG2);
                    $$ma_r{$k}{lend} = $cal_left;
                    $$ma_r{$k}{rend} = $cal_right;
                    $$flag_r .=
"the old lend and rend will be replaced by new ones in MA $k\n"
                      if ($flag_r);
                }
            }
            else {
                print $LOG1
"SEQ_MATCH_ERROR\t$header\tcan not fit seq $k in MA with that from protein table linked by the header of seq\n"
                  if ($LOG1);
                print $LOG2
"SEQ_MATCH_ERROR\t$header\tcan not fit seq $k in MA with that from protein table linked by the header of seq\nseq_MA: $a\nprot_seq: $b\nold header: $$ma_r{$k}{header}\n"
                  if ($LOG2);
                my $old_header = $$ma_r{$k}{header};
                $$ma_r{$k}{lend} = $$ma_r{$k}{rend} = 0;
                my $success =
                  &search_nraa( $$ma_r{$k}, $flag_r, $fatal_r, $header, $LOG1,
                    $LOG2 )
                  if ($blastp);
                if ($success) {
                    print $LOG1
"               \t$header\t$old_header\t$$ma_r{$k}{header_new}\tsearch against nraa is successful\n"
                      if ($LOG1);
                    print $LOG2
"               \t$header\t$old_header\t$$ma_r{$k}{header_new}\tsearch against nraa is successful\n"
                      if ($LOG2);
                }
                else {
                    print $LOG1
"SEQ_MATCH_ERROR\t$header\tcan not find the seq $k from nraa\n"
                      if ($LOG1);
                    print $LOG2
"SEQ_MATCH_ERROR\t$header\tcan not find the seq $k from nraa\n$b\n"
                      if ($LOG2);
                    $$fatal_r .= "can not find the seq $k from nraa\n"
                      if ($fatal_r);
                    next;
                }
            }
        }
    }
}

sub run_blastp {
    use strict;
    my (
        $blastp_cmd,    $seq_r,            $seq_num,
        $tmp_root,      $blastp_db_file,   $dir,
        $blastp_option, $blastp_logfile_r, $log_header,
        $log_r
    ) = @_;
    $log_header .= "--sub run_blastp";
    if ( $seq_num !~ /^\d+$/ ) {
        $$log_r .=
"***Error: $log_header--an individual sequence in the multiple sequences must be specified\n"
          if ($log_r);
        return 1;
    }
    elsif ($seq_num < 0
        or $seq_num >= $$seq_r{number} )
    {
        $$log_r .=
"***Error: $log_header--the number of the sequence specified ($seq_num) is out of range(0 to ${\($$seq_r{number}-1)})\n"
          if ($log_r);
        return 1;
    }
    my $line = &format_sequence( $seq_r, 'fasta', $seq_num, 60, 0, '', '', '' );
    $tmp_root .= '_';
    my $tmp_file = &tmp_file( $dir, $tmp_root );
    open( TMP, ">$dir/$tmp_file.fasta" );
    print TMP $line;
    close TMP;
    system(
"$blastp_cmd $blastp_db_file $dir/$tmp_file.fasta $blastp_option > $dir/$tmp_file.blast_log"
    );
    $$blastp_logfile_r = "$dir/$tmp_file.blast_log";

    #    unlink ("$dir/$tmp_file", "$dir/$tmp_file.fasta");
}

sub make_seed_from_blastp_hits {
    use strict;
    my (
        $seed_r,               $seq_name,       $blastp_db_file,
        $blastp_db_index_file, $blastp_logfile, $blastp_e_cutoff,
        $max_seed,             $log_header,     $log_r
    ) = @_;
    $log_header .= "--sub make_seed_from_blastp_hits";
    system("btab -q $blastp_logfile");
    open( FH, "$blastp_logfile.btab" );
    my %prot_acc;
    my $count = 0;
    while ( my $line = <FH> ) {
        ++$count;
        chomp $line;
        my @tmp = split /\t/, $line;
        next
          if ( $tmp[19] > $blastp_e_cutoff );
        last
          if ( $count > $max_seed );
        $prot_acc{ $tmp[5] } = $tmp[19];

        #	print "$tmp[5], $tmp[19]\n";
    }
    close FH;

    #cdbyank protein sequence
    if ( scalar keys %prot_acc == 0 ) {
        $$log_r .=
"***Error: $log_header--no sequence selected from blastp search for $seq_name, either there is not hit or the E value cutoff is too high\n"
          if ($log_r);
        return 0;
    }
    for my $prot_acc ( keys %prot_acc ) {
        open( FH,
"cdbyank $blastp_db_index_file -d $blastp_db_file -a \"$prot_acc\" | "
        );
        while (<FH>) {
            $$seed_r{$seq_name}{ori} .= $_;
        }
        close FH;
        print "***Error: $log_header--can not cdbyank protein $prot_acc\n"
          if (  $log_r
            and $$seed_r{$seq_name}{ori} !~ /\S/ );
        $$log_r .= "***Error: $log_header--can not cdbyank protein $prot_acc\n"
          if (  $log_r
            and $$seed_r{$seq_name}{ori} !~ /\S/ );
    }
    &parse_sequence( $$seed_r{$seq_name}, 'fasta', '', $log_r, $log_header, '',
        '' );

    #    unlink("$blastp_logfile.btab");
    return 0;
}

sub blastp_nraa {

    #will be obsolete, currently used by acc_conversion.pl
    use strict;
    my ( $seq_r, $seq_num, $flag_r, $fatal_r, $header, $LOG1, $LOG2 ) = @_;
    if ( $seq_num !~ /^\d+$/ ) {
        print $LOG1
"***Error: an individual sequence in the multiple sequences must be specified\n"
          if ($LOG1);
        print $LOG2
"***Error: an individual sequence in the multiple sequences must be specified\n"
          if ($LOG2);
        return;
    }
    elsif ($seq_num < 0
        or $seq_num >= $$seq_r{number} )
    {
        print $LOG1
"***Error: the number of the sequence specified ($seq_num) is out of range(0 to ${\($$seq_r{number}-1)})\n"
          if ($LOG1);
        print $LOG2
"***Error: the number of the sequence specified ($seq_num) is out of range(0 to ${\($$seq_r{number}-1)})\n"
          if ($LOG2);
        return;
    }
    my $line = &format_sequence( $seq_r, 'fasta', $seq_num, 60, 0, '', '', '' );
    my $tmp_file = &tmp_file( './', 'tmp' );
    open( TMP, ">$tmp_file.seq" );
    print TMP $line;
    close TMP;
    system(
'blastp /usr/local/db/panda/AllGroup/AllGroup.niaa $tmp_file.seq -matrix BLOSUM62 > $tmp_file.seq.blast'
    );

#   system('blastp /usr/local/db/panda/nraa/nraa $tmp_file.seq -matrix BLOSUM62 > $tmp_file.seq.blast');
    my ( $success, $percent, $header_new, $lend, $rend ) =
      &parse_blast( "$tmp_file.seq.blast", '' );
    print $LOG1
"blastp result: success = $success, percent = $percent, header = $header_new, $lend, $rend\n"
      if ($LOG1);
    print $LOG2
"blastp result: success = $success, percent = $percent, header = $header_new, $lend, $rend\n"
      if ($LOG2);
    unlink(
        "$tmp_file",           "$tmp_file.seq",
        "$tmp_file.seq.blast", "$tmp_file.seq.blast.btab"
    );

    if ($success) {
        print $LOG1
"HEADER_UPDATE\t$header\tthe header is update to $header_new that has $percent\% similarity\n"
          if ($LOG1);
        print $LOG2
"HEADER_UPDATE\t$header\tthe header is update to $header_new that has $percent\% similarity\n"
          if ($LOG2);
        $$flag_r .= "the whole header of the sequence is updated\n"
          if ($flag_r);
        my ( $t1, $t2 ) = split( /\s/, $header_new, 2 );
        $$seq_r{$seq_num}{ori} = $t1 . "\t" . $$seq_r{$seq_num}{seq_gap};
        for ( keys %{ $$seq_r{$seq_num} } ) {
            delete $$seq_r{$seq_num}{$_}
              if ( $_ ne 'ori' );
        }
        &parse_single_seq( $$seq_r{$seq_num},
            'mul', $flag_r, $fatal_r, $header, $LOG1, $LOG2 );
        $$seq_r{$seq_num}{lend} = $lend;
        $$seq_r{$seq_num}{rend} = $rend;
    }
    return ($success);
}

sub parse_blast {

    #done
    #$data[10] percentage of identities, $data[5] new header
    use strict;
    my ( $file, $save_dir ) = @_;
    my ( $success, @result, $lend, $rend, $percent, $header_new );
    $success = $percent = 0;
    (
        $save_dir ne ''
        ? system("btab -q -p $save_dir $file")
        : system("btab -q $file")
    );
    if ( open( TMP, "$file.btab" ) ) {
        while (<TMP>) {
            chomp $_;
            @result = split /\t/;
            if (   $result[10] >= 98
                && $result[6] <= 6
                && $result[7] >= $result[2] - 5 )
            {
                $success = 1;
                if (
                    $result[10] > $percent
                    or (   $result[10] == $percent
                        && $result[5] =~ /^EGAD/ )
                  )
                {
                    $percent    = $result[10];
                    $lend       = $result[8] - $result[6] + 1;
                    $rend       = $result[9] + $result[2] - $result[7];
                    $header_new = $result[5];
                }
            }
        }
        close TMP;
    }
    return ( $success, $percent, $header_new, $lend, $rend );
}

sub check_MA_prot_id {

# find and add prot_id to protein list if it is in MA but not seen in protein list yet
#done
    use strict;
    my (
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    ) = @_;
    foreach my $key ( 0 .. $$ma_r{number} - 1 ) {
        my $prot_id = $$ma_r{$key}{prot_id};
        if ( $prot_id ne '' ) {
            if ( $$prot_r{$prot_id} ) {
                next;
            }
            else {
                my $query =
                  "select prot_seq from egad..protein where prot_id = $prot_id";
                my $array_r = &do_query( $db_proc, $query );
                if ( $$array_r[0]{prot_seq} ne '' ) {
                    print $LOG2
"HEADER_MESSAGE\t$header\tprot_id $prot_id is found in protein table for sequence $key\n";
                    $$prot_r{$prot_id}{prot_seq} = $$array_r[0]{prot_seq};
                }
                else {

                    #QUESTION: how to fix this, run blastp to change header?
                    print $LOG1
"HEADER_ERROR\t$header\tprot_id $prot_id is not found in protein table for sequence $key\n";
                    print $LOG2
"HEADER_ERROR\t$header\tprot_id $prot_id is not found in protein table for sequence $key\n";
                }
            }
        }
    }
}

sub check_locus_in_mul_format {

    #change '-' to '_' for mul file (hmm_seed file currently)
    #done
    use strict;
    my (
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    ) = @_;
    foreach my $k ( 0 .. $$ma_r{number} - 1 ) {
        if (   $$ma_r{$k}{db} eq 'EGAD'
            && $$ma_r{$k}{locus} =~ /-/ )
        {
            print $LOG1
"HEADER_ERROR\t$header\tdash is found in the locus name of sequence $k which is in mul format\n";
            print $LOG2
"HEADER_ERROR\t$header\tdash is found in the locus name of sequence $k which is in mul format\n";
            print $LOG1 "HEADER_UPDATE\t$header\tUpdating dash to underscore\n";
            print $LOG2 "HEADER_UPDATE\t$header\tUpdating dash to underscore\n";
            $$ma_r{$k}{locus} =~ s/-/_/g;
            $$flag_r .=
"Updating dash to underscore in locus name in mul format(currently the hmm_seed)\n";
        }
    }
}

sub check_MA_content {

    # check individual items in MA sequences
    # ma_r: hash reference of sequence alignment
    # type: not used here
    # prot_r: not used here
    # db_proc: not used here
    # blastp: not used here
    use strict;
    my (
        $ma_r,   $type,    $prot_r, $db_proc, $blastp,
        $flag_r, $fatal_r, $header, $LOG1,    $LOG2
    ) = @_;
    if ( $$ma_r{number} == 0 ) {
        $$fatal_r .= "no sequence found\n"
          if ($fatal_r);
        return;
    }
    for my $k ( 0 .. $$ma_r{number} - 1 ) {

        # check prot_id
        if (   $$ma_r{$k}{db} eq 'EGAD'
            && $$ma_r{$k}{prot_id} =~ /\D/ )
        {
            print $LOG1
"HEADER_ERROR\t$header\twrong characters in prot_id for sequence $k $$ma_r{$k}{prot_id}\n"
              if ($LOG1);
            print $LOG2
"HEADER_ERROR\t$header\twrong characters in prot_id for sequence $k $$ma_r{$k}{prot_id}\n"
              if ($LOG2);
        }

        # check locus name
        if (   $$ma_r{$k}{db} eq 'EGAD'
            && $$ma_r{$k}{locus} =~ /[^a-zA-Z0-9._\-]/ )
        {
            print $LOG1
"HEADER_ERROR\t$header\twrong characters in locus name for sequence $k, $$ma_r{$k}{locus}\n"
              if ($LOG1);
            print $LOG2
"HEADER_ERROR\t$header\twrong characters in locus name for sequence $k, $$ma_r{$k}{locus}\n"
              if ($LOG2);
        }

        # check lend
        if ( defined $$ma_r{$k}{lend}
            && $$ma_r{$k}{lend} =~ /\D/ )
        {
            print $LOG1
"HEADER_ERROR\t$header\twrong characters in lend for sequence $k, $$ma_r{$k}{lend}\n"
              if ($LOG1);
            print $LOG2
"HEADER_ERROR\t$header\twrong characters in lend for sequence $k, $$ma_r{$k}{lend}\n"
              if ($LOG2);
        }

        # check rend
        if ( defined $$ma_r{$k}{rend}
            && $$ma_r{$k}{rend} =~ /\D/ )
        {
            print $LOG1
"HEADER_ERROR\t$header\twrong characters in rend for sequence $k, $$ma_r{$k}{rend}\n"
              if ($LOG1);
            print $LOG2
"HEADER_ERROR\t$header\twrong characters in rend for sequence $k, $$ma_r{$k}{rend}\n"
              if ($LOG2);
        }

        # check space in sequence, change dash to dot
        if ( $$ma_r{$k}{seq_gap} =~ /-/ ) {
            print $LOG1
"HEADER_ERROR\t$header\tdash is found in sequence $k representing gap\n"
              if ($LOG1);
            print $LOG2
              "HEADER_UPDATE\t$header\tsequence $k by replacing dash with dot\n"
              if ($LOG2);
            $$ma_r{$k}{seq_gap} =~ s/-/./g;
            $$flag_r .= "Updating sequence $k by replacing dash with dot\n"
              if ($flag_r);
        }

        # check sequence without gap
        if ( $$ma_r{$k}{seq} eq '' ) {
            print $LOG1 "SEQ_ERROR\t$header\tempty sequence \#$k (no gap)\n"
              if ($LOG1);
            print $LOG2 "SEQ_ERROR\t$header\tempty sequence \#$k (no gap)\n"
              if ($LOG2);
            $$fatal_r .= "empty sequence \#$k (no gap)\n"
              if ($fatal_r);
        }
        if ( ( uc $$ma_r{$k}{seq} ) =~ /([^A-IK-NP-Z])/ )
        {    # just changed to allow U for selenocysteine. DHH.
            print $LOG1
"SEQ_ERROR\t$header\twrong characters($1) in sequence \#$k (no gap)\n"
              if ($LOG1);
            print $LOG2
"SEQ_ERROR\t$header\twrong characters($1) in sequence \#$k (no gap), $$ma_r{$k}{seq}\n"
              if ($LOG2);
            $$fatal_r .= "wrong characters($1) in sequence \#$k (no gap)\n"
              if ($fatal_r);
        }

        # check sequence with gap
        if ( $$ma_r{$k}{seq_gap} eq '' ) {
            print $LOG1 "SEQ_ERROR\t$header\tempty sequence \#$k (with gap)\n"
              if ($LOG1);
            print $LOG2 "SEQ_ERROR\t$header\tempty sequence \#$k (with gap)\n"
              if ($LOG2);
            $$fatal_r .= "empty sequence \#$k (with gap)\n"
              if ($fatal_r);
        }
        if ( ( uc $$ma_r{$k}{seq_gap} ) =~ /([^A-IK-NP-Z\-.])/ )
        {    # just changed to allow U for selenocysteine. DHH.
            print $LOG1
"SEQ_ERROR\t$header\twrong characters($1) in sequence \#$k (with gap)\n"
              if ($LOG1);
            print $LOG2
"SEQ_ERROR\t$header\twrong characters($1) in sequence \#$k (with gap), $$ma_r{$k}{seq_gap}\n"
              if ($LOG2);
            $$fatal_r .= "wrong characters($1) in sequence \#$k (with gap)\n"
              if ($fatal_r);
        }
        for ( 0 .. $$ma_r{$k}{number} - 1 ) {
            if ( $$ma_r{$k}{$_} =~ /[^\w\-.\[\]]/ ) {
                print $LOG1
"HEADER_ERROR\t$header\twrong header characters in sequence $k\n"
                  if ($LOG1);
                print $LOG2
"HEADER_ERROR\t$header\twrong header characters in sequence $k, $$ma_r{$k}{header}\n"
                  if ($LOG2);
            }
        }

#QUESTION: what to do if unknown DB name is found? implement blastp to get header?
        if (   defined $$ma_r{$k}{0}
            && $$ma_r{$k}{0} ne 'EGAD'
            && $$ma_r{$k}{0} ne 'SP'
            && $$ma_r{$k}{0} ne 'GP'
            && $$ma_r{$k}{0} ne 'PIR'
            && $$ma_r{$k}{0} ne 'gi'
            && $$ma_r{$k}{0} ne 'SP_ID'
            && $$ma_r{$k}{0} ne 'OMNI' )
        {
            print $LOG1
"HEADER_ERROR\t$header\twrong header with unknown DB name $$ma_r{$k}{0}\n"
              if ($LOG1);
            print $LOG2
"HEADER_ERROR\t$header\twrong header with unknown DB name $$ma_r{$k}{0}\n"
              if ($LOG2);
        }

        # each gapped sequence in MA should have the same length.
        my $s1 = length( $$ma_r{0}{seq_gap} );
        my $s2 = length( $$ma_r{$k}{seq_gap} );
        if ( $s1 != $s2 ) {
            print $LOG1
"ALIGNMENT_ERROR\t$header\t$$ma_r{$k}{header}\tthe length($s2) of sequence $k is not equal to that($s1) of sequence 0\n"
              if ($LOG1);
            print $LOG2
"ALIGNMENT_ERROR\t$header\t$$ma_r{$k}{header}\tthe length($s2) of sequence $k is not equal to that($s1) of sequence 0\n"
              if ($LOG2);
            $$fatal_r .=
"BAD_ALIGNMENT\t$header\t$$ma_r{$k}{header}\tthe length($s2) of sequence $k is not equal to that($s1) of sequence 0\n"
              if ($fatal_r);
        }
    }
}
#############################################################################################################################################
sub process_data {
    my ( $info_r, $data_r, $flag_r, $fatal_r, $LOG1, $LOG2 ) = @_;
    &parse_sequence( \%{ $data_r->{hmm_seed} },
        "mul", \$$flag_r{hmm_seed}, $fatal_r, $$data_r{header}, $LOG1, $LOG2 )
      if ( $$data_r{hmm_seed}{ori} );
    &parse_sequence( \%{ $data_r->{align_seed} },
        "fasta", \$$flag_r{align_seed}, $fatal_r, $$data_r{header}, $LOG1,
        $LOG2 )
      if ( $$data_r{align_seed}{ori} );
    &parse_sequence( \%{ $data_r->{align_full} },
        "fasta", \$$flag_r{align_full}, $fatal_r, $$data_r{header}, $LOG1,
        $LOG2 )
      if ( $$data_r{align_full}{ori} );
}

sub parse_sequence {

# parse multiple sequences
# seq_r: hash reference for multiple sequences
# format: the format of sequences
# flag_r: flag for updating database. can be null.
# fatal_r: fatal error flag for updating database. can be null.
# header: header information for log files
# LOG1: log file for short output
# LOG2: log file for long output
# Input a scalar variable of multiple sequence in mul, fasta, or other format. Split it by
# "\n", ">", or call subroutines, based on format.  Pass each fragment to parse_single_seq.
    use strict;
    my ( $seq_r, $format, $flag_r, $fatal_r, $header, $LOG1, $LOG2 ) = @_;
    my @lines;
    if ( $format eq "mul" ) {
        my @li                    = split( /\n/, $$seq_r{ori} );
        my $count_hmmbuild_method = -1;
        my $count_mini_db         = -1;
        for my $n (@li) {
            if ( $li[0] =~ /^\#\s*STOCKHOLM/ ) {
                if ( $n =~ /^\#\s*STOCKHOLM/ ) {
                }
                elsif ( $n =~ /^\#=GF\s+hmmbuild_acc\s+(\S+)$/ ) {
                    $$seq_r{stockholm}{GF}{hmmbuild_acc} = $1;
                }
                elsif ( $n =~ /^\#=GF\s+hmmbuild_method\s+(.+)$/ ) {
                    ++$count_hmmbuild_method;
                    $$seq_r{stockholm}{GF}{hmmbuild_method}
                      {$count_hmmbuild_method} = $1;
                }
                elsif ( $n =~ /^\#=GF\s+project_id\s+(\d+)$/ ) {
                    $$seq_r{stockholm}{GF}{project_id} = $1;
                }
                elsif ( $n =~ /^\#=GF\s+iso_id\s+(\d+)$/ ) {
                    $$seq_r{stockholm}{GF}{iso_id} = $1;
                }
                elsif ( $n =~
/^\#=GF\s+MINI_DB\s+ACC\s+(\S+)\s+SEQ\s+(\S+)\s+CMD\s+\"([^\"]+)\"\s+DB\s+(\S+)\s+(\S+)\s+\"(.+)\"\s*$/
                  )
                {
                    ++$count_mini_db;
                    $$seq_r{stockholm}{GF}{mini_db}{$count_mini_db}{acc} = $1;
                    $$seq_r{stockholm}{GF}{mini_db}{$count_mini_db}{seq} = $2;
                    $$seq_r{stockholm}{GF}{mini_db}{$count_mini_db}{cmd} = $3;
                    $$seq_r{stockholm}{GF}{mini_db}{$count_mini_db}{db_name} =
                      $4;
                    $$seq_r{stockholm}{GF}{mini_db}{$count_mini_db}{db_size} =
                      $5;
                    $$seq_r{stockholm}{GF}{mini_db}{$count_mini_db}{db_date} =
                      $6;
                }
                elsif ( $n =~
/^\#=GF\s+SCORE\s+PROT\s+(\S+)\s+HMM\s+(\S+)\s+B_GLOBAL\s+(\S+)\s+E_GLOBAL\s+(\S+)\s+B_FRAG\s+(\S+)\s+E_FRAG\s+(\S+)\s+TOTAL_FRAG\s+(\S+)\s+NUM_FRAG\s(\S+)\s*$/
                  )
                {
                    $$seq_r{stockholm}{GF}{score}{$1}{hmm}        = $2;
                    $$seq_r{stockholm}{GF}{score}{$1}{b}          = $3;
                    $$seq_r{stockholm}{GF}{score}{$1}{e}          = $4;
                    $$seq_r{stockholm}{GF}{score}{$1}{total_frag} = $7;
                    if ( $$seq_r{stockholm}{GF}{score}{$1}{total_frag} > 1 ) {
                        $$seq_r{stockholm}{GF}{score}{$1}{$8}{b} = $5;
                        $$seq_r{stockholm}{GF}{score}{$1}{$8}{e} = $6;
                    }
                }
                elsif ( $n =~ /^\#=GF\s+INFO\s+(.+)$/ ) {
                    my @tmp = split /<>/, $1;
                    my %hash;
                    for my $m (@tmp) {
                        $hash{$1} = $2
                          if ( $m =~ /^(\S+)\s+(.+)$/ );
                    }
                    $$seq_r{stockholm}{GF}{info}{ $hash{PROT} }{genus} =
                      $hash{GENUS};
                    $$seq_r{stockholm}{GF}{info}{ $hash{PROT} }{species} =
                      $hash{SPECIES};
                    $$seq_r{stockholm}{GF}{info}{ $hash{PROT} }{length} =
                      $hash{LENGTH};
                    $$seq_r{stockholm}{GF}{info}{ $hash{PROT} }{name} =
                      $hash{NAME};
                    $$seq_r{stockholm}{GF}{info}{ $hash{PROT} }{first_acc} =
                      $hash{FIRST_ACC};
                }
                elsif ( $n =~ /^\#=/ ) {
                }
                else {
                    push @lines, $n;
                }
            }
            else {
                push @lines, $n;
            }
        }
        $$seq_r{stockholm}{GF}{hmmbuild_method}{number} =
          $count_hmmbuild_method + 1
          if ( $li[0] =~ /^\#\s*STOCKHOLM/ );
        $$seq_r{stockholm}{GF}{mini_db}{number} = $count_mini_db + 1
          if ( $li[0] =~ /^\#\s*STOCKHOLM/
            && $count_mini_db > -1 );
    }
    elsif ( $format eq "fasta" ) {
        @lines = split( />/, $$seq_r{ori} );
        shift(@lines)
          if ( $lines[0] eq '' );
    }
    elsif ( $format eq 'selex' ) {
        my $num = -1;
        my %link;
        @lines = split( /\n/, $$seq_r{ori} );
        for my $n (@lines) {
            if ( $n =~ /^\#=SQ\s+(\S+)/ ) {
                ++$num;
                $$seq_r{$num}{header} = $1;
                $link{$1} = $num;
                if ( $$seq_r{$num}{header} =~ /^(\S+?)(\/(\d+)-(\d+))?$/ ) {
                    $$seq_r{$num}{lend} = $3
                      if ( $3 ne '' );
                    $$seq_r{$num}{rend} = $4
                      if ( $4 ne '' );
                    my @tmp = split( /\|/, $1 );
                    if (   @tmp == 3
                        && $tmp[0] eq 'EGAD' )
                    {
                        $$seq_r{$num}{db}      = 'EGAD';
                        $$seq_r{$num}{prot_id} = $tmp[1];
                        $$seq_r{$num}{locus}   = $tmp[2];
                    }
                    else {
                        $$seq_r{$num}{db} = $tmp[0];
                        for my $m ( 1 .. @tmp - 1 ) {
                            $$seq_r{$num}{ $m - 1 } = $tmp[$m];
                        }
                        $$seq_r{$num}{number} = @tmp - 1;
                    }
                }
            }
            elsif ( $n =~ /^\#=RF\s+(\S+)/ ) {
                $$seq_r{RF} .= $1;
            }
        }
        $$seq_r{RF} =~ s/^\s+|\s+$//g;
        $$seq_r{number} = $num + 1;
        for my $n (@lines) {
            $$seq_r{ $link{$1} }{seq_gap} .= $2
              if (  $n !~ /^\#=/
                and $n =~ /^(\S+)\s+(\S+)/ );
        }
        for my $n ( 0 .. $$seq_r{number} - 1 ) {
            $$seq_r{$n}{seq} = $$seq_r{$n}{seq_gap};
            $$seq_r{$n}{seq} =~ s/[.\-]//g;
        }
    }
    elsif ( $format eq 'msf' ) {
        my ( $num, $n ) = ( -1, 0 );
        @lines = split( /\n/, $$seq_r{ori} );
        while ( $n < @lines ) {
            if ( $lines[$n] =~ /^\s+Name:\s+(\S+).*\s+Len:/ ) {
                ++$num;
                $$seq_r{$num}{header} = $1;
                if ( $$seq_r{$num}{header} =~ /^(\S+?)(\/(\d+)-(\d+))?$/ ) {
                    $$seq_r{$num}{lend} = $3
                      if ( $3 ne '' );
                    $$seq_r{$num}{rend} = $4
                      if ( $4 ne '' );
                    my @tmp = split( /\|/, $1 );
                    if (   @tmp == 3
                        && $tmp[0] eq 'EGAD' )
                    {
                        $$seq_r{$num}{db}      = 'EGAD';
                        $$seq_r{$num}{prot_id} = $tmp[1];
                        $$seq_r{$num}{locus}   = $tmp[2];
                    }
                    else {
                        $$seq_r{$num}{db} = $tmp[0];
                        for my $m ( 1 .. @tmp - 1 ) {
                            $$seq_r{$num}{ $m - 1 } = $tmp[$m];
                        }
                        $$seq_r{$num}{number} = @tmp - 1;
                    }
                }
            }
            elsif ( $lines[$n] =~ /\/\// ) {
                $$seq_r{number} = $num + 1;
            }
            elsif ($$seq_r{0}{header} ne ''
                && $lines[$n] =~ /^\s*\Q$$seq_r{'0'}{header}\E/ )
            {
                $$seq_r{0}{seq_gap} .= $';
                for my $a ( 1 .. $$seq_r{number} - 1 ) {
                    ++$n;
                    if ( $lines[$n] =~ /^\s*\Q$$seq_r{$a}{header}\E/ ) {
                        $$seq_r{$a}{seq_gap} .= $';
                    }
                    else {
                        print $LOG1
                          "SEQ_ERROR\t$header\twrong sequence in the msf file\n"
                          if ($LOG1);
                        print $LOG2
                          "SEQ_ERROR\t$header\twrong sequence in the msf file\n"
                          if ($LOG2);
                    }
                }
            }
            ++$n;
        }
        for ( 0 .. $$seq_r{number} - 1 ) {
            $$seq_r{$_}{seq_gap} =~ s/[\s\n\t]//g;
            $$seq_r{$_}{seq} = $$seq_r{$_}{seq_gap};
            $$seq_r{$_}{seq} =~ s/[.\-]//g;
        }
    }
    elsif ( $format eq 'pfam' ) {

        # the information get here is used to find EGAD prot_id later
        my $num     = 0;
        my $seq_num = -1;
        @lines = split( /\n/, $$seq_r{ori} );
        while ( $num < @lines ) {
            chomp $lines[$num];
            if ( $lines[$num] !~ /^\#=/ ) {
                ++$seq_num;
                ( $$seq_r{$seq_num}{header}, $$seq_r{$seq_num}{seq_gap} ) =
                  split( /[ \t]+/, $lines[$num], 2 );
                if ( $$seq_r{$seq_num}{header} =~ /^(\S+?)(\/(\d+)-(\d+))?$/ ) {

                    # QUESTION: Is following a locus????
                    $$seq_r{$seq_num}{0}      = $1;
                    $$seq_r{$seq_num}{number} = 1;
                    $$seq_r{$seq_num}{lend}   = $3
                      if ( $3 ne '' );
                    $$seq_r{$seq_num}{rend} = $4
                      if ( $4 ne '' );
                }
                while ( $lines[ ++$num ] =~ /^\#=/ ) {
                    chomp $lines[$num];
                    if ( $lines[$num] =~ /^\#=GS \S+\s+AC (\S+)$/ ) {
                        $$seq_r{$seq_num}{pfam_acc} = $1;
                    }
                }
                print $LOG1 "***Error: pfam_acc not found\n"
                  if ( $LOG1
                    && $$seq_r{$seq_num}{pfam_acc} eq '' );
                print $LOG2 "***Error: pfam_acc not found\n"
                  if ( $LOG2
                    && $$seq_r{$seq_num}{pfam_acc} eq '' );
            }
            else {
                ++$num;
            }
        }
        $$seq_r{number} = $seq_num + 1;
        for ( 0 .. $$seq_r{number} - 1 ) {
            $$seq_r{$_}{seq_gap} =~ s/[\s\n\t]//g;
            $$seq_r{$_}{seq} = $$seq_r{$_}{seq_gap};
            $$seq_r{$_}{seq} =~ s/[.\-]//g;

            #	print "$$ref->{$_}{header}\t$$ref->{$_}{seq}\n";
        }
    }
    else {
        print $LOG1 "SEQ_FORMAT_ERROR\t$header\tunknown format: $format\n"
          if ($LOG1);
        print $LOG2 "SEQ_FORMAT_ERROR\t$header\tunknown format: $format\n"
          if ($LOG2);
    }
    if (   $format eq 'mul'
        or $format eq 'fasta' )
    {
        my $num = -1;
        foreach (@lines) {
            if ( $_ =~ /\S/ ) {
                ++$num;
                $$seq_r{$num}{ori} = $_;
                &parse_single_seq( $$seq_r{$num}, $format, $flag_r, $fatal_r,
                    $header, $LOG1, $LOG2 );
            }
            else {
                print $LOG1
"$header\tSEQUENCE_ERR\tempty entry in multiple sequence alignment\n"
                  if ($LOG1);
                print $LOG2
"$header\tSEQUENCE_ERR\tempty entry in multiple sequence alignment\n"
                  if ($LOG2);
                $$flag_r .= "remove one empty line from MA\n"
                  if ($flag_r);
            }
        }
        $$seq_r{number} = $num + 1;
    }
}

sub parse_single_seq {

# the header is decomposed to db, prot_id, locus, lend, rend, comment for EGAD
# the header is decomposed to db, 0,1,..,n,number, lend, rend, comment for known DB name.
# the header is decomposed to 0,...n,number, lend, rend, comment for all other cases
# so the general form is (db,prot_id, locus, 0, ..., n,number, lend, rend, comment) if each field is not null.
    use strict;
    my ( $ref, $format, $flag_r, $fatal_r, $header_log, $LOG1, $LOG2 ) = @_;
    return
      if ( !defined $$ref{ori}
        or $$ref{ori} eq '' );
    my $input = $$ref{ori};
    if ( $input =~ /(\A +| +\Z)/ ) {
        $$flag_r .=
"$header_log\tHEADER_ERROR\tremove heading and tailing spaces in header: $input\n"
          if ($flag_r);
        $input =~ s/(\A +| +\Z)//g;    #remove heading and tailing spaces
    }

    # deal with a special error: 4567 EGAD|4567|XXXXX
    if ( $input =~ /^\d+ *(EGAD|SP)/ ) {
        print $LOG1
"$header_log\tHEADER_ERROR\tthe substring of the wrong header is: ${\(substr($input, 0, 20))}\n"
          if ($LOG1);
        print $LOG2
"$header_log\tHEADER_ERROR\tthe substring of the wrong header is: ${\(substr($input, 0, 20))}\n"
          if ($LOG2);
        $$flag_r .=
          "the accession header is changed from ${\(substr($input, 0, 15))} to "
          if ($flag_r);
        $input =~ s/^\d+ *//;
        print $LOG1
"$header_log ******Updating: the substring of the new header is: ${\(substr($input, 0, 20))}\n"
          if ($LOG1);
        print $LOG2
"$header_log ******Updating: the substring of the new header is: ${\(substr($input, 0, 20))}\n"
          if ($LOG2);
        $$flag_r .= "${\(substr($input, 0, 15))}\n"
          if ($flag_r);
    }
    if ( $input =~ /^\|/ ) {
        print $LOG2
"HEADER_ERROR\t$header_log\textra '|' at the beginning of the header\n"
          if ($LOG2);
        $$flag_r .= "remove extra '|' at the beginning of the header\n"
          if ($flag_r);
        $input =~ s/^\|//;
        print $LOG1
"HEADER_UPDATE\t$header_log\tremove extra '|' at the beginning of the header\n"
          if ($LOG1);
        print $LOG2
"HEADER_UPDATE\t$header_log\tremove extra '|' at the beginning of the header\n"
          if ($LOG2);
    }

    # for mul or fasta format, separate the entry to header and sequence
    my ( $header, @header );
    if ( $format eq 'mul' ) {
        ( $header, $$ref{'seq_gap'} ) = split( /[ \t]+/, $input, 2 );
        $header =~ s/^\s+|\s+$//g;
    }
    elsif ( $format eq 'fasta' ) {
        ( $header, $$ref{'seq_gap'} ) = split( /\n/, $input, 2 );
        $header =~ s/^\s+|\s+$//g;

   # check sequence length of each line in fasta file which should be 60 or less
        my @tmp = split( /\n/, $$ref{seq_gap} );
        foreach (@tmp) {
            if ( length($_) > 60 ) {
                print $LOG2
"$header_log ***Error: the length of lines in fasta file is longer than 60 characters\n"
                  if ($LOG2);
                print $LOG2
"$header_log ******Updating: the length of lines in fasta file will be adjusted to 60 characters or less\n"
                  if ($LOG2);
                $$flag_r .=
"the length of lines in fasta file will be adjusted to 60 characters or less\n"
                  if ($flag_r);
                last;
            }
        }
        $$ref{'seq_gap'} =~ s/\n//g;
    }
    else {
        print $LOG1 "$header_log ***Error: unknown file type\n"
          if ($LOG1);
        print $LOG2 "$header_log ***Error: unknown file type\n"
          if ($LOG2);
        die "$header_log ***Error: unknown file type\n";
    }
    $$ref{header} = $header;
    $$ref{'seq'} = $$ref{'seq_gap'};
    $$ref{'seq'} =~ s/[-.]//g;

#suppose the header has form like: GP|XXX|XXX/29-50 XX|XX|XXX   or  GP|XXX|XXX XX|XX|XX
    my $tmp1;
    if ( $header =~ /^(\S+)\s+(\S+.+)$/ ) {
        $$ref{comment} = $2;
        $tmp1 = $1;
        @header = split( /[|\/]/, $tmp1 );
    }
    else {
        @header = split( /[|\/]/, $header );
    }
    for (@header) {
        print $LOG1 "HEADER_ERROR\t$header_log\textra '|' or '/' in header\n"
          if ( $_ eq ''
            && $header[0] ne 'GP'
            && $LOG1 );
        print $LOG2
          "HEADER_ERROR\t$header_log\textra '|' or '/' in header: $tmp1\n"
          if ( $_ eq ''
            && $header[0] ne 'GP'
            && $LOG2 );
    }
    my $num = @header;

    #find lend, rend and comment if exist, and save all other in {0,1,2...}
    if ( $header[ $num - 1 ] =~ /^(\d+)-(\d+)$/ ) {
        ( $$ref{lend}, $$ref{rend} ) = ( $1, $2 );
        --$num;
    }
    $$ref{db} = $header[0];
    if ( $$ref{db} eq 'EGAD' ) {
        if ( $num != 3 ) {
            print $LOG1 "$header_log\tHEADER_ERROR\tstrange header $header\n"
              if ($LOG1);
            print $LOG2 "$header_log\tHEADER_ERROR\tstrange header $header\n"
              if ($LOG2);
        }
        $$ref{prot_id} = $header[1];
        $$ref{locus}   = $header[2]
          if ( $header[2] !~ /^(\d+)-(\d+)$/ );
        for ( 3 .. $num - 1 ) {
            $$ref{ $_ - 3 } = $header[$_];
        }
        $num -= 2;
    }
    else {
        for ( 1 .. $num - 1 ) {
            $$ref{ $_ - 1 } = $header[$_];
        }
    }
    $$ref{number} = $num - 1;
}
#############################################################################################################################################
sub get_data_from_db {
    use strict;
    my ( $info_r, $data_r, $flag_r, $fatal_r, $LOG1, $LOG2 ) = @_;

# get_hmm_seed_from_db must be run first to get iso_id and hmm_acc and compose $header.
    print $LOG2 "get hmm_seed from db.............................\n";
    &get_hmm_seed_from_db( $info_r, $data_r, $flag_r, $fatal_r, $LOG1, $LOG2 );
    print $LOG2 "get alignment_seed from db.......................\n";
    &get_alignment_from_db( $info_r, $data_r, "seed", $flag_r, $fatal_r, $LOG1,
        $LOG2 );
    print $LOG2 "get alignment_full from db.......................\n";
    &get_alignment_from_db( $info_r, $data_r, "all", $flag_r, $fatal_r, $LOG1,
        $LOG2 );
    print $LOG2 "get iso_link from db.....................\n";
    &get_iso_link_from_db( $info_r, $data_r, $LOG1, $LOG2 );
    print $LOG2 "get full length protein sequence from db..................\n";
    &get_protein_from_db( $$info_r{db}, $$info_r{db_proc}, "iso_link", "iso_id",
        $data_r, $LOG1, $LOG2 );
}

sub get_protein_from_db {

    #modified on 3-31-2000, not tested yet.
    use strict;
    my ( $db, $db_proc, $link_table, $link_id, $data_r, $LOG1, $LOG2 ) = @_;
    my $array_r = &do_query( $db_proc,
"select p.prot_id, p.prot_seq from $db..$link_table i, $db..protein p where i.$link_id = $$data_r{$link_id} and i.prot_id = p.prot_id order by p.prot_id"
    );
    for (@$array_r) {
        $$data_r{protein}{ $$_{prot_id} }{prot_seq} = $$_{prot_seq};
    }
}

sub get_iso_link_from_db {
    use strict;
    my ( $info_r, $data_r, $LOG1, $LOG2 ) = @_;
    my $query =
"select * from egad..iso_link where iso_id = $$data_r{iso_id} order by prot_id";
    $$data_r{iso_link} = &do_query( $$info_r{db_proc}, $query );
}

sub get_alignment_from_db {

    # $$info_r{db}, $$info_r{db_proc} must be defined
    # $$data_r{iso_id} must be defined.
    use strict;
    my ( $info_r, $data_r, $type, $flag_r, $fatal_r, $LOG1, $LOG2 ) = @_;
    if (   $type ne 'seed'
        && $type ne 'all' )
    {
        print $LOG1
"$$info_r{header}\tPROGRAMMING_ERROR\twrong align_type to search. current types are 'all' and 'seed'\n";
        print $LOG2
"$$info_r{header}\tPROGRAMMING_ERROR\twrong align_type to search. current types are 'all' and 'seed'\n";
    }
    my $query =
"select * from $$info_r{db}..alignment where iso_id = $$data_r{iso_id} and align_type = \'$type\'";
    my @result;
    my $array_r = &do_query( $$info_r{db_proc}, $query );
    if ( $type eq 'seed' ) {
        $$data_r{align_id_seed} = $$array_r[0]{align_id};
        $$data_r{align_seed}{ori} = $$array_r[0]{alignment};
    }
    elsif ( $type eq 'all' ) {
        $$data_r{align_id_full} = $$array_r[0]{align_id};
        $$data_r{align_full}{ori} = $$array_r[0]{alignment};
    }
}

sub get_hmm_seed_from_db {
    use strict;
    my ( $info_r, $data_r, $flag_r, $fatal_r, $LOG1, $LOG2 ) = @_;
    my $query = "select * from egad..hmm2 where id = $$data_r{hmm_id}";
    my @result;
    my $array_r = &do_query( $$info_r{db_proc}, $query );
    $$data_r{hmm_seed}{ori} = $$array_r[0]{hmm_seed};
    $$data_r{iso_id} = $$array_r[0]{iso_id};
    print "***Error: wrong hmm_acc, why?? "
      if ( $$data_r{hmm_acc} != $$array_r[0]{hmm_acc} );
    $$data_r{header} = "$$data_r{hmm_acc}\t$$data_r{iso_id}";
}

sub get_hmm_ids_from_db {

    # find id numbers for a range of TIGR acc number
    #take 4 variables: $start, $end, $db_proc, $DEBUG, and return an array
    # the query need to be modified if database changes
    use strict;
    my ( $family, $date, $start, $end, $db_proc ) = @_;
    my ( %id_hash, $query );
    if ( $family eq 'TIGR' ) {
        $query =
"select id, hmm_acc from egad..hmm2 where is_current = 1 and hmm_acc like \'TIGR%\' ";
    }
    elsif ( $family eq 'PF' ) {
        $query =
"select id, hmm_acc from egad..hmm2 where is_current = 1 and hmm_acc like \'PF%\' ";
    }
    elsif ( $family eq 'ALL' ) {
        $query = "select id, hmm_acc from egad..hmm2 where is_current = 1";
    }
    else {
        die
"***Error: unknown protein family class. known class names are 'ALL', 'TIGR' and 'PF'\n";
    }
    $query .= " and (mod_date > \'$date\' or entry_date > \'$date\')"
      if ( $date ne '' );
    my $array_r = &do_query( $db_proc, $query );
    foreach (@$array_r) {

        #following selection is based on inputing ranges in hmm_acc
        if ( $$_{hmm_acc} =~ /(\d+)/ ) {
            $id_hash{ $$_{id} } = $$_{hmm_acc}
              if ( $start <= $1
                && $1 <= $end );
        }
    }
    return %id_hash;
}
#############################################################################################################################################
sub read_build {
    use strict;
    my ( $data_r, $f, $LOG ) = @_;
    my $log_header = $f;
    open( TMP, "$f" );
    while (<TMP>) {
        $$data_r{hmm_build}{value} .= $_;
        chomp $_;
        if ( $_ =~ /^Number of sequences:\s+(\d+)\s*$/ ) {
            $$data_r{num_seed}{value} = $1;
        }
        elsif ( $_ =~ /^Prior strategy:\s+/ ) {
            if ( index( $', "PAM hack" ) >= 0 ) {
                $$data_r{prior}{value} = 'blossum62';
            }
            else {
                $$data_r{prior}{value} = $';
            }
        }
        elsif ( $_ =~ /^PAM prior weight:\s+/ ) {

            #QUESTION: is following correct?
            $$data_r{prior}{value} .= " 20"
              if ( ( $' eq '' || index( $', 'bits' ) >= 0 )
                && $$data_r{prior}{value} !~ /Dirichlet/ );
            $$data_r{prior}{value} .= " $'"
              if ( ( $' ne '' && index( $', 'bits' ) < 0 )
                && $$data_r{prior}{value} !~ /Dirichlet/ );
        }
        elsif ( $_ =~ /^Constructed a profile HMM \(length (\d+)\)\s*$/ ) {
            $$data_r{hmm_leng}{value} = $1;
        }
        elsif ( $_ =~ /^Average score:\s+(\d+\.\d+)\s+bits\s*$/ ) {
            $$data_r{avg_score}{value} = $1;
        }
        elsif ( $_ =~ /^Minimum score:\s+(\d+\.\d+)\s+bits\s*$/ ) {
            $$data_r{min_score}{value} = $1;
        }
        elsif ( $_ =~ /^Maximum score:\s+(\d+\.\d+)\s+bits\s*$/ ) {
            $$data_r{max_score}{value} = $1;
        }
        elsif ( $_ =~ /^Std\. deviation:\s+(\d+\.\d+)\s+bits\s*$/ ) {
            $$data_r{std_dev}{value} = $1;
        }
    }
    close TMP;
}

sub read_hmm {
    use strict;
    my ( $data_r, $f ) = @_;
    open( TMP, "$f" );
    while (<TMP>) {
        $$data_r .= $_;
    }
    close TMP;
}

sub read_sequence {

    # $hash_r: hash reference for sequence information
    # $file: file name
    # $format: either input or determine automatically.
    use strict;
    my ( $hash_r, $file, $format, $LOG1, $LOG2 ) = @_;
    if ( -s $file ) {
        open( TMP, "$file" );
        while (<TMP>) {
            $$hash_r{ori} .= $_;
        }
        close TMP;
        if ( $format eq '' ) {
            $format = &auto_detect_format( $$hash_r{ori}, $LOG1 );
            print $LOG1 "file format of $file is determined to be $format\n"
              if ($LOG1);
            print $LOG2 "file format of $file is determined to be $format\n"
              if ($LOG2);
        }
        &parse_sequence( $hash_r, $format, '', '', '', $LOG1, $LOG2 );
    }
    else {
        print $LOG1 "***Error: can not find file: $file\n"
          if ($LOG1);
        print $LOG2 "***Error: can not find file: $file\n"
          if ($LOG2);
    }
    $$hash_r{number} = 0
      if ( !defined $$hash_r{number}
        or $$hash_r{number} eq '' );
    $$hash_r{stockholm}{GF}{mini_db}{number} = 0
      if ( !defined $$hash_r{stockholm}{GF}{mini_db}{number}
        or $$hash_r{stockholm}{GF}{mini_db}{number} eq '' );
    $$hash_r{stockholm}{GF}{hmmbuild_method}{number} = 0
      if ( !defined $$hash_r{stockholm}{GF}{hmmbuild_method}{number}
        or $$hash_r{stockholm}{GF}{hmmbuild_method}{number} eq '' );
}

sub auto_detect_format {
    use strict;
    my ( $file, $LOG1 ) = @_;
    my %format;
    my $num_line = 0;
    my $num_2    = 0;
    my $num_vbar = 0;
    for ( split /\n/, $file ) {
        chomp $_;

        # detect msf format
        if ( $_ =~ /MSF:\s+\d+\s+Type:\s+\S+\s+Check:\s+\d+\s+\.\./ ) {
            $format{msf} .= "$&\n";
            next;
        }
        if ( $_ =~
            /^\s+Name:\s+\S+\s+Len:\s+\d+\s+Check:\s+\d+\s+Weight:\s+\d\.\d+/ )
        {
            $format{msf} .= "$&\n";
            next;
        }

        #detect pfam format
        if ( $_ =~ /\#=GS \S+\/\d+-\d+\s+AC \S+\s*$/ ) {
            $format{pfam} .= "$&\n";
            next;
        }

        # detect fasta format
        if ( $_ =~ /^>\w+/ ) {
            $format{fasta} .= "$&\n";
            next;
        }
        if ( $_ =~ /^\#\s*STOCKHOLM/ ) {
            $format{mul} .= "$&\n";
            next;
        }

        # used to check mul format
        ++$num_line;
        my @tmp = split( /\s+/, $_ );
        ++$num_2
          if ( @tmp == 2 );
        ++$num_vbar
          if ( $_ =~ /\|/ );

        # detect other format
    }

    #determine mul format
    if (   $num_line == $num_vbar
        && $num_line == $num_2 )
    {
        $format{mul} = 'true';
    }
    elsif (
        ( $num_line == $num_vbar && $num_line != $num_2 )
        or (   $num_line != $num_vbar
            && $num_line == $num_2 )
      )
    {
        $format{mul} = 'maybe';
        print $LOG1
"the sequence might be in mul format, but not sure since the total number of sequences,\nnumber of sequences with a vertical bar '|', and number of sequences that can be splited to 2 parts by space are: $num_line, $num_vbar, $num_2, which should be same\n"
          if ($LOG1);
    }

    #get results
    if ( keys %format > 1 ) {
        print $LOG1 "***Error: more than one format is possible\n"
          if ($LOG1);
    }
    elsif ( keys %format < 1 ) {
        print $LOG1 "***Error: can not determine file format\n"
          if ($LOG1);
        return "unknown";
    }
    else {
        return ( keys %format )[0];
    }
}
#############################################################################################################################################
sub connect_db {
    use strict;
    use Carp;
    my ( $db, $db_type, $server, $user, $password ) = @_;
    die("no user name provided ")
      if ( $user eq '' );
    die("no password provided ")
      if ( $password eq '' );
    my $db_proc =
      ( DBI->connect( "dbi:$db_type:server=$server", $user, $password ) );
    carp "$DBI::errstr\n"
      if ( $DBI::errstr ne '' );
    return
      unless ( defined $db_proc );
    $db_proc->{RaiseError} = 1;    # program dies if there is an error
    $db_proc->{PrintError} = 1;
    carp "$DBI::errstr\n"
      if ( $DBI::errstr ne '' );
    $db_proc->do("use $db");
    carp "$DBI::errstr\n"
      if ( $DBI::errstr ne '' );
    $db_proc->do("set textsize 200000000");
    carp "$DBI::errstr\n"
      if ( $DBI::errstr ne '' );
    return $db_proc;
}

sub doSQL {

# new version of do_query
# $info_r: hash reference with keys: db, db_proc, db_type, server, user, password
# $sleep: -1: die,  0: return 1, > 0: sleep seconds.
    use strict;
    use Carp;
    my ( $info_r, $query, $noFetch, $sleep, $log_header, $log_r ) = @_;
    $log_header .= "--sub doSQL";

    #check input data
    if ( $$info_r{db} !~ /\S/ ) {
        $$log_r .= "***Error: $log_header--DB is not defined\n"
          if ($log_r);
        return "error";
    }
    if ( $$info_r{db_type} !~ /\S/ ) {
        $$log_r .= "***Error: $log_header--db_type is not defined\n"
          if ($log_r);
        return "error";
    }
    if ( $$info_r{server} !~ /\S/ ) {
        $$log_r .= "***Error: $log_header--server is not defined\n"
          if ($log_r);
        return "error";
    }
    if ( $$info_r{user} !~ /\S/ ) {
        $$log_r .= "***Error: $log_header--user is not defined\n"
          if ($log_r);
        return "error";
    }
    if ( $$info_r{password} !~ /\S/ ) {
        $$log_r .= "***Error: $log_header--password is not defined\n"
          if ($log_r);
        return "error";
    }

    #check connection
    while ( !$$info_r{db_proc}->{Active} ) {
        print "${\(scalar localtime)}: disconnected when processing $query\n";
        delete $$info_r{db_proc};
        $$info_r{db_proc} = &connect_db(
            $$info_r{db},   $$info_r{db_type}, $$info_r{server},
            $$info_r{user}, $$info_r{password}
        );
        unless ( $$info_r{db_proc}->{Active} ) {
            if ( $sleep < 0 ) {
                print
"***Error: $log_header--can not log in database server $$info_r{server} for user \"$$info_r{user}\".\n";
                die;
            }
            elsif ( $sleep == 0 ) {
                $$log_r .=
"***Error: $log_header--can not log in database server $$info_r{server} for user $$info_r{user}\n"
                  if ($log_r);
                return "error";
            }
            elsif ( $sleep > 0 ) {
                print
"+++Warning: $log_header--sleep for $sleep seconds due to disconnected database connection\n";
                sleep $sleep;
                my $remark =
"***Error: $log_header--can not log in database server $$info_r{server} for user $$info_r{user}, sleep for $sleep seconds now\n";
                $$log_r .= "$remark"
                  if (  $log_r
                    and $$log_r !~ /\Q$remark\E$/ );
                redo;
            }
        }
    }

    # process job
    my $statementHandle = $$info_r{db_proc}->prepare($query);
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);

        #	print "Trouble preparing this doSQL:\n$query\n";
    }

    #
    $statementHandle->execute();
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);

        #	print "Trouble executing this doSQL:\n$query\n";
    }

    #
    my $array_r = $statementHandle->fetchall_arrayref( {} )
      unless ($noFetch);
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);
    }

    #
    $statementHandle->finish;
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);
    }
    return $array_r;
}

sub do_query {
    use strict;
    use Carp;
    my ( $db_proc, $query, $noFetch, $log_header, $log_r ) = @_;
    return
      unless ( defined $db_proc );
    $log_header .= "--sub do_query";
    my $statementHandle = $db_proc->prepare($query);
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n$query\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);

        #	print "Trouble preparing this do_query:\n$query\n";
    }
    $statementHandle->execute();
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n$query\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);

        #	print "Trouble executing this do_query:\n$query\n";
    }

    #
    my $array_r = $statementHandle->fetchall_arrayref( {} )
      unless ($noFetch);
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n$query\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);
    }

    #
    $statementHandle->finish;
    if ( $DBI::errstr ne '' ) {
        carp "$DBI::errstr\n$query\n";
        $$log_r .= "***Error: $log_header--$DBI::errstr\n"
          if ($log_r);
    }
    return $array_r;
}

sub get_file_list {

# get a list of files with specified suffix
# $file_r: array reference for file names that can be obtained by <*>.
# $case: set to 'i' if suffix is case insensitive, anything else otherwise
# $suffix_r: array reference of suffixes for search
# %file: key is file name, value is hash reference, 2nd keys are suffix, value is always 1
    use strict;
    my ( $file_r, $case, $suffix_r ) = @_;
    my %file;
    for my $s (@$suffix_r) {
        for my $name (@$file_r) {
            chomp $name;
            if (
                ( $name =~ /^(.+)\.($s)$/i and $case eq 'i' )
                or (    $name =~ /^(.+)\.($s)$/
                    and $case ne 'i' )
              )
            {
                $file{$1}{$2} = 1;
            }
        }
    }
    return %file;
}

sub load_htab_line_to_evidence {

# this subroutine currently (and rather unfortunately) requires $ENV{SGC_SCRIPTS}/sgc_library.dbi
    my ( $dbproc, $db, $line, $model ) = @_;
    my ( $SEE, %SEEN, %OldHits );

    # --- Parse data for each hit
    my @newhits;
    my $curated = 0;
    my (
        $accession,   $date,         $hmm_length,   $method,
        $hmm_db,      $feat_acc,     $hmmF,         $hmmT,
        $seqF,        $seqT,         $UNASSIGNED,   $domain_score,
        $total_score, $domain_index, $domain_count, $hmm_com_name,
        $com_name,    $trusted,      $noise,        $expect_whole,
        $expect_domain
    ) = split( /\t/, $line );
    if ($model) {
        $accession = $model;
    }
    if ( !defined $HMM{$accession} ) {
        my $hmm_q =
"select trusted_cutoff, trusted_cutoff2, noise_cutoff from egad..hmm2 where hmm_acc = '$accession'";
        $HMM{$accession} = $dbproc->selectrow_hashref($hmm_q);
    }
    $feat_acc =~ s/\|.*//;
    my $feat_q =
"select feat_name from ident where (feat_name='$feat_acc' or locus='$feat_acc')";
    $feat_name = &single_sql( $dbproc, $feat_q );
    if ( !$feat_name ) {
        $feat_q =
"select feat_name from ident where (feat_name='$feat_acc' or locus='$feat_acc')";
        $feat_name = &single_sql( $dbproc, $feat_q );
    }
    if ( !$feat_name ) {
        $feat_q =
"select feat_name from nt_ident where (feat_name='$feat_acc' or locus='$feat_acc')";
        $feat_name = &single_sql( $dbproc, $feat_q );
    }
    if ( !$feat_name ) {
        warn "Couldn't find feat_name for '$feat_acc'. Skipping...";
        return;
    }

#    print LOG "\t\t$accession, $seqF, $seqT, $hmmF, $hmmT, $total_score, $trusted, $noise, $expect_whole\n";
    if ($SEE) {
        print "\n************************\n";
        print "Date: $date\n";
        print "Program: $method\n";
        print "HMM database: $hmm_db\n";
        print "Subject: $feat_name $com_name\n";
        print "Match: $accession $hmm_com_name\n";
        print "HMM length: $hmm_length\n";
        print "Subject aa coords: $seqF/$seqT\n";
        print "HMM aa coords: $hmmF/$hmmT\n";
        print "Score (domain/total): $domain_score/$total_score\n";
        print "Domain index(total): $domain_index($domain_count)\n";
        print "Cutoff (trusted/noise): $trusted/$noise\n";
        print "Expect value (domain/total): $expect_domain/$expect_whole\n\n";
    }

# If the hit is above noise or the expect is less than 1, load the hit into the db...
    if (   $total_score >= $HMM{$accession}->{'noise_cutoff'}
        || $expect_whole < 1 )
    {
        print "HMM hit is good enough to load\n"
          if ($SEE);

        #	print LOG "\t\tscore above noise or expect less than 1.\n";
        # --- Get sequence info for assembly
        my $asmbl_id = &get_asmbl_id( $dbproc, $feat_name );
        if ( !$asmbl_id ) {
            die "Can't fetch asmbl_id. Why? $feat_name, $db.";
        }

        # --- Get asmbl_ids for project
        if ( !defined $mols{$asmbl_id} ) {
            $mols{$asmbl_id} = &get_seq( $dbproc, $asmbl_id );
        }
        if ( $asmbl_id !~ /\d+/ ) {
            print "!!! ---- $feat_name may have been deleted ---- !!!\n\n";
            return;
        }
        my ($seq) = $mols{$asmbl_id};
        print "\nAssembly sequence: " . substr( $seq, 0, 15 ) . "...\n"
          if ($DEBUG);

        # Get other info about ORF
        my $feat_type;
        if ( $feat_name =~ /^NTORF/ ) {
            $feat_type = "NTORF";
        }
        else {
            $feat_type = "ORF";
        }
        print "feat_type: $feat_type\n"
          if ($DEBUG);
        my ( $end5, $end3 ) =
          &get_coords( $dbproc, $asmbl_id, $feat_name, $feat_type );
        print "\n$feat_name end5/end3: $end5/$end3\n"
          if ($SEE);

        #    print LOG "\tend5/end3: $end5/$end3\n";
        # Grab current info from db
        if ( !( $SEEN{$feat_name} ) ) {
            %OldHits = &grab_hmm_hits( $dbproc, $feat_name );
            $SEEN{$feat_name} = 1;
        }

        # --- Set curation
        my $new_state;
        if (   $total_score >= $HMM{$accession}->{'trusted_cutoff'}
            && $domain_score >= $HMM{$accession}->{'trusted_cutoff2'} )
        {
            $curated   = 1;
            $new_state = "trusted";
        }
        else {
            $curated   = 0;
            $new_state = "suspect";
        }

        # --- Find evidence end5 and end3
        my ( $EVend5, $EVend3 ) =
          &relative_to_genome( $end5, $end3, $seqF, $seqT );
        if (   $EVend5 < 0
            || $EVend3 < 0 )
        {
            print
"\n\t$feat_name, $end5, $end3, $seqF, $seqT, $EVend5, $EVend3, $asmbl_id, $feat_name, $feat_type\n\n\n";
            exit();
        }
        print "\nCalculated evidence end5/end3: $EVend5/$EVend3\n"
          if ($SEE);

        # --- See if hit is already in db
        my @matches;
        for ( my $index = 0 ; $index < $OldHits{'count'} ; $index++ ) {
            if ( $accession eq $OldHits{$index}->{'accession'} ) {
                print "$accession already in evidence\n"
                  if ($SEE);
                if (   ( $EVend5 == $OldHits{$index}->{'end5'} )
                    && ( $EVend3 == $OldHits{$index}->{'end3'} ) )
                {
                    print
"Same hit: $OldHits{$index}->{'id'}! Pushing into \@matches  \n"
                      if ($SEE);
                    push( @matches, $index );
                }
                else {
                    my ( $lo, undef, undef, $hi ) = sort { $a <=> $b } (
                        $EVend5, $EVend3,
                        $OldHits{$index}->{'end5'},
                        $OldHits{$index}->{'end3'}
                    );
                    if (
                        $hi - $lo + 1 < abs( $EVend5 - $EVend3 ) + 1 + abs(
                            $OldHits{$index}->{'end5'} -
                              $OldHits{$index}->{'end3'}
                        ) + 1
                      )
                    {
                        print
"Overlaps old hit: old $OldHits{$index}->{'end5'}/$OldHits{$index}->{'end3'} -> new $EVend5/$EVend3\nDeleting old hit."
                          if ($SEE);
                        &remove_scores( $dbproc, $OldHits{$index}->{'id'} )
                          ; # removes all feat_score rows associated with id from evidence table
                        &remove_evidence( $dbproc, $OldHits{$index}->{'id'} );
                    }
                }
            }
        }

        # If there is a single equal old hit replace the feat_score values.
        if ( @matches == 1 ) {
            my $row = $matches[0];

            # do an update and remove from hash OldHits
            print "updating evidence table\n"
              if ($SEE);
            &update_evidence_date( $db, $dbproc, $OldHits{$row}->{'id'} );

            # has there been a score state change?
            my $old_state;
            if ( $OldHits{$row}->{'score'} >=
                $HMM{$accession}->{'trusted_cutoff2'} )
            {
                $old_state = "trusted";
            }
            else {
                $old_state = "suspect";
            }
            if ( $new_state eq $old_state ) {
                $curated = $OldHits{$row}->{'curated'};
            }
            my ($curated_query) =
"update evidence set curated=$curated where id=$OldHits{$row}->{'id'}";
            &RunMod( $dbproc, $curated_query );
            print "Updating curated toggle\n"
              if ($SEE);
            print "updating feat_score table\n"
              if ($SEE);
            &remove_scores( $dbproc, $OldHits{$row}->{'id'} )
              ; # removes all feat_score rows associated with id from evidence table
            &insert_score( $dbproc, $OldHits{$row}->{'id'}, 51, $domain_score );
            &insert_score( $dbproc, $OldHits{$row}->{'id'}, 50,
                $expect_domain );

            if ( $domain_count > 1 ) {
                &insert_score( $dbproc, $OldHits{$row}->{'id'},
                    143, $total_score );
                &insert_score( $dbproc, $OldHits{$row}->{'id'},
                    144, $expect_whole );
            }

            # remove the old hit from the hash now that it's new
            delete( $OldHits{$row} );
        }
        elsif ( @matches == 0 ) {
            print "No identical old hits for $feat_name::$accession\n"
              if ($SEE);
            print "New hit.  Pushing into \@newhits  \n"
              if ($SEE);

        #	    print LOG "\t\tNew hit.  Pushing into \@newhits for insertion.\n";
            push( @newhits, $feat_name );
        }
        else {
            print
"major ERROR: $accession hits $feat_name more than once at the same coords.\n";
            print "major ERROR: Evidence 'id's are:";
            foreach my $row (@matches) {
                print " $OldHits{$row}->{'id'}";
            }
            print "\n";

#	    print LOG "\t\tERROR: $accession hits $feat_name more than once at the same coords.\n";
        }

        # For new hits, insert evidence and feat_score rows
        if ( @newhits == 1 ) {
            print "Inserting evidence for $feat_name::$accession\n"
              if ($SEE);

            # --- Insert evidence row
            my $id = &insert_evidence(
                $dbproc,    $feat_name, $hmmF,    $hmmT,
                $EVend5,    $EVend3,    $seqF,    $seqT,
                $accession, "sgc3",     $curated, "HMM2",
                "htab_loader.dbi"
            );

            # --- Generate peptide for score_text and insert
            my $fragment =
              &translate( &cut_seq2( $seq, $EVend5, $EVend3 ), $db );
            my $length = length($seq);

            #	    print LOG "\t\tInserting fragment\n";
            if ( length($fragment) < 1 ) {
                print "$feat_name, $asmbl_id:$length, $EVend5, $EVend3\n\n\n";
                exit();
            }
            &insert_alignment( $dbproc, $id, $feat_name, $fragment );

            # --- Insert feat_score rows
            #	    print LOG "\t\tInserting scores.\n";
            &insert_score( $dbproc, $id, 51, $domain_score );
            &insert_score( $dbproc, $id, 50, $expect_domain );
            if ( $domain_count > 1 ) {
                &insert_score( $dbproc, $id, 143, $total_score );
                &insert_score( $dbproc, $id, 144, $expect_whole );
            }
        }
        elsif (@newhits == 0
            && @matches < 1 )
        {
            print
"major ERROR: $accession::$feat_name No new hits!!  No old matches  FUBAR!!\n";

#	    print LOG "\t\tERROR: $accession::$feat_name No new hits!!  No old matches  FUBAR!!\n";
        }
        elsif ( @newhits > 1 ) {
            print
"major ERROR: $accession::$feat_name too many new hits  FUBAR!!\n";

#	    print LOG "\t\tERROR: $accession::$feat_name too many new hits  FUBAR!!\n";
        }
        elsif (@newhits == 0
            && @matches > 0 )
        {
            print "We be cool.\n"
              if ($SEE);
        }
        else {
            print "major ERROR: $accession::$feat_name What up? \@newhits:"
              . @newhits . "\t\n";

#	    print LOG "\t\tERROR: $accession::$feat_name What up? \@newhits:" . @newhits . "\t\n";
        }
    }
    else {

        #	print LOG "\t\tHit not above noise cutoff\n";
    }    # if ($total_score >= $noise || $expect_whole < 1)
}
####################################################################################
# read extend HMM hit with details from egad database and write data in tab-delimited format 
sub print_parse_hmm_hit {
	use strict;
    our $errorMessage;
    
	my (
        $db_proc,    $db,           $output,         $method,
        $seq_r,      $query_obj,    $prot_desc,      $query_db,
        $debug,      $noise_cutoff, $trusted_cutoff, $tb1_cutoff,
        $tb2_cutoff, $te1_cutoff,   $te2_cutoff,     $inverse,
        $expand_output
      )
      = @_;

	my $expanded_format = 0;
	if ( defined $expand_output ) {
		$expanded_format = $expand_output;
	}

    my ( $details, $com_name, $hmm_len );

    if ( $db_proc ne '' ) {
        if ( $method eq 'hmmsearch' ) {
            my $ref = &get_hmm_db_info( $db_proc, $query_obj );
			if ( !defined $ref ) {
				$errorMessage = "print_parse_hmm_hit: " . $errorMessage;
				return undef;
			}
			$hmm_len      = $$ref{hmm_len};
            $noise_cutoff = $$ref{noise_cutoff} if ( $noise_cutoff eq '' );
            $trusted_cutoff = $$ref{trusted_cutoff} if ( $trusted_cutoff eq '' );
            $com_name = $$ref{com_name};
            $com_name =~ s/\t/ /g;
            $details = $$ref{details};
            print
              "hmm_len, trusted, noise, hmm_com_name: $hmm_len, $trusted_cutoff, $noise_cutoff, $com_name\n"
              if ($debug);
        }
        elsif ( $method eq 'hmmpfam' ) {
            for my $n ( 0 .. $$seq_r{number} - 1 ) {
                my $ref = &get_hmm_db_info( $db_proc, $$seq_r{$n}{header} );
				if ( !defined $ref ) {
					$errorMessage = "print_parse_hmm_hit: " . $errorMessage;
					return undef;
				}
                $$seq_r{$n}{hmm_len} = $$ref{hmm_len};
                $$seq_r{$n}{noise_cutoff} = $$ref{noise_cutoff};
                $$seq_r{$n}{trusted_cutoff} = $$ref{trusted_cutoff};
                $$seq_r{$n}{com_name} = $$ref{com_name};
                $$seq_r{$n}{com_name} =~ s/\t/ /g;
                $$seq_r{$n}{details} = $$ref{details};
                print
                  "seq \#, hmm_len, trusted, noise, hmm_com_name: $n, $$seq_r{$n}{trusted_cutoff}, $$seq_r{$n}{noise_cutoff}, $$seq_r{$n}{comment}\n"
                  if ($debug);
            }
        }
    }

    if ( $$seq_r{number} > 0 ) {
        my ( $day, $month, $year ) =
          ( (localtime)[3], (localtime)[4] + 1, (localtime)[5] + 1900 );
        print $output
          "1 HMM\t2 htab_date\t3 hmm_length\t4 method\t5 query_DB\t6 protein\t7 hmm_f\t8 hmm_t\t9 protein_f\t10 protein_t\t11 NULL\t12 domain_score\t13 total_score\t14 domain \#\t15 \# of domains\t16 HMM description\t17 protein description\t18 trusted_cutoff\t19 noise_cutoff\t20 total_E_value\t21 domain_E_value\n"
          if ($debug);
        for my $m ( 0 .. $$seq_r{number} - 1 ) {
            if ($inverse) {
                next
                  if ( $tb1_cutoff < $$seq_r{$m}{score}
                    or $te1_cutoff > $$seq_r{$m}{e_value} );
            }
            else {
                next
                  if ( $tb1_cutoff > $$seq_r{$m}{score}
                    or $te1_cutoff < $$seq_r{$m}{e_value} );
            }
            for my $n ( 1 .. $$seq_r{$m}{number} ) {
                if ( !$$seq_r{$m}{$n}{seq_f} ) {
                    next;
                }
                if ($inverse) {
                    next
                      if ( $tb2_cutoff < $$seq_r{$m}{$n}{score}
                        or $te2_cutoff > $$seq_r{$m}{$n}{e_value} );
                }
                else {
                    next
                      if ( $tb2_cutoff > $$seq_r{$m}{$n}{score}
                        or $te2_cutoff < $$seq_r{$m}{$n}{e_value} );
                }
                if ( $method eq 'hmmsearch' ) {
                	if ( $query_obj =~ /_rev$/ ) {			# check for strand specific HMM (_fwd or _rev)
                		$query_obj =~ s/_rev$//;
                		my $seq_t = $$seq_r{$m}{$n}{seq_t};
                		$$seq_r{$m}{$n}{seq_t} = $$seq_r{$m}{$n}{seq_f};
                		$$seq_r{$m}{$n}{seq_f} = $seq_t;
                	} elsif ( $query_obj =~ /_fwd$/ ) {
                		$query_obj =~ s/_fwd$//;
                	}
					my @tab = ($query_obj, "$month-$day-$year", $hmm_len, $method, $query_db, $$seq_r{$m}{header}, $$seq_r{$m}{$n}{hmm_f}, $$seq_r{$m}{$n}{hmm_t}, $$seq_r{$m}{$n}{seq_f}, $$seq_r{$m}{$n}{seq_t}, undef, $$seq_r{$m}{$n}{score}, $$seq_r{$m}{score}, $n, $$seq_r{$m}{number}, $com_name, $$seq_r{$m}{comment}, $trusted_cutoff, $noise_cutoff, $$seq_r{$m}{e_value}, $$seq_r{$m}{$n}{e_value} );
					if ( $expanded_format ) {
						push ( @tab, $details );
					}
					print $output join("\t", @tab ) . "\n";
                }
                elsif ( $method eq 'hmmpfam' ) {
                	if ( $$seq_r{$m}{header} =~ /_rev$/ ) {	# check for strand specific HMM (_fwd or _rev)
                		$$seq_r{$m}{header} =~ s/_rev$//;
                		my $seq_t = $$seq_r{$m}{$n}{seq_t};
                		$$seq_r{$m}{$n}{seq_t} = $$seq_r{$m}{$n}{seq_f};
                		$$seq_r{$m}{$n}{seq_f} = $seq_t;
                	} elsif ( $$seq_r{$m}{header} =~ /_fwd$/ ) {
                		$$seq_r{$m}{header} =~ s/_fwd$//;
                	}
                	my @tab = ( $$seq_r{$m}{header}, "$month-$day-$year", $$seq_r{$m}{hmm_len}, $method, $query_db, $query_obj, $$seq_r{$m}{$n}{hmm_f}, $$seq_r{$m}{$n}{hmm_t}, $$seq_r{$m}{$n}{seq_f}, $$seq_r{$m}{$n}{seq_t}, undef, $$seq_r{$m}{$n}{score}, $$seq_r{$m}{score}, $n, $$seq_r{$m}{number}, $$seq_r{$m}{com_name}, $prot_desc, $$seq_r{$m}{trusted_cutoff}, $$seq_r{$m}{noise_cutoff}, $$seq_r{$m}{e_value}, $$seq_r{$m}{$n}{e_value} );
                	if ( $expanded_format ) {
                		push( @tab, $$seq_r{$m}{details} );
                	}
                    print $output join("\t", @tab) . "\n";
                }
            }
        }
    }
}

####################################################################################
# read HMM's details from egad database
sub get_hmm_db_info {
	use strict;
	my ( $db_proc, $acc ) = @_;
	our $errorMessage;

# get HMM details from egad database
# if database not available, return empty details
	if ( ! defined $db_proc ) {
		return &emptyHMMDetail;
	}

# first try hmm2 table
	my $info =
		$db_proc->selectrow_hashref(
			"select hmm_acc, hmm_len, trusted_cutoff, noise_cutoff, gathering_cutoff, hmm_com_name as com_name, "
				. "gene_sym, ec_num, iso_type, trusted_cutoff2 "
				. "from hmm2 "
				. "where hmm_acc = ? and is_current = 1",
			undef, $acc );
	if ( defined $$info{hmm_acc} && $$info{hmm_acc} eq $acc ) {
		$$info{details} = "gene_sym=" . $$info{gene_sym}
			. "~~ec_num=" . $$info{ec_num}
			. "~~gathering_cutoff=" . $$info{gathering_cutoff}
			. "~~trusted_cutoff2=" . $$info{trusted_cutoff2}
			. "~~iso_type=" . $$info{iso_type};
		delete ( $$info{gene_sym} );
		delete ( $$info{ec_num} );
		delete ( $$info{gathering_cutoff} );
		delete ( $$info{trusted_cutoff2} );
		delete ( $$info{iso_type} );
		return $info;
	}

# if no results, try rfam table
	my $acc2 = $acc;
	$acc2 =~ s/_rev$//;
	$acc2 =~ s/_fwd$//;
	$info =
		$db_proc->selectrow_hashref(
			"select accession as hmm_acc, null as hmm_len, trusted_cutoff, noise_cutoff, gathering_thresh, com_name, "
				. "window_size, feat_type, feat_class, gene_sym "
				. "from rfam "
				. "where accession = ? and iscurrent = 1",
			undef, $acc2 );
	if ( defined $$info{hmm_acc} && $$info{hmm_acc} eq $acc2 ) {
		$$info{details} = "gene_sym=" . $$info{gene_sym}
			. "~~feat_type=" . $$info{feat_type}
			. "~~feat_class=" . $$info{feat_class}
			. "~~gathering_thresh=" . $$info{gathering_thresh}
			. "~~window_size=" . $$info{window_size};
		delete ( $$info{gene_sym} );
		delete ( $$info{feat_type} );
		delete ( $$info{feat_class} );
		delete ( $$info{gathering_thresh} );
		delete ( $$info{window_size} );
		return $info;
	}

# not in db, return empty details
	return &emptyHMMDetail;
}

sub emptyHMMDetail {
	use strict;
# return empty HMM details, used when db query failed
	my $empty;
	$$empty{hmm_acc} = undef;
	$$empty{hmm_len} = undef;
	$$empty{trusted_cutoff} = undef;
	$$empty{noise_cutoff} = undef;
	$$empty{com_name} = undef;
	$$empty{details} = undef;
	return $empty;
}
	
####################################################################################
# read HMM results file,
# parse and extend details from egad database,
# and write tab-delimited output
sub parse_hmm_hits {
    use strict;
    my (
        $file,       $db,           $db_proc,        $usage,
        $debug,      $align_output, $align_file,     $b1_cutoff,
        $b2_cutoff,  $e1_cutoff,    $e2_cutoff,      $format,
        $quiet,      $noise_cutoff, $trusted_cutoff, $output_dir,
        $tb1_cutoff, $tb2_cutoff,   $te1_cutoff,     $te2_cutoff,
        $inverse,    $multi, $expand_output
      )
      = @_;
	

	my $expanded_format = 0;
	if ( defined $expand_output ) {
		$expanded_format = $expand_output;
	}

    #   following are entry-wide variables
    if ( $inverse ne '' ) {
        $tb1_cutoff = (
              $tb1_cutoff ne ''
            ? $tb1_cutoff
            : 2000
        );
        $tb2_cutoff = (
              $tb2_cutoff ne ''
            ? $tb2_cutoff
            : 2000
        );
        $te1_cutoff = (
              $te1_cutoff ne ''
            ? $te1_cutoff
            : -100
        );
        $te2_cutoff = (
              $te2_cutoff ne ''
            ? $te2_cutoff
            : -100
        );
    }
    else {
        $inverse    = 0;
        $tb1_cutoff = (
              $tb1_cutoff ne ''
            ? $tb1_cutoff
            : -2000
        );
        $tb2_cutoff = (
              $tb2_cutoff ne ''
            ? $tb2_cutoff
            : -2000
        );
        $te1_cutoff = (
              $te1_cutoff ne ''
            ? $te1_cutoff
            : 100
        );
        $te2_cutoff = (
              $te2_cutoff ne ''
            ? $te2_cutoff
            : 100
        );
    }
    my ( @line, $method, $query_obj, $query_db, %seq, $output, $prot_desc,
        $fhtab );
    my $version;

    # drink in the output from file or stdin
    if ( $file ne '' ) {
        chomp $file;
        open( FH, "$file" )
          || die "Can't open $file for reading: $!\n";
        @line = <FH>;
        close FH;
        $fhtab = "$file.htab";
        $fhtab =~ s:^(.*/)?:$output_dir/:
          if ($output_dir);
        open( FHTAB, ">$fhtab" )
          || die "Can't open $fhtab for writing: $!\n";
        $output = \*FHTAB;
    }
    else {
        @line   = <STDIN>;
        $output = \*STDOUT;
    }

    #parse file content
    my $n = 0;
    while ( $n < @line ) {
		if ( $line[$n] =~ /^ *RF [ Xx]+$/ ) {
        #get program used
		} elsif ( $line[$n] =~ /hmmsearch/ ) {
            $method = 'hmmsearch';
            print "Method: $method\n"
              if ($debug);
        }
        elsif ( $line[$n] =~ /hmmpfam/ ) {
            $method = 'hmmpfam';
            print "Method: $method\n"
              if ($debug);
        }
        else {
            if ( $method eq '' ) {
                if ( \*FHTAB ) {
                    close FHTAB;
                    unlink "$file.htab"
                      if ( !-s "$file.htab" );
                }
                die "***Error: unknown search program used: $line[$n]\n";
            }
        }
        if ( $line[$n] =~ /^HMMER (\d)\.(\S+) \(/ ) {
            my ( $maj_ver, $min_ver ) = ( $1, $2 );
            print "Version: $maj_ver.$min_ver\n"
              if ($debug);
            if ( $maj_ver == 2 ) {
                if ( $min_ver eq '1.1' ) {
                    $version = 1;
                }
                else {
                    $version = 2;
                }

#                } elsif ($min_ver eq '2g') {
#                    $version = 2;
#                } else {
#                    if (\*FHTAB) {
#                        close FHTAB;
#                        unlink "$file.htab" if(!-s "$file.htab");
#                    }
#                    die "***Error: invalid HMMer version '$maj_ver.$min_ver' used.\n";
#                }
            }
            else {
                if ( \*FHTAB ) {
                    close FHTAB;
                    unlink "$file.htab"
                      if ( !-s "$file.htab" );
                }
                die
                  "***Error: invalid HMMer version '$maj_ver.$min_ver' used.\n";
            }
        }
        elsif ( $line[$n] =~ /^Logical Depth HMMER (\d+)\.(\S+)/ ) {
            my ( $maj_ver, $min_ver ) = ( $1, $2 );
            print "Version: $maj_ver.$min_ver\n"
              if ($debug);
            $version =
              2;  # Logical Depth HMMer currently uses 2.2g output format
        }
        # New code for CLC hmmer handling 2009-04-08 by aklump
		  elsif ( $line[$n] =~ /^CLC\w*\s+(\d+)\.(\S+)/ ){
            my ( $maj_ver, $min_ver ) = ( $1, $2 );
            print "Version: $maj_ver.$min_ver\t"
              if ($debug);
            print "Assuming CLC is still HMMer2 compliant..if things break, verify this first.\n"
              if ($debug);
				$version = 2;
		  }

        #get sequence db name for hmmsearch
        if ( $line[$n] =~ /^Sequence database:\s+([\w\-.]+)/ ) {
            $query_db = $1;
            print "Sequence database: $query_db\n"
              if ($debug);
        }

        #get HMM db name for hmmpfam
        if ( $line[$n] =~ /^HMM file:\s+/ ) {
            if ( $' ne '' ) {
                $query_db = $';
                chomp $query_db;
                $query_db =~ s/\t//g;
                $query_db =~ s/^\s+//;
                $query_db =~ s/\s+$//;
                print "HMM database: $query_db\n"
                  if ($debug);
            }
        }
        if ( $version == 1 ) {
            if (
                (
                       $line[$n] =~ /^Query HMM: ([\w.\-]+)\|\|/
                    or $line[$n] =~ /^Query:\s+(\S+)/
                )
                && $seq{number} > 0
              )
            {

                #get hmm length from either DB or HMM file
                &print_parse_hmm_hit(
                    $db_proc,      $db,             $output,
                    $method,       \%seq,           $query_obj,
                    $prot_desc,    $query_db,       $debug,
                    $noise_cutoff, $trusted_cutoff, $tb1_cutoff,
                    $tb2_cutoff,   $te1_cutoff,     $te2_cutoff,
                    $inverse, $expanded_format
                  )
                  unless ($quiet);
                $query_obj = $prot_desc = "";
                %seq = ();
            }

            #get query hmm name for hmmsearch
            if ( $line[$n] =~ /^Query HMM: ([\w.\-]+)\|/ ) {
                $query_obj = $1;
                print "query HMM: $query_obj\n"
                  if ($debug);
            }

            #get query sequence name for hmmpfam
            if ( $line[$n] =~ /^Query:\s+(\S+)/ ) {
                $query_obj = $1;
                $prot_desc = '';
                if ( $' ne '' ) {
                    $prot_desc = $';
                    chomp $prot_desc;
                    $prot_desc =~ s/\t/ /g;
                    $prot_desc =~ s/^\s+//;
                    $prot_desc =~ s/\s+$//;
                }
                print "query seq: $query_obj\n"
                  if ($debug);
            }
        }
        elsif ( $version == 2 ) {

            #get query sequence name for hmmpfam
            if ( $line[$n] =~ /^Query sequence:\s+(\S+)/ ) {
                my $tmp_obj = $1;
                &print_parse_hmm_hit(
                    $db_proc,      $db,             $output,
                    $method,       \%seq,           $query_obj,
                    $prot_desc,    $query_db,       $debug,
                    $noise_cutoff, $trusted_cutoff, $tb1_cutoff,
                    $tb2_cutoff,   $te1_cutoff,     $te2_cutoff,
                    $inverse, $expanded_format
                  )
                  unless ($quiet);
                $query_obj = $tmp_obj;
                $prot_desc = "";
                %seq       = ();
                print "query seq: $query_obj\n"
                  if ($debug);
            }

            #get query hmm name for hmmsearch
            if ( $line[$n] =~ /^Query HMM:\s+(\S+)/ ) {
                $query_obj = $1;
                print "query HMM: $query_obj\n"
                  if ($debug);
            }
            if ( $line[$n] =~ /^Description:\s+(.+)/ ) {
                $prot_desc = $1;
                chomp $prot_desc;
                $prot_desc =~ s/\t/ /g;
                $prot_desc =~ s/^\s+//;
                $prot_desc =~ s/\s+$//;
            }
        }

#get header, total score, e_value, number of domains of sequences for hmmsearch or hmmpfam
# 2007-05-08 RAR altered regex to accept either N or #D as number of domains column header
        if (   $line[$n] =~ /^(Sequence\s+)Description\s+Score\s+E-value\s+(?:N|\#D)/
            || $line[$n] =~ /^(Model\s+)Description\s+Score\s+E-value\s+(?:N|\#D)/ )
        {
            my $m            = 0;
            my $seqid_length = length($1);

            # if these are results from searches on a multi-sequence input file
            # this return should not happen
            unless ($multi) {

                # create one line file if there are no hits.
                if ( $line[ $m + $n + 2 ] =~ /\[no hits above thresholds\]/ )
                {
                    print $output "No hits above thresholds.\n" ;
                    ###return ();   #### <---- this line was causing the program to exit even if other files followed the one with no hits!!!
		                    #### Removed by selengut 6/11/07
                }
            }
            while ( $line[ $m + $n + 2 ] =~ /^\S+/ ) {
                my @data = split /\s+/, $line[ $m + $n + 2 ];
                $seq{$m}{header}  = shift @data;
                $seq{$m}{number}  = pop @data;
                $seq{$m}{e_value} = pop @data;
                $seq{$m}{score}   = pop @data;
                $seq{$m}{comment} = join " ", @data;
                print
                  "$m --- $seq{$m}{header}\t$seq{$m}{comment}\t$seq{$m}{score}\t$seq{$m}{e_value}\t$seq{$m}{number}\n"
                  if ($debug);
                ++$m;
            }
            $seq{number} = $m;
        }

#get seq-f, seq-t, hmm-f, hmm-t, score, e_value for each domain for hmmsearch or hmmpfam
        if ( $line[$n] =~
            /^Sequence\s+Domain\s+seq-f\s+seq-t\s+hmm-f\s+hmm-t\s+score\s+E-value/i
            or $line[$n] =~
            /^Model\s+Domain\s+seq-f\s+seq-t\s+hmm-f\s+hmm-t\s+score\s+E-value/i
          )
        {
            my $m = 0;
            while ( $line[ $m + $n + 2 ] =~
                /^(\S+)\s+(\d+)\/(\d+)\s+(\d+)\s+(\d+).+?(\d+)\s+(\d+).+?(-?\d\S*)\s+(\d\S*)/
              )
            {
                for my $k ( 0 .. $seq{number} - 1 ) {
                    if ( $seq{$k}{header} eq $1 ) {
                        $seq{$k}{$2}{seq_f}   = $4;
                        $seq{$k}{$2}{seq_t}   = $5;
                        $seq{$k}{$2}{hmm_f}   = $6;
                        $seq{$k}{$2}{hmm_t}   = $7;
                        $seq{$k}{$2}{score}   = $8;
                        $seq{$k}{$2}{e_value} = $9;
                        print "***Error: number of domain do not match\n"
                          if (  $3 != $seq{$k}{number}
                            and $debug );
                        print "$m --- $1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\n"
                          if ($debug);
                    }
                }
                ++$m;
            }
        }

        # get aligned sequences ---------------------------------------------------
        if (   $align_output
            && $line[$n] =~ /^Alignments of top-scoring domains:/
            && $method eq 'hmmsearch' )
        {
            my ( %align, %posi,  %align_seq );
            my ( $s_hmm, $e_hmm, $hmm_max_len ) = ( 10000, 0, 0 );

            #get the maximum range of hmm sequnces, required for fragment model
            for my $k1 ( 0 .. $seq{number} - 1 ) {
                for my $k2 ( 1 .. $seq{$k1}{number} ) {
                    $s_hmm = $seq{$k1}{$k2}{hmm_f}
                      if ( $seq{$k1}{$k2}{hmm_f} < $s_hmm );
                    $e_hmm = $seq{$k1}{$k2}{hmm_t}
                      if ( $seq{$k1}{$k2}{hmm_t} > $e_hmm );
                }
            }

            #get hmm and prot sequences
            my $m = 1;
            while ( $line[ $n + $m ] !~ /^Histogram of all scores:/ ) {
                if ( $line[ $n + $m ] =~ /^\S+/ ) {
                    for my $k ( 0 .. $seq{number} - 1 ) {
                        next
                          if ( $b1_cutoff > $seq{$k}{score}
                            or $e1_cutoff < $seq{$k}{e_value} );
                        if ( $line[ $n + $m ] =~
                            /^\Q$seq{$k}{header}\E:\s+domain\s+(\d+)/ )
                        {
                            my $domain = $1;
 
  							# ignore extra line generated when HMM's RF flag is on
							if ( $line[ $n + $m + 1 ] =~ /^ *RF [ Xx]+$/ ) { $m++ }

                            next
                              if ( $b2_cutoff > $seq{$k}{$domain}{score}
                                or $e2_cutoff < $seq{$k}{$domain}{e_value} );
                            while ( $line[ $m + $n + 1 ] =~ /^\s+/
                                and $line[ $m + $n + 3 ] =~
                                /^\s+([^\|\.\:\/\s]+((\||\.|\:|\/)[^\|\.\:\/\s]+)?)/
                                && index( $seq{$k}{header}, $1 ) == 0 )
                            {
                                $align{$k}{$domain}{prot} .= $3
                                  if ( $line[ $m + $n + 3 ] =~
                                    /^\s+(\S+)\s+(\d+|-)\s+(\S+)\s+(\d+|-)/
                                  );
                                $align{$k}{$domain}{hmm} .=
                                  $line[ $m + $n + 1 ];
                                $m += 4;
								if ( $line[ $n + $m + 1 ] =~ /^ *RF [ Xx]+$/ ) { $m++ }
                            }
                            $align{$k}{$domain}{hmm} =~ s/[\s\n]//g;
                            $align{$k}{$domain}{hmm} =~ s/\*->|<-\*//g
                              ; # can not merge since <-* might be broken to two sections
                            $align{$k}{$domain}{show} = 1
                              ; # means this line is above all cutoff and should be printed out later
                            die
                              "***Error: wrong character(s) \"$&\" in HMM sequence $align{$k}{$domain}{hmm} for file $file\n"
                              if (
                                $align{$k}{$domain}{hmm} =~ /[^a-zA-Z.]/ );
                            last;
                        }
                    }
                }
                ++$m;
            }

#align hmm and prot sequences to the maximum range of hmm sequence as defined by $s_hmm, $e_hmm, gap in prot
# sequences is dot, gap in hmm sequence is dash since dot means insertion which needs to be processed later.
            for my $k1 ( 0 .. $seq{number} - 1 ) {
                for my $k2 ( 1 .. $seq{$k1}{number} ) {

                    # what if header has coords in it???????????
                    if (   $seq{$k1}{$k2}{hmm_f} > 0
                        && $seq{$k1}{$k2}{hmm_t} > $seq{$k1}{$k2}{hmm_f} )
                    {
                        $align{$k1}{$k2}{new_prot} =
                            '.' x ( $seq{$k1}{$k2}{hmm_f} - $s_hmm )
                          . $align{$k1}{$k2}{prot}
                          . '.' x ( $e_hmm - $seq{$k1}{$k2}{hmm_t} );
                        $align{$k1}{$k2}{new_hmm} =
                            '-' x ( $seq{$k1}{$k2}{hmm_f} - $s_hmm )
                          . $align{$k1}{$k2}{hmm}
                          . '-' x ( $e_hmm - $seq{$k1}{$k2}{hmm_t} );
                        $align{$k1}{$k2}{new_prot} =~ s/-/\./g;

# following can not be replaced by ($e_hmm - $s_hmm +1) since hmm sequence contains dor as gap
                        $hmm_max_len = length $align{$k1}{$k2}{new_hmm}
                          if (
                            length $align{$k1}{$k2}{new_hmm} >
                            $hmm_max_len );

                        #			print "A $align{$k1}{$k2}{new_prot}\n";
                        #			print "B $align{$k1}{$k2}{new_hmm}\n";
                    }
                }
            }

            #get position of gaps in hmm sequences and the indices of the hmm sequences
            for my $n ( 1 .. $hmm_max_len - 1 ) {  # $n start from 1 here
                for my $k1 ( 0 .. $seq{number} - 1 ) {
                    for my $k2 ( 1 .. $seq{$k1}{number} ) {
                        if (
                            substr( $align{$k1}{$k2}{new_hmm}, $n, 1 ) eq
                            '.' )
                        {
                            ++$posi{ $n - 1 }{$k1}{$k2};
                            $posi{ $n - 1 }{max} = $posi{ $n - 1 }{$k1}{$k2}
                              if ( $posi{ $n - 1 }{max} <
                                $posi{ $n - 1 }{$k1}{$k2} );
                            $align{$k1}{$k2}{new_hmm} =
                                substr( $align{$k1}{$k2}{new_hmm}, 0, $n )
                              . substr( $align{$k1}{$k2}{new_hmm}, $n + 1 );
                            redo;
                        }
                    }
                }
            }

            #re-align protein sequences based on gaps in hmm sequences
            my $count = 0;
            for my $n ( sort { $a <=> $b } keys %posi ) {
                for my $k1 ( 0 .. $seq{number} - 1 ) {
                    for my $k2 ( 1 .. $seq{$k1}{number} ) {
                        if ( $posi{$n}{max} > 0 ) {
                            $align{$k1}{$k2}{new_prot} =
                              substr( $align{$k1}{$k2}{new_prot},
                                0, $n + $count + 1 )
                              . '.' x
                              ( $posi{$n}{max} - $posi{$n}{$k1}{$k2} )
                              . substr( $align{$k1}{$k2}{new_prot},
                                $n + $count + 1 );
                        }
                    }
                }
                $count += $posi{$n}{max};
            }

            #remove a gap if all sequences have it
            my $n = 0;
            while ( $n < length $align{0}{1}{new_prot} ) {
                my $dele = 1;
                for my $k1 ( 0 .. $seq{number} - 1 ) {
                    for my $k2 ( 1 .. $seq{$k1}{number} ) {
                        $dele = 0
                          if (
                            substr( $align{$k1}{$k2}{new_prot}, $n, 1 ) ne
                            '.' );
                    }
                }
                if ($dele) {
                    for my $k1 ( 0 .. $seq{number} - 1 ) {
                        for my $k2 ( 1 .. $seq{$k1}{number} ) {
                            $align{$k1}{$k2}{new_prot} =
                                substr( $align{$k1}{$k2}{new_prot}, 0, $n )
                              . substr( $align{$k1}{$k2}{new_prot}, $n + 1 );
                        }
                    }
                }
                else {
                    ++$n;
                }
            }

            #output
            $m = -1;
            for my $k1 ( 0 .. $seq{number} - 1 ) {
                for my $k2 ( 1 .. $seq{$k1}{number} ) {
                    if ( $align{$k1}{$k2}{show} ) {
                        ++$m;
                        $align_seq{$m}{ori} =
                          "$seq{$k1}{header}\t$align{$k1}{$k2}{new_prot}";
                        &parse_single_seq( $align_seq{$m}, 'mul', '', '', '',
                            '', '' );
                        $align_seq{$m}{lend} = $seq{$k1}{$k2}{seq_f};
                        $align_seq{$m}{rend} = $seq{$k1}{$k2}{seq_t};
                    }
                }
            }
            $align_seq{number} = $m + 1;

            #get more header info for each protein
            my %prot;
            for my $n ( 0 .. $align_seq{number} - 1 ) {
                (
                      $align_seq{$n}{db} eq "EGAD"
                    ? $prot{"$align_seq{$n}{db}|$align_seq{$n}{prot_id}"} =
                      {}
                    : $prot{"$align_seq{$n}{db}|$align_seq{$n}{0}"} = {}
                );
            }

            #	    &download_nraa(\%prot, '', '', 1);
            #	    &yank_prot(\%prot, '', '');
            for my $n ( 0 .. $align_seq{number} - 1 ) {
                if ( $align_seq{$n}{db} eq "EGAD" ) {
                    $align_seq{$n}{locus} = $1
                      if (
                        $prot{"$align_seq{$n}{db}|$align_seq{$n}{prot_id}"}
                        {first_acc} =~ /^EGAD\|\d+\|([^\|\/\s]+)/ );
                }
                else {
                    $align_seq{$n}{1} = $1
                      if ( $prot{"$align_seq{$n}{db}|$align_seq{$n}{0}"}
                        {first_acc} =~ /^[A-Za-z]+\|[^\|]+\|([^\|\/\s]+)/ );
                    $align_seq{$n}{number} = 2;
                }
            }

            #output results
            my $out =
              &format_sequence( \%align_seq, $format, '', 60, 1, '', '',
                '' );
            if ( $align_file ne '' ) {
                $align_file .= ".msf"
                  if ( $format eq 'msf' );
                $align_file .= ".fa"
                  if ( $format eq 'fasta' );
                $align_file .= ".mul"
                  if ( $format eq 'mul' );
                open( FH, ">$align_file" );
                print FH $out;
                close FH;
                unlink "$align_file"
                  if ( -s "$align_file" == 0 );
            }
            else {
                open( FH, ">$file.aligned" );
                print FH $out;
                close FH;
                unlink "$file.aligned"
                  if ( -s "$file.aligned" == 0 );
            }

            #end of job
            $n += $m;
        }
        ++$n;
    }
    &print_parse_hmm_hit(
        $db_proc,    $db,           $output,         $method,
        \%seq,       $query_obj,    $prot_desc,      $query_db,
        $debug,      $noise_cutoff, $trusted_cutoff, $tb1_cutoff,
        $tb2_cutoff, $te1_cutoff,   $te2_cutoff,     $inverse,
        $expanded_format
      )
      unless ($quiet);
    close FHTAB;
    unlink "$file.htab"
      if ( !-s "$file.htab" );
    return "success";
}
1;
