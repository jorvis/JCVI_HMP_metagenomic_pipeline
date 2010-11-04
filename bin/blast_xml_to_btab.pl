#!/usr/local/bin/perl
use strict;

# constant used to cnvert evalue to pvalue
my $euler = 2.718281828;

# define output
my @hitcols = split( /,/,
	    "query_id,run_date,query_length,blastpgm,blastdb,subject_id,"
	  . "query_end5,query_end3,subject_left,subject_right,"
	  . "pct_identity,pct_similarity,bit_score,"
	  . "skip_offset_begin,skip_offset_end,subject_definition,"
	  . "skip_frame,query_strand,subject_length,evalue,pvalue" );

# define XML parsing
my %xmlattr;
$xmlattr{"BlastOutput_program"} = "blastpgm";
$xmlattr{"BlastOutput_db"}      = "blastdb";
$xmlattr{"Iteration_query-def"} = "query_definition";
$xmlattr{"Iteration_query-len"} = "query_length";
$xmlattr{"Hit_def"}             = "subject_definition";
$xmlattr{"Hit_len"}             = "subject_length";
$xmlattr{"Hsp_bit-score"}       = "bit_score";
$xmlattr{"Hsp_evalue"}          = "evalue";
$xmlattr{"Hsp_query-from"}      = "query_end5";
$xmlattr{"Hsp_query-to"}        = "query_end3";
$xmlattr{"Hsp_hit-from"}        = "subject_left";
$xmlattr{"Hsp_hit-to"}          = "subject_right";
#$xmlattr{"Hsp_query-frame"}     = "query_frame";
#$xmlattr{"Hsp_hit-frame"}       = "subject_frame";
$xmlattr{"Hsp_identity"}        = "num_identical";
$xmlattr{"Hsp_positive"}        = "num_similar";
#$xmlattr{"Hsp_gaps"}            = "num_gaps";
$xmlattr{"Hsp_align-len"}       = "alignment_length";

# parse XML
my $blastpgm;
my $blastdb;
my $hsp = newHsp();
while ( my $line = <STDIN> ) {
	$line =~ s/[\n\r]/ /g;
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	$line =~ s/&lt;/</g;
	$line =~ s/&gt;/>/g;
	$line =~ s/&amp;/&/g;
	$line =~ s/&quot;/"/g;
	$line =~ s/&apos;/'/g;

	# check XML tag
	if ( substr( $line, 0, 1 ) eq "<" ) {
		$line = substr( $line, 1 );
		my ($xml_tag) = split( ">", $line );
		my $attrname = $xmlattr{"$xml_tag"};
		my $attrval;

		# get tag value
		if ( defined $attrname ) {
			$line = substr( $line, length($xml_tag) + 1 );
			my $close_tag = "</" . $xml_tag . ">";
			my $eod = index( $line, $close_tag );
			if ( $close_tag == -1 ) {
				$attrval = $line;
			}
			else {
				$attrval = substr( $line, 0, $eod );
			}

			# save tag value
			if ( $attrname eq "blastpgm" ) {
				$blastpgm = $attrval;
			}
			elsif ( $attrname eq "blastdb" ) {
				$blastdb = $attrval;
			}
			elsif ( $attrname eq "subject_definition" ) {
				my $subject_definition = $attrval;
				$$hsp{$attrname} = $attrval;
				my ($subject_id) = split( /\s/, $subject_definition );
#				$subject_id = parseBlastSeqId($subject_id);
				$$hsp{"subject_id"} = $subject_id;
			}
			elsif ( $attrname eq "query_definition" ) {
				my $subject_definition = $attrval;
				my $query_definition = $attrval;
				my ($query_id) = split( /\s/, $query_definition );
#				$query_id = parseBlastSeqId($query_id);
				$$hsp{"query_id"} = $query_id;
			}
			else {
				$$hsp{$attrname} = $attrval;
			}
		}

		# save hit
		if ( $attrname eq "alignment_length" ) {
			
			# save blast program and database
			$$hsp{blastpgm} = $blastpgm;
			$$hsp{blastdb} = $blastdb;

			# calculate orientation
			my $orientation = 1;
			if ( $$hsp{query_frame} < 0 ) {
				$orientation = -$orientation;
			}
			if ( $$hsp{subject_frame} < 0 ) {
				$orientation = -$orientation;
			}
			if ( $orientation < 0 ) {
				$$hsp{query_strand} = "Minus";
			} else {
				$$hsp{query_strand} = "Plus";
			}

			# calculate percentages
			$$hsp{pvalue} = 1. - $euler**( -$$hsp{evalue} );
			$$hsp{pct_identity} = int(
				1000. * $$hsp{num_identical} / ( $$hsp{alignment_length} ) ) / 10.;
			$$hsp{pct_similarity} = int(
				1000. * $$hsp{num_similar} / ( $$hsp{alignment_length} ) ) / 10.;

			# write hsp to STDOUT
			my @row = extractArrayFromHash( $hsp, @hitcols );
			print STDOUT join( "\t", @row ) . "\n";

			# start a new hsp
			$hsp = newHsp($hsp);
		}
	}
}
exit(0);

# attempt to find best identifier in sequence identifier string
sub parseBlastSeqId {
	my ($rawid) = @_;

	$rawid =~ s/,\s*$//;
	$rawid =~ s/\|\s*$//;
	my @id = split( /\|/, $rawid );
	if ( scalar @id < 2 ) { return $rawid }
	my $i = 0;
	while ( $i < scalar @id - 1 ) {
		if (
			index( ".gi.gb.rf.emb.dbj.pir.prf.sp.ref.",
				"." . lc( $id[$i] ) . "." ) >= 0
		  )
		{
			return lc( $id[$i] ) . "|" . $id[ $i + 1 ];
		}
		$i++;
	}
	return $rawid;
}

sub newHsp {
	my ($oldhsp) = @_;
	my %hsp;
	if ( defined $oldhsp ) {
		%hsp = %$oldhsp;
	}
	else {
		$hsp{query_id} = "";
	}
	$hsp{num_gaps} = 0;
	$hsp{run_date} = today();

	return \%hsp;
}

sub extractArrayFromHash{
	my ( $inhash, @keys ) = @_;
	
	my @outarray = ();
	for my $key ( @keys ) {
		if ( defined $inhash->{$key} ) {
			push(@outarray,$inhash->{$key});
		} else {
			push(@outarray,undef);
		}
	}
	return @outarray;
}

sub today {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	
	$mday = &lpad($mday,2,"0");
	$mon = &lpad($mon+1,2,"0");
	$year = &lpad($year+1900,4,"0");
	$year+=100;
	my $today = "$year-$mon-$mday";
	return $today;
}

sub lpad {
	my ( $text, $pad_len, $pad_char) = @_;

	if ( $pad_len<=0 ) {
		return "";
	} elsif ( $pad_len<length($text) ) {
		return substr($text,0,$pad_len);
	}

	if ( !defined $pad_char ) {
		$pad_char = " ";
	} elsif ( length($pad_char)>1 ) {
		$pad_char = substr($pad_char,0,1);
	}

	if ( $pad_len>length($text) ) {
		$text = $pad_char x ( $pad_len - length( $text ) ). $text;
	}
	
	return "$text"; 
}
