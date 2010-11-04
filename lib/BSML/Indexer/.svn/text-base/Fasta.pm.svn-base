package BSML::Indexer::Fasta;

=head1 NAME

BSML::Indexer::Fasta;

=head1 VERSION

This document refers BSML::Indexer::Fasta.pm.

=head1 SYNOPSIS

  use BSML::Indexer::Fasta;

  my $indexdir = "/tmp";
  my $file = "/home/jdoe/bacteria.nraa";

  my $indexer = BSML::Indexer::Fasta->new($file, $indexdir);

  # Check the health of the various indices for the data file.
  my @check = $indexer->check_indices;

  # Create the indices if necessary...
  if ($check[0] == 1) { $indexer->index_entries };
  if ($check[1] == 1) { $indexer->index_headers };

  # Get the name of the entry index file.
  my $entry_index = $indexer->entry_index;

  # Get the name of the header index file.
  my $header_index = $indexer->header_index;

  my $filepath = $indexer->filepath;
  my $filename = $indexer->filename;

=head1 DESCRIPTION

This class provides methods to check whether index files for Yank data
files are up to date, and if not, to generate them.

=head2 Class and object methods

=over 4

=cut

use strict;
use base qw(BSML::Indexer);
use Data::Dumper;
use File::Basename;
use File::Path;
BEGIN {
use BSML::Logger;
}

my @valid_sources = qw(GB GP PANDA PIR SP);

# Logger instance used by this class
my $logger = BSML::Logger::get_logger(__PACKAGE__);


# Permissions for the index files.
my $mode = 0664;

## internal variables and identifiers
our $VERSION = (qw$Revision: 1.5 $)[1];
our @DEPEND = qw(Data::Dumper Fcntl File::Basename
                 File::Path Log::Log4perl TIGR::Yank::Indexer);

# Eliminate annoying warnings
if ($^W) {  # Boolean for the global warnings flag.
    $VERSION = $VERSION;
    @DEPEND = @DEPEND;
}


=item BSML::Indexer::Fasta->new($file, $indexdir);

B<Description:> This is the constructor which creates an instance
of the class.

B<Parameters:> $file - The datafile (a multi-FASTA file) that is
to be indexed, or checked for updated indices.
 $indexdir - The directory where the indices are to be stored.

 Example:
 my $indexer = BSML::Indexer->new($file, $indexdir);
                                                    

B<Returns:> $self, an instance of the class.

=cut

sub new {
    my ($class, $file, $indexdir, @args) = @_;
    $logger->info("In constructor.");

    my $self = $class->SUPER::new($file, $indexdir, @args);

    return $self;
}


=item $obj->header_index();

B<Description:> This method returns the name of the header index file
for a datafile.

B<Parameters:> None.

B<Returns:> The full path to the header index file.

=cut

sub header_index {
    my ($self, @args) = @_;

    my $filename = $self->filename;
    my $indexdir = $self->indexdir;
    my $header_index = "$indexdir/${filename}.header.cdb";
    return $header_index;
}


=item $obj->check_indices();

B<Description:> Check if the indices for the data source for this search
exist and are up-to-date. 

B<Parameters:> None.

B<Returns:> In scalar context, a reference to a 3 element array containing
flags (0 or 1) as to whether a particular set of indices exists and is
up-to-date.  The first element is indicates whether the "entry" index is okay.
The second and third elements indicate whether the organism and accession
indices are okay. In list context, the array itself is returned.

=cut

sub check_indices {
    my $self = shift;
    $logger->info("In check_indices.");

    # Check that the entry index is created and that it is current
    # by examining the timestamp and comparing it to the datafile's
    # timestamp.
    my ($needs_entry_index, $needs_header_index) = (0, 0);
    my $indexdir = $self->indexdir;

    if (-d $indexdir) {
        my $filepath = $self->filepath;

        # Index filenames
        my $entry_idx = $self->entry_index;
        my $header_idx = $self->header_index;

        my $timestamp_data = (stat($filepath))[9];

        # Check the entry index file.
        if (-e $entry_idx) {
            my $timestamp_entry_idx = (stat($entry_idx))[9];
            if ($timestamp_data > $timestamp_entry_idx) {
                $needs_entry_index = 1;
                $logger->warn("Entry index file is out of date.");
            } elsif ($logger->is_debug()) {
                $logger->debug("Entry index file is okay.");
            }
        } else {
            $needs_entry_index = 1;
            $logger->warn("Entry index file '$entry_idx' is missing.");
        }

        # Check the header index file.
        if (-e $header_idx ) {
            my $timestamp_header_idx = (stat($header_idx))[9];
            if ($timestamp_data > $timestamp_header_idx) {
                $needs_header_index = 1;
                $logger->warn("Header index file is out of date.");
            } elsif ($logger->is_debug()) {
                $logger->debug("Header index file is okay.");
            }
        } else {
            $needs_header_index = 1;
            $logger->warn("Header index file '$header_idx' is missing.");
        }
    } else {
        $logger->warn("Index directory not present. ",
                      "All indices are missing.");
        # Obviously if the directory wasn't there, then the index files
        # also need to be created.
        ($needs_entry_index, $needs_header_index) = (1, 1);
    }

    my @check = ($needs_entry_index, $needs_header_index);
    return wantarray ? @check : \@check;
}


=item $obj->index_headers();

B<Description:> Used to the various header index. For each
data file, the first word in the header (delimited by white space)
is indexed for subsequent searching. For example, in a header of
the form:

  ">abc def..."
       or
  "> abc def..."

the indexed word is "abc".

B<Parameters:> None. 

B<Returns:> None.

=cut

sub index_headers {
    my $self = shift;
    $logger->info("In index_headers.");

    # We don't want this to be interruptible, so ignore CNTROL-C.
    local $SIG{INT} = 'IGNORE';

    my %headers;
    my $filepath = $self->filepath;

    my $indexdir = $self->indexdir;
    $self->_make_indexdir unless (-d $indexdir);

    my $entry_index = 1;

    # Set the record separator to newline followed by a "greater than" sign.
    # We cannot simply use the "greater than" sign because it may be elsewhere
    # inside the header.
    local $/ = "\n>";

    if(-e "$filepath.gz"){
	open(DATA, "<:gzip", "$filepath.gz")  or $logger->logdie("Could not open data file ",
								 "$filepath.gz for reading: $!");
    }
    else{
	open (DATA, "<", $filepath)
	    or $logger->logdie("Could not open data file ",
			       "$filepath for reading: $!");
    }

    while (<DATA>) {
        my $segment = $_;
        my @lines = split('\n', $segment);
        my $header = $lines[0];

        # The first line will not have the > removed from the scan because
        # there was no preceeding newline, so we remove it ourselves.
        $header =~ s/^>// if ($entry_index == 1);

        # Remove leading white space in preparation for
        # the white space split, since we don't want have white space as the
        # indexed term, first_word.
        $header =~ s/^\s+//; 
        my $first_word = (split(/\s+/, $header))[0];
        push (@{ $headers{$first_word} }, $entry_index);
        $entry_index++;   # Increase the FASTA entry count;
    }

    close (DATA)
        or $logger->logdie("Could not close filehandle for $filepath.");

    my $header_index = $self->header_index;

    $logger->debug("Creating header index.");

    my $temporary = "${header_index}.new";

    my $db = CDB_File->new($header_index, $temporary)
        or $logger->logdie("Could not tie hash to $header_index.");

    if ( scalar(keys %headers) ) {
        foreach my $header (sort keys %headers) {
            # The following check is to prevent a dereferencing error on the
            # next line.
            $headers{$header} =
               ( defined($headers{$header}) ) ? $headers{$header} : [];
            my $indices = join(',', sort @{ $headers{$header} });
            $db->insert($header, $indices);
        }
    } else {
        $logger->warn("Nothing to index, index will be empty.");
    }

    $db->finish;

    chmod $mode, $header_index;
}


1;

__END__

=back

=head1 ENVIRONMENT

This module does not set or read any environment variables.

=head1 DIAGNOSTICS

This module has a unit testing script in the Yank installation test
directory.

=head1 BUGS

There are no known bugs at this time. Please contact the ANTware team
at antware@tigr.org, or submit a bug report to bits.antware@tigr.org
if you discover one.

=head1 SEE ALSO

 CDB_File
 Data::Dumper - To print data structures.
 File::Basename - to help determine index filenames and locations.
 File::Path - for the mkpath function, to create missing directories.
 Log::Log4perl - For logging messages and debugging.
 TIGR::Yank::Indexer - Parent class.

=head1 AUTHOR(S)

 The Institute for Genomic Research
 9712 Medical Center Drive
 Rockville, MD 20850

=head1 COPYRIGHT

Copyright (c) 2003-2004, The Institute for Genomic Research. All Rights Reserved.
