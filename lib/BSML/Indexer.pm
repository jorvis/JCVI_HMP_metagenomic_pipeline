package BSML::Indexer;

=head1 NAME

BSML::Indexer;

=head1 VERSION

This document refers to $Revision: 1.8 $ of TIGR::Yank::Indexer.pm.

=head1 SYNOPSIS

  use BSML::Indexer;

  my $indexdir = "/tmp";
  my $file = "/home/jdoe/bacteria.nraa";

  my $indexer = TIGR::Yank::Indexer->new($file, $indexdir);

  my $entry_index = $indexer->entry_index;
  my $org_index = $indexer->organism_index;
  my $acc_index_pir_1 = $indexer->accession_index("PIR", 1);

  my $filepath = $indexer->filepath;
  my $filename = $indexer->filename;

  # Check the health of the various indices for the data file.
  my @check = $indexer->check_indices;

  if ($check[0] == 1) { $indexer->index_entries };
  if ($check[1] == 1) { $indexer->index_organisms };
  if ($check[2] == 1) { $indexer->index_accessions };

=head1 DESCRIPTION

This class provides methods to check whether index files for Yank data
files are up to date, and if not, to generate them.

=head2 Class and object methods

=over 4

=cut

use strict;
use CDB_File 'create';
use Cwd 'abs_path';
use Data::Dumper;
use File::Basename;
use File::Path;
use BSML::Logger;

my @valid_sources = qw(GB GP PANDA PIR SP);

# Logger instance used by this class
my $logger = BSML::Logger::get_logger(__PACKAGE__);


# Permissions for the index files.
my $mode = 0664;

## internal variables and identifiers
our $VERSION = (qw$Revision: 1.8 $)[1];
our @DEPEND = qw(Data::Dumper CDB_File Cwd File::Basename 	 
                  File::Path Log::Log4perl);

# Eliminate annoying warnings
if ($^W) {  # Boolean for the global warnings flag.
    $VERSION = $VERSION;
    @DEPEND = @DEPEND;
}


=item BSML::Indexer->new($file, $indexdir);

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



    my ($package, $filename, $line) = caller;

#    print STDERR "Indexer.pm line 99: package '$package' filename '$filename' line '$line'\n";

    $logger->debug("Datafile: $file; Indexdir: $indexdir.");

    my $self = bless {}, ref($class) || $class;

    # Call the init method to initialize the lookup hash.
    $self->_init($file, $indexdir);

    return $self;
}


=item $obj->_init($file, $indexdir);

B<Description:> This is the object initializer, a private method used to set
object attributes.

B<Parameters:> Both the $file, and $indexdir parameters described for the
constructor "new" are required.
B<Returns:> None.

=cut

sub _init {
    my ($self, $file, $indexdir) = @_;
    $logger->info("In _init.");
    $logger->logdie("\$file and \$indexdir parameters are required.")        unless ((defined $file) && defined($indexdir));

    my $filename = basename($file);

    $logger->debug("Datafile: $file; Indexdir: $indexdir.");
    $self->{filepath} = $file;
    $self->{indexdir} = $indexdir || abs_path(dirname($file));
    $self->{filename} = $filename;
}


=item $obj->accession_index($source, $position);

B<Description:> This method returns the name of the accession index file
for a particular source type and position, which are specified when calling
the method.

B<Parameters:> $source - a valid source type, which must be one of the
following: GB, GP, PANDA, PIR, SP.
 $position - Must be either 1 or 2, which are the only valid positions
 in an accession. The positions are separated by the pipe symbol "|".

B<Returns:> The full path to the accession index file.

=cut

sub accession_index {
    my ($self, $source, $position, @args) = @_;
    unless ($position == 1 || $position == 2) {
        $logger->logdie("Incorrect position parameter. Must be 1 or 2.");
    }

    # Should validate $source against @valid_sources here.

    my $filename = $self->filename;
    my $indexdir = $self->indexdir;
    my $acc_index = "$indexdir/${filename}.acc.${source}.${position}.cdb";
    return $acc_index;
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
    my ($needs_entry_index, $needs_org_index, $needs_acc_index) = (0, 0, 0);
    my $indexdir = $self->indexdir;

    if (-d $indexdir) {
        my $filepath = $self->filepath;

        # Index filenames
        my $entry_idx = $self->entry_index;
        my $org_idx = $self->organism_index;

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
            $logger->warn("Entry index file is missing.");
        }

        # Check the organism index file.
        if (-e $org_idx ) {
            my $timestamp_org_idx = (stat($org_idx))[9];
            if ($timestamp_data > $timestamp_org_idx) {
                $needs_org_index = 1;
                $logger->warn("Organism index file is out of date.");
            } elsif ($logger->is_debug()) {
                $logger->debug("Organism index file is okay.");
            }
        } else {
            $needs_org_index = 1;
            $logger->warn("Organism index file is missing.");
        }

        # There are numerous accession indices, so we have a separate method
        # to do the checking.
        $needs_acc_index = $self->_check_accession_indices;
    } else {
        $logger->warn("Index directory not present. ",
                      "So all indices are missing...");
        # Obviously if the directory wasn't there, then the index files
        # also need to be created.
        ($needs_entry_index, $needs_org_index, $needs_acc_index) = (1, 1, 1);
    }

    my @check = ($needs_entry_index, $needs_org_index, $needs_acc_index);
    return wantarray ? @check : \@check;
}


=item $obj->filepath();

B<Description:> An accessor used to retrieve the full path to the filename
used as the data source for the search.

B<Parameters:> None.

B<Returns:> A scalar holding the path to the datafile.

=cut

sub filepath { $_[0]->{filepath} }


=item $obj->entry_index();

B<Description:> Get the path to the entry index file.

B<Parameters:> None.

B<Returns:> None.

=cut

sub entry_index {
    my $self = shift;
    $logger->info("In entry_index.");
    unless (exists ($self->{entry_index})) {
        my $indexdir = $self->indexdir;
        my $filename = $self->filename;
        my $entry_index = "$indexdir/${filename}.entry.cdb";
        $self->{entry_index} = $entry_index;
    }
    return $self->{entry_index};
}


=item $obj->filename();

B<Description:> An accessor used to retrieve the name of the file
to be used as the data source for the search. This is simply the
basename of the data file as returned by File::Basename's basename
function.

B<Parameters:> None.

B<Returns:> A scalar holding the name of the data file.

=cut

sub filename { $_[0]->{filename} };


=item $obj->indexdir();

B<Description:> Returns the directory where the various index files are
to be stored.

B<Parameters:> None.

B<Returns:> The scalar path to the directory.

=cut

sub indexdir { $_[0]->{indexdir} }


=item $obj->index_accessions();

B<Description:> Used to create the various accession indices. For each
data file, there are 2 indices created for each source type. The first is 
to index entry index number with the accession found in the
first position (separated by a pipe), and the second is for the accession
found in the second position. 

B<Parameters:> None. 

B<Returns:> None.

=cut

sub index_accessions {
    my $self = shift;
    $logger->info("In index_accessions.");

    # We don't want this to be interruptible, so ignore CNTROL-C.
    local $SIG{INT} = 'IGNORE';

    my %accessions;

    my $filepath = $self->filepath;
    my $accession_regex = qr/^([A-Za-z]{2,})\|([\w\.]+)\|([\w\.]+)/;
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
	    or $logger->logdie("Could not open data file '$filepath' for reading: $!");
    }
    while (<DATA>) {
        my $segment = $_;
        my @lines = split('\n', $segment);
        my $header = $lines[0];

        # Since the number of entries could be very high, only log
        # every 100,000th one encountered.
        $logger->debug("Reading entry index number: $entry_index.")
            if (($entry_index % 100000) == 0);

        # The first line will not have the > removed from the scan because
        # there was no preceeding newline, so we remove it ourselves.
        $header =~ s/^>// if ($entry_index == 1);

        # Headers use the following delimiter: ^|^
        my @parts = split (/\^\|\^/, $header);
        my $source;
        foreach my $part (@parts) {
            if ($part =~ m/$accession_regex/) {
                # $1 = Source (SP, GB, GI, PIR, PANDA)
                # $2 = unique numeric ID
                # $3 = accession identifier (dependent on Source).

                # Ensure keys are upper cased, as later search terms
                # from the user will be uppercased as well.
                $source = uc($1);

                my ($first, $second) = my ($first_v, $second_v) = ($2, $3);

                # We are hashing both $2 and $3 because we don't know which
                # one the user may ask for. We therefore, hash both. In
                # addition, we have to hash (or index), the versioned, and
                # unversioned idents, because specifying the version is
                # optional, and if unspecified, the user should get all the
                # results. This approach makes the hash, and resulting index file,
                # much larger.

                # For the first position index.
                if ($first_v) {
                    $first =~ s/\..*$//;
                    push (@{ $accessions{$source}->{1}->{$first_v} },
                        $entry_index);
                    push (@{ $accessions{$source}->{1}->{$first} },
                        $entry_index) if ($first ne $first_v);
                }
                # For the second position index    
                if ($second_v) {
                    $second =~ s/\..*$//;
                    push (@{ $accessions{$source}->{2}->{$second_v} },
                        $entry_index);
                    push (@{ $accessions{$source}->{2}->{$second} },
                        $entry_index) if ($second ne $second_v);
                }
            } else {
                $logger->info("Index number $entry_index has a header that does",
                              "not match the selection criteria.");
            }
        }
        $entry_index++;   # Increase the FASTA entry count;
    }

    close (DATA)
        or $logger->logdie("Could not close filehandle for $filepath.");

    # $identifier represents the positions in the whole identifier:
    # $source|1|2.
    my ($identifier, $acc_index);
    my $indexdir = $self->indexdir;
    my $filename = $self->filename;

    $self->_make_indexdir unless (-d $indexdir);

    foreach $identifier qw( 1 2 ) {
        $logger->debug("Creating indices for accession $identifier.");

        # Not all the elements of valid_sources may be listed as actual keys
        # in the %accessions hash. But those that are not there will have dbm
        # files created anyway, even if empty, so that future index checks
        # will not trip everytime (if one is missing all are recreated, which
        # is expensive).
        foreach my $source (@valid_sources) {
            $logger->debug("Creating indices for accession $identifier, ",
                           "source $source.");

            $acc_index = $self->accession_index($source, $identifier);
            my $temporary = "${acc_index}.new";

            my %subhash = (exists $accessions{$source}->{$identifier})
                ?  %{ $accessions{$source}->{$identifier} } : ();

            if ( scalar(keys %subhash) ) {
                my $db = CDB_File->new($acc_index, $temporary)
                    or $logger->logdie("Could not tie hash to $acc_index.");
                foreach my $key (keys %subhash) {
                    my $value = join(',', @{ $subhash{$key} } );
                    $db->insert($key, $value);
                }
                $db->finish;
            } else {
                $logger->debug("Subhash for $source was empty. ",
                               "Tie-ing emtpy hash.");
                # If %subhash is empty, this should make an empty DBM file.
                my %accessions_idx = ();
                create %accessions_idx, $acc_index, $temporary;
            }

            chmod $mode, $acc_index;
            # Since we're done with this branch of the hash, let's try to
            # conserve memory by deleting it.
            delete $accessions{$source}->{$identifier};
        }
    }
}


=item $obj->index_entries();

B<Description:> Method to create the entry index file, which holds the
information on how to access the entries in a FASTA file by entry index
number. The index contains information on the starting position and length
of each entry.

B<Parameters:> None.

B<Returns:> None.

=cut

sub index_entries {
    my $self = shift;

    $logger->debug("In index_entries.") if $logger->is_debug;

    # We don't want this to be interruptible, so ignore CONTROL-C.
    local $SIG{INT} = 'IGNORE';

    # Set the record separator to newline followed by a "greater than" sign.
    # We cannot simply use the "greater than" sign because it may be elsewhere
    # inside the header.
    local $/ = "\n>";

    my $filepath = $self->filepath;
    my $entry_index = 1;
    my $start = 0;
    my $total_length = 0;
    my $end = 0;
    my $length;

    $logger->debug("filepath '$filepath'") if $logger->is_debug;



    my $indexdir = $self->indexdir;
    $self->_make_indexdir unless (-d $indexdir);

    my $entry_index_file = $self->entry_index;

    # Let's get rid of it, if there's already an outdated copy there.
    unlink $entry_index_file if (-e $entry_index_file);

    my $temporary = "${entry_index_file}.new";

    my $db;

    eval {
	$db = CDB_File->new($entry_index_file, $temporary) || 
	$logger->logdie("Cannot create entry index file ",
			"'$entry_index_file': $!");
    };
    if ($@){
	$logger->logdie("Caught exception while attempting to instantiate a CDB_File ".
	"with following parameters: entry_index_file '$entry_index_file' temporary '$temporary'.  ".
	"The eval exception was: $!");
    }

    # The index number starts at 1, and we will index the start position
    # in the file, and the length of the string to the index number.
    # The start and length are sufficient to then use syseek to extract
    # the entry from the file.

    if(-e "$filepath.gz"){
	open (DATA, "<", "$filepath.gz")
	    or $logger->logdie("Could not open data file '$filepath.gz' for reading.: $!");

    }
    else{
	open (DATA, "<", $filepath)
	    or $logger->logdie("Could not open data file ",
			       "$filepath for reading: $!");
    }
    while (<DATA>) {
        my $segment = $_;
        $length = ($entry_index == 1) ?
            length($segment) - 1 : length($segment);
        $start = $total_length;
        $total_length += $length;
        $end = $total_length;
        # Insert the entry into the index.
        $db->insert($entry_index, "$start,$length");

        $entry_index++;   # Increase the FASTA entry count;
    }

    close (DATA)
        or $logger->logdie("Could not close filehandle for $filepath.");

    $db->finish;
    chmod $mode, $entry_index_file;
}


=item $obj->index_organisms();

B<Description:> Method to create the organism index file.

B<Parameters:> None. 

B<Returns:> None.

=cut

sub index_organisms {
    my $self = shift;
    $logger->info("In organism_index.");

    # We don't want this to be interruptible, so ignore CNTROL-C.
    local $SIG{INT} = 'IGNORE';

    my $filepath = $self->filepath;
    my $indexdir = $self->indexdir;
    my $org_index = $self->organism_index;

    # The organisms hash keeps track of which indices each organism
    # appears in. The key is the organism name, and the the value is
    # an array reference containing entry indices.
    my %organisms;

    my $entry_index = 1;

    # Precompile the regular expression for greater speed.
    my $organism_regex = qr/\{(.+)\}/;

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

        my @parts = split (/\^\|\^/, $header);
        my $source;
        foreach my $part (@parts) {
            if ($part =~ m/$organism_regex/) {
                my $organism = $1;
                # Remove excess whitespace.
                $organism =~ s/\s{2,}/ /g;
                # Remove leading whitespace.
                $organism =~ s/^\s+//;
                push (@{ $organisms{$organism} }, $entry_index);
            }
        }
        $entry_index++;   # Increase the FASTA entry count;
    }

    close (DATA)
        or $logger->logdie("Could not close filehandle for $filepath.");

    $self->_make_indexdir unless (-d $indexdir);

    # We may have a problem if multiple users try to create the index at the
    # same time. Needs investigation.
    my $temporary = "${org_index}.new";
    my $db = CDB_File->new($org_index, $temporary)
        or $logger->logdie("Could not tie organisms index filehandle.");

    foreach my $org (sort keys %organisms) {
        # The following check is to prevent a dereferencing error on the
        # next line.
        $organisms{$org} =
            ( defined($organisms{$org}) ) ? $organisms{$org} : [];
        my $indices = join(',', sort by_number @{ $organisms{$org} });
        $db->insert($org, $indices);
    }
    $db->finish;
    chmod $mode, $org_index;
}


=item $obj->organism_index();

B<Description:> Get the organism index file.

B<Parameters:> None.

B<Returns:> The filename of the organism index.

=cut

sub organism_index {
    my $self = shift;
    $logger->debug("In organism_index.");

    unless (exists $self->{organism_index}) {
        my $filename = $self->filename;
        my $indexdir = $self->indexdir;
        my $org_index = "$indexdir/${filename}.org.cdb";
        $self->{organism_index} = $org_index;
    }
    return $self->{organism_index};
}


=item $obj->_check_accession_indices();

B<Description:> This private method is used to examine for the existence and
timeliness of the various accession indices for the data file to be used in a
particular search. If any of the file are missing or out of date, a true value
(1) is returned to indicate that the files should be rebuilt.

B<Parameters:> None.

B<Returns:> $needs_acc_index, a flag denoting whether the accession indices
need to be recreated.

=cut

sub _check_accession_indices {
    my $self = shift;
    $logger->info("In _check_accession_indices.");

    # Begin by assuming we don't need to re-index, and if we detect something
    # out of date or missing, set to 1, and return.
    my $needs_acc_index = 0;

    my $filepath = $self->filepath;
    my $filename = $self->filename;
    my $indexdir = $self->indexdir;
    my $timestamp_data = (stat($filepath))[9];
    
    OUTER: foreach my $source (@valid_sources) {
        for (my $i=1; $i<=2; $i++) {
            my $acc_idx = $self->accession_index($source, $i);
            if (-e $acc_idx) {
                my $timestamp_acc_idx = (stat($acc_idx))[9];
                if ($timestamp_data >  $timestamp_acc_idx) {
                    $needs_acc_index = 1;
                    $logger->warn("Accession index file for $source id $i ",
                                  "is out of date.");
                    last OUTER;
                }
            } else {
                $logger->warn("Accession index for $source id $i missing.");
                $needs_acc_index = 1;
                last OUTER;
            }    
        }
    }
    if ($needs_acc_index == 1) {
        $logger->debug("Accession indices need to be regenerated.");
    } elsif ($needs_acc_index == 0) {
        $logger->debug("Accession indices appear to be updated.");
    }
    return $needs_acc_index;
}


=item $obj->_make_indexdir;

B<Description:> This private method is used to create the directory that
will hold the various index files if does not exist at the time the
indexing is initiated.

B<Parameters:> None.

B<Returns:> None.

=cut

sub _make_indexdir {
    my $self = shift;
    $logger->info("In _make_indexdir.");

    my $indexdir = $self->indexdir;
    unless (-d $indexdir) {
        $logger->warn("Index directory not present. Creating it...");
        mkpath($indexdir, 0, 0777)
            or $logger->logdie("Could not create index directory: $!");
    }
}

sub by_number {
    my ($a, $b);
    $a <=> $b;
}


1;

__END__

=back

=head1 ENVIRONMENT

This module does not read or set any environment variables.

=head1 DIAGNOSTICS

This module has a unit testing script in the Yank installation test
directory.

=head1 BUGS

There are no known bugs at this time. Please contact the Antware team
at antware@tigr.org, or submit a bug report to bits.antware@tigr.org
if you discover one.

=head1 SEE ALSO

 CDB_File - The platform independent database file format.
 Data::Dumper - To print data structures.
 File::Basename - To help determine index filenames and locations.
 File::Path - For the mkpath function, to create missing directories.
 Log::Log4perl - For logging messages and debugging.

=head1 AUTHOR(S)

 The Institute for Genomic Research
 9712 Medical Center Drive
 Rockville, MD 20850

=head1 COPYRIGHT

Copyright (c) 2003, The Institute for Genomic Research. All Rights Reserved.
