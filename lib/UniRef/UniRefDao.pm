#!/usr/local/bin/perl
package UniRef::UniRefDao;
use UniRef::UniRef;
use strict;
use warnings;
use DBI();


#constructor
sub new {
	my ($class,$sqliteDb) = @_;
	my $self = bless({}, $class);
	
	
	
	## connect to mysql database
	$self->{DBH} = DBI->connect( "dbi:SQLite:$sqliteDb", "", "", {PrintError=>1,RaiseError=>1,AutoCommit=>0} );
	if ( !defined $self->{DBH} ) {
		die "could not connect to sqlite database $sqliteDb: " . $DBI::errstr;
	}	
	
#	$self->{DBH} = DBI->connect("DBI:mysql:ifx_commonDB;host=mysql51-lan-pro.jcvi.org",
#		"access", "access", { 'RaiseError' => 1 });		
	return $self;
}

sub close {
	my ($self) = @_;
	$self->{DBH}->disconnect;
}

#generate a mock entry
sub getMockEntry() {
	my ($self,$uniprotId) = @_;
	
	my $unirefEntry = undef;
	
	if($uniprotId eq 'C0QRJ2') {
		my @goIds = ('GO:0004650','GO:0005975','GO:0007047');
		my @ecIds = ('2.3.1.41','4.2.1.61');
		my @cazyIds = ();
		my @taxonIds = (123214);
		$unirefEntry = new UniRef::UniRef('C0QRJ2','Protein translocase subunit secA','secA','Persephonella marina',\@ecIds,\@goIds,\@cazyIds,\@taxonIds,125,1) ;
	}
	elsif($uniprotId eq 'D3SNW0') {
		my @goIds = ('GO:0005737','GO:0005886','GO:0005524');
		my @ecIds = ('2.3.1.38');
		my @cazyIds = ();
		my @taxonIds = (638303);
		$unirefEntry = new UniRef::UniRef('A3HWR0','Protein translocase subunit secA','secA','Thermocrinis albus (strain DSM 14484 / JCM 11386 / HI 11/12)',\@ecIds,\@goIds,\@cazyIds,\@taxonIds,125,1) ;
	}
	elsif($uniprotId eq 'C3Q5L0') {
		my @goIds = ('GO:0004650','GO:0005975','GO:0007047');
		my @ecIds = ();
		my @cazyIds = ();
		my @taxonIds = (457395);
		$unirefEntry = new UniRef::UniRef('C3Q5L0','Glycoside hydrolase family 28 protein','','Bacteroides sp. 9_1_42FAA',\@ecIds,\@goIds,\@cazyIds,\@taxonIds,125,1) ;		
	}
	else {
		my @goIds = ('GO:0004650','GO:0005975','GO:0007047');
		my @ecIds = ();
		my @cazyIds = ();
		my @taxonIds = (457395);
		$unirefEntry = new UniRef::UniRef('C3Q5L0','Glycoside hydrolase family 28 protein','','Bacteroides sp. 9_1_42FAA',\@ecIds,\@goIds,\@cazyIds,\@taxonIds,125,1) ;		
	}
	return $unirefEntry;
}

## fetch uniprot entry
sub get() {
	my ($self, $uniprotId) = @_;

	my $geneSymbol = '';
	my %goIds = ();
	my %ecIds = ();
	my %cazyIds = ();
	my $isReviewed = 0;

	#my $uniprotQuery ="SELECT external_db_name,external_db_identifier FROM uniref_clusters_map where uniref_identifier=?";
	my $uniprotQuery ="SELECT key,value FROM annotation where accession=?";
	
	## execute query
	my $sth = $self->{DBH}->prepare($uniprotQuery);
	$sth->execute($uniprotId);

	my (
		$type,$value
	);

	$sth->bind_col(1, \$type);
	$sth->bind_col(2, \$value);
	
	while ($sth->fetch) {
		if($type eq 'GN') {$geneSymbol = $value};
		if($type eq 'GO') {$goIds{$value} = undef};		
		## union of BRENDA (EC_DR) and Uniprot (EC_DE) EC assignments
		if($type eq 'EC_DR') {$ecIds{$value} = undef};
		if($type eq 'EC_DE') {$ecIds{$value} = undef};
		if($type eq 'CAZY') {$cazyIds{$value} = undef};	
		if($type eq 'REVIEWED') {if($value eq 'TRUE'){$isReviewed = 1;}};
	}
	
	my @uniqueGoIds = keys %goIds;
	my @uniqueEcIds = keys %ecIds;
	my @uniqueCazyIds = keys %cazyIds;
	
	return new UniRef::UniRef($geneSymbol,\@uniqueEcIds,\@uniqueGoIds,\@uniqueCazyIds,$isReviewed);		
}
1;