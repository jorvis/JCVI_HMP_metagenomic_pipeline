##################################################################
# Description: PNU.pm
# --------------
# Author: jgoll 
# Email:  jgoll@jcvi.org
# Date:   Sep 29, 2010  
##################################################################

package CAMERA::Parser::SynonymsDao;

use strict;
use warnings;

#constructor
sub new() {
	my ($class,$fileName) = @_;
    my $synonymsHash = undef; 	
      	
    if($fileName)	 {
	   	open FILE, $fileName or die "could not open $fileName"; 
		
	   	while(<FILE>) {
	   		chomp $_;
	   		my ($badName,$goodName) = split('\t',$_);
	   		$synonymsHash->{lc($badName)} = $goodName;
	   	}
	   
	   	close FILE;
    }
    else {
    	$synonymsHash = ();
    }
   	
    my $self = {
        _synonymsHash => $synonymsHash,
    };
    
    bless $self, $class;
    return $self;
}

#getters
sub rename {
	my ($self,$name) = @_;
	
   if(exists $self->{_synonymsHash}->{lc($name)}) {
   		$name = $self->{_synonymsHash}->{lc($name)};
   }
   return $name;
}
1

