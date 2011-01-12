##################################################################
# Description: UniRef.pm
#
# Data Tranfer Object used by the UnirefDao to return back to the 
# the UnirefDao client.
# --------------
# Author: jgoll 
# Email:  jgoll@jcvi.org
# Date:   Sep 29, 2010  
##################################################################

package UniRef::UniRef;
use strict;
use warnings;

#constructor
sub new() {
	my $class = shift;
   
    my $self = {
    	 _geneSymbol => shift,
        _ecIds => shift,
        _goIds  => shift,
        _cazyIds  => shift,
        _reviewed => shift,
    };
 
    bless $self, $class;
    return $self;
}

#getters
sub getGeneSymbol {
	my $self = shift;
    return $self->{_geneSymbol};
}
sub getEcIds {
	my $self = shift;
    return $self->{_ecIds};
}
sub getGoIds {
	my $self = shift;
    return $self->{_goIds};
}
sub getCazyIds {
	my $self = shift;
    return $self->{_cazyIds};
}
sub isReviewed {
	my $self = shift;
    return $self->{_reviewed};
}

#rerturns true if at least one annotation type has been defined
sub hasAnnotation() {
	my $self = shift;
	
	if($self->{_geneSymbol} eq '' && (scalar @{$self->{_ecIds}}) == 0 && (scalar @{$self->{_goIds}})  == 0 &&
	   (scalar @{$self->{_cazyIds}})  == 0 && $self->{_reviewed} == 0) {
		return 0;
	}
	else {
		return 1;
	}
}

## print object
sub print() {
	my $self = shift(@_);
	print "Gene Symbol:".$self->getGeneSymbol().":\n";
   	print "GO IDs:".join('||',@{$self->getGoIds()}).":\n";
	print "EC IDs:".join('||',@{$self->getEcIds()}).":\n"; 
	print "CAZY IDs:".join('||',@{$self->getCazyIds()}).":\n"; 
    print "Is Reviewed:".$self->isReviewed().":\n";
    print "Has Annotation:".$self->hasAnnotation().":\n";
}
1