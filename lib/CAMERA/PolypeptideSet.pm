#!/usr/local/bin/perl

package CAMERA::PolypeptideSet;

use strict;
use warnings;
use Carp;

my $set = {};

sub new {
    my ($class, %args) = @_;

    ## allow a hash ref to be passed in and used for storing polypeptide set
    if (ref($args{'HASH'}) eq 'HASH' || ref($args{'HASH'}) eq 'DBM::Deep::Hash') {
        $set = $args{'HASH'};
    }
    
    return bless {'_set' => $set}, $class;
}

sub exists {
    my ($self, $id) = @_;

    unless ($id) {
        confess("exists called with empty string");
    }
    
    if (exists($self->{'_set'}->{$id})) {
        return 1;
    } else {
        return 0;
    }
}

sub add {
    my ($self, $pep) = @_;
   
    unless (ref($pep) eq 'CAMERA::Polypeptide') {
        return 0;
    } else {
        my $id = $pep->{'id'};
        $self->{'_set'}->{$id} = $pep;
        return 1;
    }
}

sub get {
    my ($self, $id) = @_;

    unless ($id) {
        confess "get called with empty id string";
    }
    
    return $self->{'_set'}->{$id};
}

sub get_keys {
    my ($self) = @_;

    return keys(%{$self->{'_set'}});
}

1;
