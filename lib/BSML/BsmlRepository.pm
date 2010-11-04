package BSML::BsmlRepository;

# $Id: BsmlRepository.pm,v 1.24 2006/01/24 21:31:57 angiuoli Exp $

# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.

=head1 NAME

BsmlRepository.pm - A module for managing a BSML repository

=head1 VERSION

This document refers to version $Name: bsml-v2r10b1 $ of frontend.cgi, $Revision: 1.24 $. 
Last modified on $Date: 2006/01/24 21:31:57 $

=head1 SYNOPSIS

=head1 DESCRIPTION

my $bsmlrepository = new BSML::BsmlRepository('PATH'=>$repositorypath);

=over 4

=cut


use strict;
use Data::Dumper;
BEGIN {
use BSML::Logger;
}

=item new

B<Description:> The module constructor.

B<Parameters:> %arg, a hash containing attribute-value pairs to
initialize the object with. Initialization actually occurs in the
private _init method.

my $bsmlrepository = new BSML::BsmlRepository('PATH'=>$repositorypath);

B<Returns:> $self (A BSML::BsmlRepository object).

=cut

sub new {
    my ($class) = shift;
    my $self = bless {}, ref($class) || $class;
    $self->{_logger} = BSML::Logger::get_logger(__PACKAGE__);
    $self->{_BSML_FILE_EXT} = "bsml";
    $self->_init(@_);
    $self->{"_PATH"} = $self->{"_BSML_repository"};
    $self->{_logger}->debug("Setting repository path $self->{_PATH}") if($self->{_logger}->is_debug());
    return $self;
}


=item $obj->_init([%arg])

B<Description:> Initialize object variables

B<Parameters:> %arg, a hash containing attributes to initialize the testing
object with. Keys in %arg will create object attributes with the same name,
but with a prepended underscore.

B<Returns:> None.

=cut

sub _init {
    my $self = shift;
    
    my %arg = @_;
    $self->{_logger}->debug(Dumper(@_)) if($self->{_logger}->is_debug());
    foreach my $key (keys %arg) {
	$self->{_logger}->debug("Parsing argument $key=$arg{$key}") if($self->{_logger}->is_debug());
        $self->{"_$key"} = $arg{$key};
    }
    if(!($self->{"_BSML_repository"})){
	$self->{_logger}->logdie("Required parameter BSML_repository not passed to object constructor");
    }
    if(!(-d $self->{"_BSML_repository"})){
	$self->{_logger}->logdie("Required parameter BSML_repository is not a valid directory");
    }
}

#return the name of the BSML repository
sub get_dirname{
    my $self = shift;
    return $self->{"_PATH"};
}

sub get_assembly_doc{
    my ($self,$asmbl_id) = @_;
    my $file = $self->{_PATH}."/$asmbl_id".".".$self->{_BSML_FILE_EXT};
    if(-e $file){
	return $file;
    }
    else{
	return undef;
    }
}
    

#pull list of assemblies from a glob for now.
#This will work for now because bsml files are named consistently based 
#on the assembly name
sub list_assemblies{
    my $self = shift;
    $self->{_logger}->debug("Listing assemblies from directory $self->{_PATH}") if($self->{_logger}->is_debug());
    opendir BSMLDIR, "$self->{_PATH}" or $self->{_logger}->logdie("Can't read directory $self->{_PATH}");
    my @asmblfiles = grep /\.$self->{_BSML_FILE_EXT}$/, readdir BSMLDIR;
    my @asmbllist;
    foreach my $asmbl (@asmblfiles ){
	$asmbl =~ s/\.$self->{_BSML_FILE_EXT}//;
	$self->{_logger}->debug("Parsed assembly, $asmbl, from filename") if($self->{_logger}->is_debug());
	push @asmbllist, $asmbl;
    }
    return \@asmbllist;    
}

sub list_bsml_files{
    my ($self,$subdir) = shift;
    my $dir = $self->{"_PATH"};
    $self->{_logger}->debug("Listing bsml files from directory $dir") if($self->{_logger}->is_debug());
    opendir BSMLDIR, "$dir" or $self->{_logger}->logdie("Can't read directory $dir");
    my @bsmlfiles = grep /\.$self->{_BSML_FILE_EXT}$/, readdir BSMLDIR;
    my @bsmllist;
    foreach my $bsml (@bsmlfiles ){
	push @bsmllist, "$dir/$bsml";
    }
    return \@bsmllist;    
}
1;
