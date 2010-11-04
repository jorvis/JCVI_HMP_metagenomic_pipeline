#!/usr/local/bin/perl

package CAMERA::Parser::Hypothetical;

use strict;
use warnings;

use CAMERA::Polypeptide;
use CAMERA::PolypeptideSet;
use CAMERA::AnnotationData::Polypeptide;
use Carp;

my $to_text;
my $text_fh;

sub new {
    my ($class, %args) = @_;

    my $self = {};

    ## pass in a reference to a CAMERA::PolypeptideSet object
    ## can be null if no other annotation has been parsed yet
    if (ref($args{'polypeptides'}) eq 'CAMERA::PolypeptideSet') {
        $self->{'polypeptides'} = $args{'polypeptides'};
    } else {
        $self->{'polypeptides'} = new CAMERA::PolypeptideSet();
    }

    return bless $self, $class;
}

sub parse {
    my ($self, $file, $text_mode, $text_fh) = @_;

    if ($text_mode) {
        $to_text = 1;
    }

    my $infh = _open_file_read($file);
 
    my $result = '';
    my $counter = 0; 
    while (<$infh>) {
        if (/>(\S+)/) {
            my $pep_id = $1;

            my $pep = undef;
            if ($self->{'polypeptides'}->exists($pep_id)) {
                $pep = $self->{'polypeptides'}->get($pep_id);
            } else {
                $pep = new CAMERA::Polypeptide('id' => $pep_id);
                $self->{'polypeptides'}->add($pep);
            }
            
            my $annot = new CAMERA::AnnotationData::Polypeptide(
                            'id'           => $pep_id,
                            'type'         => 'Hypothetical',
                            'source'       => 'CAMERA',
                            'rank'         => 1,
                                                           );
            
            #$annot->add_attribute('common_name', 'protein of unassigned function');
            $annot->add_attribute('common_name', 'hypothetical protein');

            if ($to_text) {
                if ($text_fh) {
                    print $text_fh $annot->to_string();
                } else {
                    $result .= $annot->to_string();
                }
            } else {
                $pep->add_annotation($annot);
            }

        } else {
            next;
        }

    }
    return $result;
}

sub _open_file_read {
    my ($file) = @_;

    open (my $infh, $file) || die "Failed to open file '$file' for reading: $!";

    return $infh;
}

1;
