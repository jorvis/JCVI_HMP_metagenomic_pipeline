#!/usr/local/bin/perl

package CAMERA::AnnotationRules::Util;

use strict;
use warnings;
use Carp;

sub ec_collapse {
    my ($ec_arr_ref) = @_;

    unless (ref($ec_arr_ref) eq 'ARRAY') {
        confess "ec_collapse called with incorrect parameters";
    }
    
    my $ec = {};
   
    foreach my $ec_num(@{$ec_arr_ref}) {
        ## just in case
        chomp $ec_num;

        my @ec_pos = split(/\./, $ec_num);

        ## clean up parent level ec numbers
        for (my $i = 0; $i < 4; $i++) { 
            unless ($ec_pos[$i] && $ec_pos[$i] =~ /^\d+$/) {
                $ec_pos[$i] = '-';
            }
        }
                
        if ($ec_pos[0] ne '-' && ! exists $ec->{$ec_pos[0]}) {
            $ec->{$ec_pos[0]} = {};
        }
        if ($ec_pos[1] ne '-' && ! exists $ec->{$ec_pos[0]}->{$ec_pos[1]}) {
            $ec->{$ec_pos[0]}->{$ec_pos[1]} = {};
        }
        if ($ec_pos[2] ne '-' && ! exists $ec->{$ec_pos[0]}->{$ec_pos[1]}->{$ec_pos[2]}) {
            $ec->{$ec_pos[0]}->{$ec_pos[1]}->{$ec_pos[2]} = {};
        }
        if ($ec_pos[3] ne '-' && ! exists $ec->{$ec_pos[0]}->{$ec_pos[1]}->{$ec_pos[2]}->{$ec_pos[3]}) {
            $ec->{$ec_pos[0]}->{$ec_pos[1]}->{$ec_pos[2]}->{$ec_pos[3]} = {};
        }
    }
   
    my @ec_results;
    my @keys1 = keys(%{$ec});
    my $ec_string = '';
    foreach my $key1(@keys1) {
        $ec_string = "$key1";
        my @keys2 = keys(%{$ec->{$key1}});
        unless (@keys2) {
            $ec_string .= '.-.-.-';
        }
        foreach my $key2(@keys2) {
            $ec_string .= ".$key2";
            my @keys3 = keys(%{$ec->{$key1}->{$key2}});
            unless (@keys3) {
                $ec_string .= '.-.-';
            }
            foreach my $key3(@keys3) {
                $ec_string .= ".$key3";
                my @keys4 = keys(%{$ec->{$key1}->{$key2}->{$key3}});
                unless (@keys4) {
                    $ec_string .= '.-';
                }
                foreach my $key4(@keys4) {
                    my $ec_string .= ".$key4";
                }
            }
        }
        push (@ec_results, $ec_string); 
    }
    
    return \@ec_results;
}

1;
