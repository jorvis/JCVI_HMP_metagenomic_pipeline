#!/usr/local/bin/perl

## contains class variables that define recognized annotation types

package CAMERA::AnnotationData::PolypeptideDataTypes;

use strict;
use warnings;

our @type_list = (
                'TIGRFAM::FullLength::Equivalog',
                'TIGRFAM::FullLength::HypotheticalEquivalog',
                'TIGRFAM::FullLength::Exception',
                'TIGRFAM::FullLength::Subfamily',
                'TIGRFAM::FullLength::Superfamily',
                'TIGRFAM::FullLength::EquivalogDomain',
                'TIGRFAM::FullLength::HypotheticalEquivalogDomain',
                'TIGRFAM::FullLength::SubfamilyDomain',
                'TIGRFAM::FullLength::Domain',
                'TIGRFAM::FRAG::Equivalog',
                'TIGRFAM::FRAG::HypotheticalEquivalog',
                'TIGRFAM::FRAG::Exception',
                'TIGRFAM::FRAG::Subfamily',
                'TIGRFAM::FRAG::Superfamily',
                'TIGRFAM::FRAG::EquivalogDomain',
                'TIGRFAM::FRAG::HypotheticalEquivalogDomain',
                'TIGRFAM::FRAG::SubfamilyDomain',
                'TIGRFAM::FRAG::Domain',
                'PFAM::FullLength::Equivalog',
                'PFAM::FullLength::HypotheticalEquivalog',
                'PFAM::FullLength::Subfamily',
                'PFAM::FullLength::Superfamily',
                'PFAM::FullLength::EquivalogDomain',
                'PFAM::FullLength::HypotheticalEquivalogDomain',
                'PFAM::FullLength::SubfamilyDomain',
                'PFAM::FullLength::Domain',
                'PFAM::FullLength::Uncategorized',
                'PFAM::FRAG::Equivalog',
                'PFAM::FRAG::HypotheticalEquivalog',
                'PFAM::FRAG::Subfamily',
                'PFAM::FRAG::Superfamily',
                'PFAM::FRAG::EquivalogDomain',
                'PFAM::FRAG::HypotheticalEquivalogDomain',
                'PFAM::FRAG::SubfamilyDomain',
                'PFAM::FRAG::Domain',
                'PFAM::FRAG::Uncategorized',
#                'PandaBLASTP::Characterized',
#                'PandaBLASTP::HighConfidence',
#                'PandaBLASTP::Putative',
#                'PandaBLASTP::ConservedDomain',
                'UnirefBLASTP::Reviewed',
                'UnirefBLASTP::HighConfidence',
                'UnirefBLASTP::Putative',
                'UnirefBLASTP::ConservedDomain',    
                'UnirefBLASTP::LowConfidence',                
                'PublicSourceAnnotation',
                'TMHMM',
                'LipoproteinMotif',
                'PRIAM',
                'Hypothetical',
                 );

our @attribute_list = (
                        'common_name',
                        'gene_symbol',
                        'GO',
                        'EC',
                        'CAZY',
                        'TAXON'
                      );                   
1;
