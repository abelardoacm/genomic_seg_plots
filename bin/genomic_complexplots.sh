#!/bin/bash
# sample use ./genomic_complexplot.sh $1 12 1.9 2.1
#
#
#
#
perl Genbank_to_peptide_locations.pl $1.gb
perl Gb_2_faa.pl $1.gb
./Just_a_seg_envelope.sh $1 $2 $3 $4
Rscript Complexity_plots.R $1

