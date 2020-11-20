#!/bin/bash

inBed="$1"

cat 05_bed_files/repeatmasker.bed |
grep -v "Simple_repeat" |
grep -v "Low_complex" |
grep -v "tRNA=" |
grep -v "rRNA=" |
grep -v "snRNA=" |
grep -v "scRNA=" |
grep -v "Satellite=" > 05_bed_files/TEs.bed

grep -v "Unknown=" |

bedtools intersect -a $inBed -b 05_bed_files/TEs.bed -wb > ${inBed%.*}_rm_overlap.bed

