#!/bin/bash

inNMI="$1"

bedtools closest -a 05_bed_files/tss.bed -b ${inNMI} -d |
awk '{if ($13!=0) {
        print $4, "NA"
    } else if ($13==0) {
        if ($6 == "+") {
            print $4, $2 - $8
        } else if ($6 == "-") {
            print $4, $9 - $3
        }
    }
}' > 05_bed_files/upstream_NMI_widths_$(basename ${inNMI%.*}).txt

