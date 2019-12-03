#!/bin/bash
# dml_in_dmr_summary.sh

if [ $# -ne 2 ]
then
    echo "Usage: $0 <DSS_ouput_basename> <DSS_ouput_settings>"
    echo "e.g. tab complete using *dml* or *dmr* file, replace 'dml' or 'dmr' with a space..."
    exit
fi

file=$1
settings=$2

dml=$(gunzip -c "${file}dml${settings}" | awk 'END{print NR -1}')

gunzip -c "${file}dmr${settings}" | 
awk -v dml=$dml '(NR != 1) {
    sum += $5
} END {
    print NR - 1 " DMRs contain " sum " / " dml " (" (sum/dml*100) "%) DMLs"
}'

