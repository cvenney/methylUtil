#!/usr/bin/env python3
# split_bedGraph_by_chr.py

import sys
import gzip
import re

if len(sys.argv) != 3:
    print("Usage: " + sys.argv[0] + " <*.bedGraph.gz> <chr_prefix [e.g. NC_0273]>")
    sys.exit()

infile = sys.argv[1]
chr_prefix = sys.argv[2] + ".*"
print(infile)
outbase = infile.replace(".bedGraph.gz", "_")
cchr = ""
outfile = ""
chrfile = open("02_reference/chrs.txt", "w")


with gzip.open(infile, "rt") as file:
    header = file.readline()
    for line in file:
        lchr = line.split()[0]
        if lchr != cchr and re.match(chr_prefix, lchr):
            chrfile.write(lchr)
            if outfile != "":
                outfile.close()
            cchr = lchr
            outname = outbase + lchr + ".bedGraph.gz"
            outfile = gzip.open(outname, "at")
            outfile.write(header)
            outfile.write(line)
        elif lchr != cchr and re.match(chr_prefix, cchr) and (re.match("NW_.*", lchr) or re.match("NC_00.*", lchr)):
            if outfile != "":
                outfile.close()
            cchr = lchr
            outname = outbase + "contigs.bedGraph.gz"
            outfile = gzip.open(outname, "at")
            outfile.write(header)
            outfile.write(line)
        else:
            outfile.write(line)
outfile.close()
chrfile.write("contigs")
chrfile.close()
