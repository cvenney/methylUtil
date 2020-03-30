#!/usr/bin/env python

import sys
import re
from Bio import SeqIO

if len(sys.argv) != 4:
    print("Usage: " + sys.argv[0] + " <fasta file> <RE motif (e.g. C^CGG)> <insert size range (e.g. 300-500)>")
    raise SystemExit

infile = sys.argv[1]

if not re.search("\^", sys.argv[2]):
    print("ERROR: Invalid restriction site provided! Must include the postion of the cut site using \"^\".")
    print("Usage: " + sys.argv[0] + " <fasta file> <RE motif (e.g. C^CGG)> <insert size range (e.g. 300-500)>")
    raise SystemExit

motif_5prime, motif_3prime = sys.argv[2].split("^")
motif = ''.join(motif_5prime + motif_3prime)

insert_min, insert_max = (int(x) for x in sys.argv[3].split("-"))

with open(infile, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        end=0
        fragments=re.split(motif, str(record.seq), flags=re.IGNORECASE)
        for frag in range(len(fragments)-1):
            frag_len = len(fragments[frag]) + len(motif_5prime) + len(motif_3prime)
            end = end + frag_len
            if (frag_len >= insert_min) & (frag_len <= insert_max):
                print('\t'.join((record.id, str(end - len(motif_5prime) - len(motif_3prime)), str(end))))
            
