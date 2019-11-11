# metylUtil

Scripts to 1) mask C/T polymorphisms in BSseq data, and 2) process genome and annotation to identify transcriptional start sites (TSS), promoter regions, CpG islands and shores

## Pre-requisites

### Software
samtools v1.8+: https://samtools.github.io/
bedtools v2.29.0+: https://bedtools.readthedocs.io/en/latest/
EMBOSS v6.6.6+: http://emboss.sourceforge.net/index.html

### Data

-Reference genome or reference sequences in fasta format (required)
-Gene annotation information in GFF format (only required for creating accessory BED files)
-BSseq results in bedGraph format (similar to that output by Bismark or MethylDackel)
-BED file of C/T and A/G SNP polymorphisms 

## Usage
1) Copy or download the genome in fasta format into the `./02_reference` folder and rename it "genome.fasta"

2) Optional: Copy or download the gene annotation into the `./02_reference` folder and rename it "genes.gff" (can be gzipped)

3) Copy bedGraph files (gzipped) into the `./03_methylation_bedGraphs` folder 

4) Execute scripts of interest and collect results

## Contact

Kyle Wellband - kyle.wellband [at] gmail.com

## License

CC share-alike

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">CT-poly-wgbs</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Kyle Wellband</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
