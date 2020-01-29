# methylUtil

Scripts to:  
1) mask C/T polymorphisms in BSseq data 
2) process genome and annotation to identify:
   - transcriptional start sites (TSS)
   - promoter regions
   - CpG islands and shores
3) test for differential methylation using R packages "_DSS_" and "_dmrseq_"

## Pre-requisites

#### Software

bedtools: https://bedtools.readthedocs.io/en/latest/  
R: https://www.r-project.org/  
DSS: http://bioconductor.org/packages/release/bioc/html/DSS.html  

_Optional:_  
dmrseq: http://bioconductor.org/packages/release/bioc/html/DSS.html  
EMBOSS: http://emboss.sourceforge.net/index.html  

#### Data

- Reference genome or reference sequences in fasta format (required)
- Gene annotation information in GFF format (only required for creating accessory BED files)
- BSseq results in bedGraph format (similar to that output by Bismark or MethylDackel)
- BED file of C/T and A/G SNP polymorphisms 

## Usage

1) Copy or download the genome in fasta format into the `./02_reference` folder and rename it "genome.fasta"

2) Optional: Copy or download the gene annotation into the `./02_reference` folder and rename it "genes.gff" (can be gzipped)

3) a) Copy bedGraph files (gzipped) into the `./03_raw_bedGraphs` folder and execute `./01_scripts/01_remove_CT_SNPs_from_bedGraphs.sh`

_OR_

3) b) Place analysis ready bedGraph files (gzipped)into the `./04_filtered_bedGraphs` folder 

4) Create sample information file (a script to make a crude file based on the available bedGraphs is found in the `./01_scripts/util/` directory)

5) Create a copy of the `./config.yml` file and adjust the parameters to suit your analysis 

5) Execute scripts of interest and collect results

## Contact

Kyle Wellband - kyle.wellband [at] gmail.com

## License

CC share-alike

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">methylUtil</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Kyle Wellband</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
