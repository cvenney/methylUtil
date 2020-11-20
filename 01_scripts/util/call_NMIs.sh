#!/bin/bash
# call_NMIs.sh

INDIR="04_filtered_bedGraphs"

files=($(ls -1 $INDIR | grep "bedGraph"))

if [[ ! ${#files[@]} -gt 1 ]]
then
    echo "Error: no bedGraph files found!"
    exit
fi

echo "Found ${#files[@]} bedGraph files, summarizing coverage"

mkdir -p $INDIR/tmp

for f in ${files[@]}
do

    if [[ $(echo $f | grep "\.gz$") != "" ]]
    then
        gunzip -c $INDIR/$f | awk -vOFS="\t" '{cov = $5 + $6; print $1, $2, $3, cov}' > $INDIR/tmp/${f%.*}_cov.bg
        gunzip -c $INDIR/$f | awk -vOFS="\t" '{print $1, $2, $3, $5}' > $INDIR/tmp/${f%.*}_M.bg
    else
        awk -vOFS="\t" '{cov = $5 + $6; print $1, $2, $3, cov}' $INDIR/$f > $INDIR/tmp/${f}_cov.bg
        awk -vOFS="\t" '{print $1, $2, $3, $5}' $INDIR/$f > $INDIR/tmp/${f}_M.bg
    fi
       
done

for i in cov M
do
    
    echo "Combining $i..."
    bedtools unionbedg -i $INDIR/tmp/*${i}.bg | \
    awk -vOFS="\t" '{tot = 0; for(c=4;c<=NF;c++) tot += $c; print $1, $2, $3, tot}' > $INDIR/tmp/combined_${i}.bedGraph
    echo "\b done!"
    
done

bedtools unionbedg -i $INDIR/tmp/combined_M.bedGraph $INDIR/tmp/combined_cov.bedGraph | \
    awk -vOFS="\t"  '{N = $5 - $4; P = int(100 * $4 / $5); print $1, $2, $3, P, $4, N}' | gzip > $INDIR/combined.bedGraph.gz

rm -r $INDIR/tmp


