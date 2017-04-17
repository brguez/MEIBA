#!/bin/bash

# Create one miniBAM per insertion, taking all the reads in +/-1000 bp around the 
# insertion position, AND their mates

# $1: tab-delimited file with no header and columns: chrom, pos, class, projectCode, donorId
# $2: BAM file path
# $3: output directory

input=$1
BAM=$2
outDir=$3 

# For each insertion
while read CHROM POS CLASS PROJECT DONOR; do
    
    echo -e "\nProcessing: $CHROM $POS $CLASS $PROJECT $DONOR"
    
    # Get path to BAM file
    echo "BAM file: $BAM"
    
    # Extract read IDs in interval of +/-1000 bp
    START=$(( POS - 1000 ))
    END=$(( POS + 1000 ))
    echo "Running samtools view in interval ${CHROM}:${START}-${END}"
    echo "samtools view $BAM ${CHROM}:${START}-${END} | cut -f1 > $outDir/${CHROM}_${POS}_${CLASS}.reads"
    samtools view $BAM ${CHROM}:${START}-${END} | cut -f1 > $outDir/${CHROM}_${POS}_${CLASS}.reads
    NUMREADS=`cat $outDir/${CHROM}_${POS}_${CLASS}.reads | wc -l`
    NUMREADS=$(( NUMREADS * 2 ))
    
    # Extract all the reads and their mates
    echo "Extracting all found reads and their mates ($NUMREADS)"
    echo "samtools view -H $BAM > $outDir/${CHROM}_${POS}_${CLASS}.sam"    
    samtools view -H $BAM > $outDir/${CHROM}_${POS}_${CLASS}.sam
    echo "samtools view $BAM | grep -m $NUMREADS -Ff $outDir/${CHROM}_${POS}_${CLASS}.reads >> $outDir/${CHROM}_${POS}_${CLASS}.sam" 
    samtools view $BAM | grep -m $NUMREADS -Ff $outDir/${CHROM}_${POS}_${CLASS}.reads >> $outDir/${CHROM}_${POS}_${CLASS}.sam
    echo "Compressing to BAM"
    echo "samtools view -Sb $outDir/${CHROM}_${POS}_${CLASS}.sam > $outDir/${CHROM}_${POS}_${CLASS}.bam"    
    samtools view -Sb $outDir/${CHROM}_${POS}_${CLASS}.sam > $outDir/${CHROM}_${POS}_${CLASS}.bam
    rm $outDir/${CHROM}_${POS}_${CLASS}.reads $outDir/${CHROM}_${POS}_${CLASS}.sam
    echo "Done"
    
done < $input

echo -e "\nALL DONE"
