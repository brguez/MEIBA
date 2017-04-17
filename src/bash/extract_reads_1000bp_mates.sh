#!/bin/bash

# Create one miniBAM per insertion, taking all the reads in +/-1000 bp around the 
# insertion position, AND their mates

# $1: tab-delimited file with no header and columns: chrom, pos, class, projectCode, donorId
# $2: tab-delimited file with no header and columns: donorId, bamFile


# For each insertion
while read CHROM POS CLASS PROJECT DONOR; do
    
    echo -e "\nProcessing: $CHROM $POS $CLASS $PROJECT $DONOR"
    
    # Get path to BAM file
    BAM=`grep -P "^$DONOR\t" $2 | cut -f2`
    echo "BAM file: $BAM"
    
    # Extract read IDs in interval of +/-1000 bp
    START=$(( POS - 1000 ))
    END=$(( POS + 1000 ))
    echo "Running samtools view in interval ${CHROM}:${START}-${END}"
    samtools view $BAM ${CHROM}:${START}-${END} | cut -f1 > ${CHROM}_$POS.reads
    NUMREADS=`cat ${CHROM}_$POS.reads | wc -l`
    NUMREADS=$(( NUMREADS * 2 ))
    
    # Extract all the reads and their mates
    echo "Extracting all found reads and their mates ($NUMREADS)"
    samtools view -H $BAM > ${CHROM}_$POS.sam
    samtools view $BAM | grep -m $NUMREADS -Ff ${CHROM}_$POS.reads >> ${CHROM}_$POS.sam
    echo "Compressing to BAM"
    samtools view -Sb ${CHROM}_$POS.sam > FDR_${CLASS}_${CHROM}_$POS.bam
    rm ${CHROM}_$POS.reads ${CHROM}_$POS.sam
    echo "Done"
    
done < $1

echo -e "\nALL DONE"
