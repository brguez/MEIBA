#!/bin/bash

mkdir -p logs
mkdir -p tmp
tail -n +2 $1 | while read CHROM POS CLASS PROJECT DONOR; do
    echo -e "$CHROM\t$POS\t$CLASS\t$PROJECT\t$DONOR" > tmp/${CLASS}_${CHROM}_${POS}.input
    bsub -o $PWD/logs/log.%J -e $PWD/logs/err.%J -q basement -n 1 -M15000 -R"span[hosts=1] select[mem>=15000] rusage[mem=15000]" "/nfs/dog_n_devil/adrian/software/scripts/MEI_genotyping/extract_reads_1000bp_mates.sh $PWD/tmp/${CLASS}_${CHROM}_${POS}.input /lustre/scratch109/sanger/jt14/Berni/Projects/Pancancer/MEI_germline/Analyses/MEI_genotyping/donorId_bamPath.txt"
done
