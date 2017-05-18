#!/bin/bash

# Description
##############
# Extract reads overlapping insertion breakpoints from a BAM file.
# Read IDs are collected from a TraFiC insertions file and corresponding
# reads are output as FASTA (in order of appeareance in BAM).

# Usage
#########
# bash extract_reads_bam2fa insertions.tsv sample.bam /path/to/outdir

# Input
########
# 1) TraFiC output file (TSV) containing insertion infos
# 2) BAM read alignments
# 3) output directory

# Output
#########
# 1) FASTA file containing reads referenced in TraFiC insertions file.

scriptname=$(basename $0)

function writeUsage {
    echo "" &>2
    echo "*** reads extractor ***" >&2
    echo "usage: $scriptname trafic_insertions.tsv reads.bam out_dir" >&2
    echo "" &>2
}
function writeError {
    printf "[ERROR] ($scriptname) $1\n" >&2
}
function writeInfo {
    printf "[INFO] ($scriptname) $1\n"
}


# Check input params
if [[ $# < 3 ]]; then
    writeUsage
    exit -1
elif [[ ! -f $1 ]]; then
    writeError "file does not exist: $1"
    exit -1
elif [[ ! -f $2 ]]; then
    writeError "file does not exist: $2"
    exit -1
elif [[ ! -d $3 ]]; then
    writeError "directory does not exist: $3"
    exit -1
fi

# get input/output options
inInsertions=$1
inBam=$2
outDir=$3
outSam=${outDir}/targetReads.sam
outBam=${outDir}/targetReads.bam
outTxt=${outDir}/readIDs.txt
outFasta=${outDir}/targetReads.fasta

writeInfo "reading read pair IDs from: $inInsertions"

# get read IDs from TSV file (store them in memory) and extract reads from the BAM
# conditions test the following:
#   (NR==FNR):   Reading first file (insertions TSV)
#   (!/^#/):     line does not start with '#' (e.g. headers, comments)
#   (/^>/):      Reading FASTA header line
#   (a[1] in r): part of FASTA header is in array r (stores read pair IDs)

awkCmd='
(NR==FNR) && (!/^#/) {
    split($6","$12, a, ",");
    for(x in a) r[a[x]];
    next
}
(/^>/) {
    split(substr($1,2), a, "/");
    if (a[1] in r) {
        print > outFasta;
        getline;
        print > outFasta
        c++
    }
}
END {
    for (x in r) e++;
    print "[INFO] ('$scriptname') reads expected:  "e*2
    print "[INFO] ('$scriptname') reads extracted: "c
}
'
awk -v outFasta=$outFasta "$awkCmd" \
 $inInsertions \
 <(samtools fasta $inBam)

writeInfo "wrote output to: $outFasta"
