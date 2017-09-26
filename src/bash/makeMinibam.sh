#!/bin/bash

##### makeMinibam
## makeMinibam is aimed to obtain a minibam containing the reads corresponding to a list of regions and their mates for a given sample. It takes as input three arguments:

# 1) Bed file containing the list of regions of interest
# 2) Bam file
# 2) Output file name
# 3) Output directory

## Get input parameters
bed=$1
bam=$2
fileName=$3
outDir=$4

echo
echo "******* Configuration *******"
echo "bed: " $bed
echo "bam: " $bam
echo "fileName: " $fileName
echo "outDir: " $outDir
echo "*****************************"
echo


#### Dependencies
## Samtools
SAMTOOLS="/Users/brodriguez/Research/Apps/samtools/samtools-1.3.1/samtools" # Laptop
#SAMTOOLS="/software/vertres/bin-external/samtools-1.4" # Sanger

## Picard 
PICARD="java -jar /Users/brodriguez/Research/Apps/Picard/2.12.1/picard.jar FilterSamReads" # Laptop
#PICARD="java -Xms10G -Xmx10G -jar /software/CGP/external-apps/picard-tools-1.80/lib/FilterSamReads.jar" # Sanger


### 1) Generate a text file containing the read-pair ids 
#########################################################
# from all the regions of interest
###################################

echo
echo "** 1) Generate a text file containing the read-pair ids from all the regions of interest **"

awk '$0 !~ /^#/{chrom=$1; beg=$2; end=$3; region=chrom":"beg"-"end; print region;}' $bed | while read region
do
    #echo "Extracting read-pair ids from $region region"
    $SAMTOOLS view $bam $region | cut -f1 >> $outDir/${fileName}_readIds.tmp
done


### 2) Generate a minibam file only containing the read-pairs from 1)
######################################################################

echo
echo "** 2) Generate a minibam file only containing the read-pairs from 1) **"

$PICARD I=$bam O=$outDir/${fileName}_minibam.bam READ_LIST_FILE=$outDir/${fileName}_readIds.tmp FILTER=includeReadList WRITE_READS_FILES=false VALIDATION_STRINGENCY=SILENT

### Make cleanup
rm $outDir/${fileName}_readIds.tmp


