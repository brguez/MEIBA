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


### 0) If 0 regions generate an empty bam file containing header and exit
#############################################################################

nbRegions=`grep -v '^#' $bed | wc -l`

if [ $nbRegions == 0 ]
then
    echo "** 0 input regions. Generate an empty bam file containing header and exit ** " 
    $SAMTOOLS view -Hb -o $outDir/${fileName}_minibam.bam $bam
    echo "FINISHED!!"
    exit 0
fi



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


### 3) Sort minibam 
####################
echo
echo "** 3) Sort minibam **"

$SAMTOOLS sort $outDir/${fileName}_minibam.bam > $outDir/${fileName}_minibam.sorted.bam


### 4) Index minibam file
#########################
echo
echo "** 4) Index minibam file **"

$SAMTOOLS index $outDir/${fileName}_minibam.sorted.bam >$outDir/${fileName}_minibam.sorted.bam.bai


### Make cleanup
rm $outDir/${fileName}_readIds.tmp
rm $outDir/${fileName}_minibam.bam

echo "FINISHED!!"

