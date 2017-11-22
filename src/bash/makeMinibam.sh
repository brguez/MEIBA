#!/bin/bash

##### makeMinibam
## makeMinibam is aimed to obtain a minibam containing the reads corresponding to a list of regions for a given sample. It takes as input three arguments:

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

### 1) Generate sam file with header
#####################################
echo
echo "** 1) Generate sam file with header **"

$SAMTOOLS view -H $bam > $outDir/${fileName}_minibam.sam


### 2) Add the alignments in the target regions to the sam generated in 1)
###########################################################################
echo
echo "** 2) Add the alignments in the target regions to the sam generated in 1) **"

awk '$0 !~ /^#/{chrom=$1; beg=$2; end=$3; region=chrom":"beg"-"end; print region;}' $bed | while read region
do
    #echo "Extracting read-pair ids from $region region"
    $SAMTOOLS view $bam $region >> $outDir/${fileName}_minibam.sam
done


### 3) Convert sam into bam
############################
echo
echo "** 3) Convert sam into bam **"

$SAMTOOLS view -Sb $outDir/${fileName}_minibam.sam > $outDir/${fileName}_minibam.bam


### 4) Sort minibam 
####################
echo
echo "** 4) Sort minibam **"

$SAMTOOLS sort $outDir/${fileName}_minibam.bam > $outDir/${fileName}_minibam.sorted.bam


### 5) Index minibam file
#########################
echo
echo "** 5) Index minibam file **"

$SAMTOOLS index $outDir/${fileName}_minibam.sorted.bam >$outDir/${fileName}_minibam.sorted.bam.bai

### Make cleanup
rm $outDir/${fileName}_minibam.sam
rm $outDir/${fileName}_minibam.bam


echo
echo "** END **"
