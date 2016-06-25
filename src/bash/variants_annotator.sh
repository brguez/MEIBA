#!/bin/bash


<<authors
******************************************************************************
	
	variants_annotator.sh
	
	Copyright (c) 2016 Bernardo Rodríguez-Martín
	
	Mobile Genomes & Disease Lab.
	Universidad de Vigo (Spain)

	Licenced under the GNU General Public License 3.0 license.
******************************************************************************
authors

# Description
##############
# 

# usage
#######
# Usage:    bash variants_annotator.sh sample.vcf repeatsDb.bed sampleId [outDir]  
# Example:  bash variants_annotator.sh e52ffa79-557a-4024-81f3-f3826c227ec5.vcf repeats_repeatMasker_hg19.bed e52ffa79-557a-4024-81f3-f3826c227ec5 variant_annotation_dir

# Input
########
# 1) VCF
# 2) Sample identifier
# 3) Output directory

# Output
###########
# 1) VCF with annotated MEI (overlapping region, gene, repeat, satellite region...)

### will exit if there is an error or in a pipe
set -e -o pipefail

# In case the user does not provide any input file
###################################################
if [[ ! -e "$1" ]] || [[ ! -e "$2" ]] || [[ "$3" == "" ]] 
then
    echo "" >&2
    echo "*** variants_annotator ***" >&2
    echo "" >&2
    echo "Usage:    bash variants_annotator.sh sample.vcf repeatsDb.bed sampleId [outDir] " >&2
    echo "" >&2
    echo "Example:  bash variants_annotator.sh e52ffa79-557a-4024-81f3-f3826c227ec5.vcf repeats_repeatMasker_hg19.bed e52ffa79-557a-4024-81f3-f3826c227ec5 variant_annotation_dir" >&2
    echo "" >&2
    echo "...description..." >&2
    echo "" >&2
    echo "Input:" >&2
    echo "1) VCF" >&2
    echo "2) RepeatMasker repeats database in bed format" >&2	
    echo "3) Sample identifier" >&2
    echo "4) Output directory" >&2
    echo "" >&2
    echo "Output:" >&2 
    echo "1) VCF with annotated variants (overlapping region, gene, repeat, satellite region...)" >&2
    echo "" >&2
    exit -1
fi


# In case the user does not provide windowSize or outDir
#########################################################
# provide default values
########################
inputVCF=$1
repeatsDb=$2
donorId=$3


if [ ! -n "$4" ]
then
	outDir=.
else
	outDir=$4
fi	


# Directories 
#############
## Set root directory
path="`dirname \"$0\"`"              # relative path
rootDir="`( cd \"$path\" && pwd )`"  # absolute path

if [ -z "$rootDir" ] ; 
then
  # error; for some reason, the path is not accessible
  # to the script
  log "Path not accessible to the script\n" "ERROR" 
  exit 1  # fail
fi

## Set awk scripts directory
awkDir=$rootDir/../awk

## Set python scripts directory
pyDir=$rootDir/../python

## Annovar directory
annovarDir=$rootDir/../../apps/annovar

# Programs, scripts and databases
##################################

## Awk
VCF2ANNOVAR=$awkDir/vcf2annovar.MEI.awk
VCF2BED=$awkDir/vcf2bed.MEI.awk

## Python
ADD_GENE2VCF=$pyDir/addGnAnnot2VCF.py
ADD_REPEAT2VCF=$pyDir/addRepeatAnnot2VCF.py

## Annovar:
ANNOVAR=$annovarDir/annotate_variation.pl
ANNOVAR_DB=$annovarDir/humandb/

## DISPLAY PROGRAM CONFIGURATION  
##################################
printf "\n"
header="CONFIGURATION FOR $donorId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n" "***** MANDATORY ARGUMENTS *****"
printf "  %-34s %s\n" "inputVCF:" "$inputVCF"
printf "  %-34s %s\n" "repeatsDb:" "$repeatsDb"
printf "  %-34s %s\n\n" "donorId:" "$donorId"
printf "  %-34s %s\n" "***** OPTIONAL ARGUMENTS *****"
printf "  %-34s %s\n\n" "outDir:" "$outDir"
printf "\n\n"


##########
## START #
##########

header="Executing variants_annotator for $donorId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"


######################
# 1. GENE ANNOTATION # 
######################

echo "1. Gene annotation" >&1

# 1.1 Convert variants into annovar input format
#################################################

## Output:
# - $outDir/annovar.input.txt
annovarInputPath=$outDir/annovar.input.txt

echo "1.1 Convert variants into annovar input format" >&1

echo "awk -f $VCF2ANNOVAR $inputVCF > $annovarInputPath" >&1
awk -f $VCF2ANNOVAR $inputVCF > $annovarInputPath


# 1.2 Annotate variants with annovar
######################################
## Output:
# - $outDir/annovar.variant_function

annovarOut="$outDir/annovar.variant_function"

echo "1.2 Annotate variants with annovar" >&1

echo "perl $ANNOVAR -build hg19 -out $outDir/annovar -dbtype wgEncodeGencodeBasicV19 $annovarInputPath $ANNOVAR_DB 1> $outDir/annovar.out 2> $outDir/annovar.err"  >&1
perl $ANNOVAR -build hg19 -out $outDir/annovar -dbtype wgEncodeGencodeBasicV19 $annovarInputPath $ANNOVAR_DB 1> $outDir/annovar.out 2> $outDir/annovar.err

## 1.3 Add annotation information to the VCF
#############################################
## Output:
# -  $outDir/$donorId".gnAnnot.vcf"

echo "1.3 Add annotation information to the VCF" >&1

echo "python $ADD_GENE2VCF $inputVCF $annovarOut $donorId $outDir" >&1
python $ADD_GENE2VCF $inputVCF $annovarOut $donorId --outDir $outDir 1> $outDir/addGnAnnot.out 2> $outDir/addGnAnnot.err


#########################
# 2. REPEATS ANNOTATION # 
#########################

echo "2. Repeat annotation" >&1

# 2.1 Convert MEI calls into bed format
########################################

## Output:
# - $outDir/$donorId.bed
insertionsBed=$outDir/$donorId.bed

echo "2.1 Convert MEI calls into bed format" >&1

echo "awk -f $VCF2BED $outDir/$donorId.gnAnnot.vcf  > $insertionsBed" >&1
awk -f $VCF2BED $outDir/$donorId.gnAnnot.vcf  > $insertionsBed


# 2.2  Interserct MEI with repeats database
#####################################################
## Output:
# - $outDir/insertions_repeatAnnot.txt

repeatAnnot=$outDir/insertions_repeatAnnot.txt

echo "2.2  Interserct MEI with repeats database" >&1

echo "bedtools intersect -wao -a $insertionsBed -b $repeatsDb | awk -v OFS='\t' '{print $1, $2, $3, $4, $8, $9}' > $repeatAnnot"  >&1
bedtools intersect -wao -a $insertionsBed -b $repeatsDb | awk -v OFS='\t' '{print $1, $2, $3, $4, $8, $9}' > $repeatAnnot


## 2.3 Add repeats annotation information to the VCF
#####################################################
## Output:
# -  $outDir/$donorId.annotated.vcf

echo "2.3 Add repeats annotation information to the VCF" >&1

echo "python $ADD_REPEAT2VCF $outDir/$donorId.gnAnnot.vcf $repeatAnnot $donorId $outDir 1> $outDir/addRepeatAnnot.out 2> $outDir/addRepeatAnnot.err" >&1
python $ADD_REPEAT2VCF $outDir/$donorId.gnAnnot.vcf $repeatAnnot $donorId --outDir $outDir 1> $outDir/addRepeatAnnot.out 2> $outDir/addRepeatAnnot.err

######################
# 3. CLEANUP AND END #
######################
echo "3. Cleanup and end" >&1
echo >&1
# echo "rm " >&1




