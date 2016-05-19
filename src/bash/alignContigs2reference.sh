#!/bin/bash


<<authors
******************************************************************************
	
	alignContigs2reference.sh
	
	Copyright (c) 2016 Bernardo Rodríguez-Martín
	
	Mobile Genomes & Disease Lab.
	Universidad de Vigo (Spain)

	Licenced under the GNU General Public License 3.0 license.
******************************************************************************
authors

# Description
##############
# Aligns a set of contigs on the TE inserted region and the TE consensus sequence.
# Contigs produced through the assembly of TraFiC + or - cluster supporting read-pairs

# usage
#######
# Usage:    bash alignContigs2reference.sh contigs.fa insertionId genome.fa TE.fa [windowSize outDir]  
# Example:  bash alignContigs2reference.sh L1:8_37122847_37122890:-.contigs.fa L1:8_37122847_37122890:- hs37d5.fa L1_consensus.fa 1000 blat_test

# Input
########
# 1) Assembled contigs fasta file.
# 2) Transposable element insertion identifier. Format: '${family}:${chr}_${beg}_${end}:${orientation}.fa'
## where:
#	- family: TE family (L1, ALU, SVA...)
#	- chr: insertion chromosome
# 	- beg: insertion beginning
#	- end: insertion end
#	- orientation: cluster (+ or -)
# 3) Reference genome fasta file. Mandatory.
# 4) Consensus transposable element sequence (L1, Alu or SVA depending on the TE's family).

# Output
###########
# .psl file with contig blat aligments.
# 286	0	0	0	0	0	0	0	+	NODE_2_length_614_cov_6.804560	634	0	286	L1	na	37127581	37127867	1	286,	0,	37127581,
# 330	0	0	0	0	0	0	0	+	NODE_2_length_614_cov_6.804560	634	304	634	8	na	37122889	37123219	1	330,	304,	37122889,


### will exit if there is an error or in a pipe
set -e -o pipefail

# In case the user does not provide any input file
###################################################
if [[ ! -e "$1" ]] || [[ "$2" == "" ]] || [[ ! -e "$3" ]] || [[ ! -e "$4" ]] 
then
    echo "" >&2
    echo "*** alignContigs2reference ***" >&2
    echo "" >&2
    echo "Usage:    bash alignContigs2reference.sh contigs.fa insertionId genome.fa TE.fa [windowSize outDir]  " >&2
    echo "" >&2
    echo "Example:  bash alignContigs2reference.sh L1:8_37122847_37122890:-.contigs.fa L1:8_37122847_37122890:- hs37d5.fa L1_consensus.fa 1000 blat_test" >&2
    echo "" >&2
    echo "Aligns a set of contigs on the TE inserted region and the TE consensus sequence." >&2
    echo "Contigs produced through the assembly of TraFiC + or - cluster supporting read-pairs" >&2
    echo "" >&2
    echo "Input:" >&2
    echo "1) Assembled contigs fasta file. Mandatory" >&2
    echo "2) Transposable element insertion identifier. Format: '${family}:${chr}_${beg}_${end}:${orientation}.fa'. Mandatory" >&2
    echo "3) Reference genome fasta file. Mandatory" >&2
    echo "4) Consensus transposable element sequence (L1, Alu or SVA depending on the TE's family). Mandatory"  >&2
    echo "5) Window size around insertion breakpoints to define blat target region. Default=1000"  >&2
    echo "6) Output directory. Default: current working directory"  >&2
    echo "" >&2
    echo "Output:" >&2 
    echo "1) .psl file with contig blat aligments." >&2
    echo "" >&2
    exit -1
fi


# In case the user does not provide windowSize or outDir
#########################################################
# provide default values
########################
contigs=$1
insertionId=$2
genome=$3
TEseq=$4

if [ ! -n "$5" ]
then
	windowSize=1000
	outDir=.
else
	windowSize=$5
	
	if [ ! -n "$6" ]
	then
		outDir=.
	else
		outDir=$6
	fi
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

# Programs and scripts
########################

## Awk
ADDOFFSET=$awkDir/addOffset2template.awk

## DISPLAY PROGRAM CONFIGURATION  
##################################
printf "\n"
header="CONFIGURATION FOR $insertionId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n" "***** MANDATORY ARGUMENTS *****"
printf "  %-34s %s\n" "contigs:" "$contigs"
printf "  %-34s %s\n" "insertion identifier:" "$insertionId"
printf "  %-34s %s\n" "ref. genome:" "$genome"
printf "  %-34s %s\n\n" "TE ref. sequence :" "$TEseq"
printf "  %-34s %s\n" "***** OPTIONAL ARGUMENTS *****"
printf "  %-34s %s\n" "window size:" "$windowSize"
printf "  %-34s %s\n\n" "outDir:" "$outDir"
printf "\n\n"


##########
## START #
##########

header="Executing alignContigs2reference for $insertionId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"


###########################################################
# 1. MAKE FASTA WITH TARGET DNA REGION FOR BLAT ALIGNMENT # 
###########################################################
## Output:
# - $outDir/insertion_region.fa
targetRegionPath=$outDir/insertion_region.fa 

echo "1. Make fasta with target dna region for blat alignment" >&1

read family chr beg end cluster <<<$(echo $insertionId | awk '{split($1, info, ":"); family=info[1]; cluster=info[3]; split(info[2], coord, "_"); chr=coord[1]; beg=coord[2]; end=coord[3]; print family, chr, beg, end, cluster;}')

targetBeg=`expr $beg - $windowSize` 
targetEnd=`expr $end + $windowSize`
targetInterval=$chr":"$targetBeg"-"$targetEnd
offset=$targetBeg

echo "samtools faidx $genome $targetInterval > $targetRegionPath" >&1
samtools faidx $genome $targetInterval > $targetRegionPath


###########################################################################
# 2. CONCATENATE TARGET REGION FASTA WITH THE L1 REFERENCE SEQUENCE FASTA #
###########################################################################
## Output:
# - $outDir/target_sequences.fa
targetSeqPath=$outDir/target_sequences.fa

echo "2. Concatenate target region fasta with TE reference sequence fasta" >&1

echo "cat $outDir/insertion_region.fa $TEseq > $targetSeqPath" >&1
cat $outDir/insertion_region.fa $TEseq > $targetSeqPath


##################################################
# 3. BLAT CONTIGS INTO THE FASTA GENERATED IN 2. #
##################################################
## Output:
# - $outDir/$insertionId".tmp.psl"
tmpBlatPath=$outDir/$insertionId".tmp.psl" 

echo "3. Blat contigs into the fasta generated in 2." >&1

echo "blat -t=dna -q=dna -stepSize=5 -minScore=20 -out=psl -noHead $outDir/target_sequences.fa $contigs $tmpBlatPath" >&1
blat -t=dna -q=dna -stepSize=5 -minScore=20 -out=psl -noHead $outDir/target_sequences.fa $contigs $tmpBlatPath


####################################################################
# 4. CONVERT PSL TEMPLATE COORDINATES TO GENOMIC ONES (ADD OFFSET) #
####################################################################
## Output:
# - $outDir/$insertionId".psl"
blatPath=$outDir/$insertionId".psl"

echo "4. Convert psl template coordenates to genomic coordenates (add offset) " >&1

echo "awk -v OFS='\t' -v offset=$offset -f $ADDOFFSET $outDir/$insertionId".tmp.psl" > $blatPath" >&1
awk -v OFS='\t' -v offset=$offset -f $ADDOFFSET $outDir/$insertionId".tmp.psl" > $blatPath


######################
# 5. CLEANUP AND END #
######################
echo "5. Cleanup and end" >&1
echo >&1
rm $targetRegionPath $targetSeqPath $tmpBlatPath 



