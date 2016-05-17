#!/bin/bash


<<authors
******************************************************************************
	
	TE-bkpAssembler.sh
	
	Copyright (c) 2016 Bernardo Rodríguez-Martín
	
	Mobile Genomes and Disease group.
	Universidad de Vigo (Spain)

	Licenced under the GNU General Public License 3.0 license.
******************************************************************************
authors


# Function 1. Print basic usage information
############################################
function usageDoc
{
cat <<help
	

**** TE-bkpAssembler.sh	version $version ****
Execute  for one dataset (sample).
	
*** USAGE

	$0 -i <TraFic_insertions> -f <FASTA> -s <sample_identifier> [OPTIONS]

*** MANDATORY 
		
	-i	<TXT>			TraFiC TE somatic insertion calls for a given sample.		

	-f	<FASTA>			Fasta containing TE insertions supporting reads. 

	-g 	<GENOME>		Reference genome (RG). Please make sure you provide the same RG version you used to run TraFiC. 
					Also, make sure the same chromosome naming conventions are used.
	
	-s	<STRING>		Sample id. Output file will be named accordingly.	
		
*** [OPTIONS] can be:
* General:
 
	-k 	<INTEGER>        	K-mer length of the words being hashed for the assembly. Default: N=21.

	-o	<PATH>			Output directory. Default current working directory. 
	
	-h				Display usage information.
		

help
}

# Function 2. Parse user's input
################################
function getoptions {

while getopts ":i:f:g:s:k:o:h" opt "$@"; 
do
   case $opt in   	
      
      ## MANDATORY ARGUMENTS
      i)
	  if [ -n "$OPTARG" ];
	  then
              input=$OPTARG
	  fi
	  ;;
      
      f)
	  if [ -n "$OPTARG" ];
	  then
              fasta=$OPTARG
	  fi
	  ;;

      g)
	  if [ -n "$OPTARG" ];
	  then
              genome=$OPTARG
	  fi
	  ;;      

      s)
	  if [ -n "$OPTARG" ];
	  then
              sampleId=$OPTARG
	  fi
	  ;;
      
      ## OPTIONS
      k)
	  if [ -n "$OPTARG" ];
	  then
              kmerLen=$OPTARG
	  fi
	  ;;

      o)
	  if [ -n "$OPTARG" ];
	  then
       	      outDir=$OPTARG
	  fi
	  ;;

      h)
	  usageDoc;
	  exit 1
	  ;;
      
      :)
          echo "Option -$OPTARG requires an argument." >&2
          exit 1
  esac
done
}

# Function 3. Print log information (Steps and errors)
#######################################################
function log {
    string=$1
    label=$2
    if [[ ! $ECHO ]];then
        if [[ "$label" != "" ]];then
            printf "[$label] $string"
        else
            printf "$string"
        fi
    fi
}

# Function 4. Print a section header for the string variable
##############################################################
function printHeader {
    string=$1
    echo "`date` ***** $string *****"
}

# Function 5. Print a subsection header for the string variable
################################################################
function printSubHeader {
    string=$1
    echo "`date` * $string *"
}

# Function 6. Execute and print to stdout commands 
###################################################
function run {
    command=($1)
    if [[ $2 ]];then
         ${2}${command[@]}
    else
        echo -e "\t"${command[@]}""
        eval ${command[@]}
    fi
}


# SETTING UP THE ENVIRONMENT
############################

# TE-bkpAssembler version 
version=v1

# Enable extended pattern matching 
shopt -s extglob

# 1. Root directory
##############################
# to set the path to the bin and python directories. 

path="`dirname \"$0\"`"              # relative path
rootDir="`( cd \"$path\" && pwd )`"  # absolute path

if [ -z "$rootDir" ] ; 
then
  # error; for some reason, the path is not accessible
  # to the script
  log "Path not accessible to the script\n" "ERROR" 
  exit 1  # fail
fi


# 2. Parse input arguments with getopts  
########################################

# A) Display help and exit if no input argument is provided
if [ $# -eq 0 ];
then
    usageDoc
    exit 0
else
    getoptions $@ # call Function 2 and passing two parameters (name of the script and command used to call it)	
fi

# 3. Check input variables 
##########################

## Mandatory arguments
## ~~~~~~~~~~~~~~~~~~~

if [[ ! -e $input ]]; then log "The TraFiC TE insertion calls file does not exist. Mandatory argument -i\n" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -e $fasta ]]; then log "The TE insertion supporting reads fasta file does not exist. Mandatory argument -f" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -e $genome ]]; then log "The reference genome fasta file does not not exist. Mandatory argument -g" "ERROR" >&2; usageDoc; exit -1; fi
if [[ $sampleId == "" ]]; then log "The sample id is not provided. Mandatory argument -s\n" "ERROR" >&2; usageDoc; exit -1; fi


## Optional arguments
## ~~~~~~~~~~~~~~~~~~

# K-mer length for assembly
if [[ "$kmerLen" == "" ]]; 
then 
	kmerLen='21'; 
fi

# Output directory
if [[ "$outDir" == "" ]]; 
then 
	outDir=${SGE_O_WORKDIR-$PWD};
else
	if [[ ! -e "$outDir" ]]; 
	then
		log "Your output directory does not exist. Option -o\n" "ERROR" >&2;
		usageDoc; 
		exit -1; 
	fi	
fi


# 4. Directories
################
## binaries and scripts
srcDir=$rootDir/src
binDir=$rootDir/bin

## Output files directories
fastaDir=$outDir/Fasta
contigsDir=$outDir/Contigs
assemblyLogsDir=$contigsDir/logs
blatDir=$outDir/Blat

# 5. Programs/Scripts
######################
CLUSTERS2FASTA=$srcDir/clusters2fasta.py
VELVETH=$binDir/velveth
VELVETG=$binDir/velvetg

## DISPLAY PROGRAM CONFIGURATION  
##################################
printf "\n"
header=" TE-bkpAssembler CONFIGURATION FOR $sampleId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n\n" "TE-BkpAssembler $version"
printf "  %-34s %s\n" "***** MANDATORY ARGUMENTS *****"
printf "  %-34s %s\n" "input:" "$input"
printf "  %-34s %s\n" "fasta:" "$fasta"
printf "  %-34s %s\n" "genome:" "$genome"
printf "  %-34s %s\n\n" "sampleId:" "$sampleId"
printf "  %-34s %s\n" "***** OPTIONAL ARGUMENTS *****"
printf "  %-34s %s\n" "K-mer length:" "$kmerLen"
printf "  %-34s %s\n\n" "outDir:" "$outDir"
	 
	
##########
## START #
##########
header="Executing TE-bkpAssembler $version for $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
start=$(date +%s)

logFile=$outDir/TE-bkpAssembler_${sampleId}.log

#######################    	
# 0) PRELIMINARY STEPS #
#######################

## 0.1) Make directories
#######################
# Fastas for assembly
if [[ ! -d $fastaDir ]]; then mkdir $fastaDir; fi

# Assembled contigs of TE insertion breakpoints and assembly logs
if [[ ! -d $contigsDir ]]; then mkdir $contigsDir; fi
if [[ ! -d $assemblyLogsDir ]]; then mkdir $assemblyLogsDir; fi

# Contigs alignment with Blat:
if [[ ! -d $blatDir ]]; then mkdir $blatDir; fi


# 1) Produces for each TE insertion two fasta (one for read pairs supporting + cluster and another one for read pairs supporting - cluster)
############################################################################################################################################
# These fasta files will be used to assemble the 5' and 3' TE insertion breakpoint sequences 
#############################################################################################
# Each fasta will be named according to this convention: 
########################################################
# - ${fastaDir}/${family}:${chr}_${beg}_${end}:${orientation}.fa
# where:
#	- family: TE family (L1, ALU, SVA...)
#	- chr: insertion chromosome
# 	- beg: insertion beginning
#	- end: insertion end
#	- orientation: cluster (+ or -)
 	
step="CLUSTERS2FASTA"
startTime=$(date +%s)
printHeader "Producing per TE insertion two fasta for insertion bkp assembly"  
run "python $CLUSTERS2FASTA $input $fasta --outDir $fastaDir > $logFile" "$ECHO"	
endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


# 2) Assemble the 5' and 3' TE insertion breakpoints with velvet. An independent assembly is performed for each insertion breakpoint 
######################################################################################################################################
# It will produce a fasta containing the assembled contigs for each cluster insertion
######################################################################################
# Fasta naming convention:
########################### 
# - ${contigsDir}/${family}:${chr}_${beg}_${end}:${orientation}.contigs.fa
# where:
#	- family: TE family (L1, ALU, SVA...)
#	- chr: insertion chromosome
# 	- beg: insertion beginning
#	- end: insertion end
#	- orientation: cluster (+ or -)

step="BKP-ASSEMBLY"
startTime=$(date +%s)
printHeader "Assembling the 5' and 3' TE insertion breakpoints with velvet"  

ls $fastaDir | grep '.*fa' | while read bkpFasta; 
do 	
	bkpId=${bkpFasta%.fa}	
	fastaPath=${fastaDir}/${bkpFasta}	
	contigPath=${contigsDir}/${bkpId}".contigs.fa"
		
	log "** ${bkpId} breakpoint **\n" $step
	log "1. Preparing files for assembly\n" $step
	run "$VELVETH $contigsDir $kmerLen -fasta -short $fastaPath >> $logFile " "$ECHO"
	log "2. Breakpoint assembly with velvet\n" $step
	run "$VELVETG $contigsDir -exp_cov auto -cov_cutoff auto >> $logFile" "$ECHO"
	log "3. Rename output files\n" $step
	run "mv ${contigsDir}/Log $assemblyLogsDir/${bkpId}.log" "$ECHO"
	run "mv ${contigsDir}/contigs.fa ${contigsDir}/${bkpId}.contigs.fa" "$ECHO"
	log "4. Cleaning\n" $step
	run "rm ${contigsDir}/Sequences ${contigsDir}/Roadmaps ${contigsDir}/PreGraph ${contigsDir}/stats.txt ${contigsDir}/LastGraph ${contigsDir}/Graph2" "$ECHO"
done

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


# 3) Align the assembled bkp contigs into the reference genome with blat
#########################################################################
# It will produce a psl with the blat alignments for the assembled contigs 
###########################################################################
# for each cluster insertion
##############################
# Psl naming convention:
########################### 
# - ${blatPath}/${family}:${chr}_${beg}_${end}:${orientation}.psl
# where:
#	- family: TE family (L1, ALU, SVA...)
#	- chr: insertion chromosome
# 	- beg: insertion beginning
#	- end: insertion end
#	- orientation: cluster (+ or -)

step="BLAT"
startTime=$(date +%s)
printHeader "Aligning the assembled bkp contigs into the reference genome with blat"  

## 3.1 Pool all the fasta in a single one
# Ouput:
# - allContigsPath=${blatDir}/"allContigs.fa"

allContigsPath=${blatDir}/"allContigs.fa"
echo -n "" > $allContigsPath

log "1. Producing a single fasta with the contigs from all insertions\n" $step
ls $fastaDir | grep '.*fa' | while read bkpFasta; 
do 	
	bkpId=${bkpFasta%.fa}	
	contigPath=${contigsDir}/${bkpId}".contigs.fa"
		
	awk -v bkpId=$bkpId '! /^>/{row=$0; print $0} /^>/{sub(/>/, "", $1); row=">"bkpId"::"$1; print row;}' $contigPath >> $allContigsPath
done 

## 3.2 Align the contigs with Blat into the reference genome and consensus L1 sequence
# Ouput:
# - blatPath=${blatDir}/"allContigs.psl"

## Notes about blat alignment configuration:
# Default blat configuration does not work well for 5-prime informative contigs (those spanning TE - genomic dna bkp)
# As the TE piece of sequence was mapping multiple times  with better score the genomic dna alignment was not reported.
# Problem solved decreasing the -repMatch value:
# -repMatch=N    Sets the number of repetitions of a tile allowed before
#                it is marked as overused.  Typically this is 256 for tileSize
#                12, 1024 for tile size 11, 4096 for tile size 10.

blatPath=${blatDir}/"allContigs.psl"

log "2. Align the contigs with Blat into the reference genome (and consensus L1 sequence)\n" $step
run "blat -t=dna -q=dna -stepSize=5 tileSize=11 -minScore=20 -repMatch=256 -out=psl -noHead $genome $allContigsPath $blatPath  >> $logFile" "$ECHO"

## 3.3 Split blat output in a single file per insertion and cluster
# Output:
# - a psl for each cluster and insertion
# - ${blatPath}/${family}:${chr}_${beg}_${end}:${orientation}.psl
log "3. Split blat output in a single file per insertion and cluster\n" $step
awk -v OFS="\t" -v outDir=${blatDir} '{split($10,id,"::"); $10=id[2]; outFile=outDir"/"id[1]".psl"; print $0 >> outFile; close(outFile)}' $blatPath

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


# 4) TE insertion breakpoints analysis from assembled and aligned contigs
##########################################################################
# It will produce a single file containing...
##############################################
# $outDir 

## 4.1 Make a list with the TE insertion ids:
# Output:
# - $outDir/insertions_list.txt
insertionsList=$outDir/insertions_list.txt

ls $blatDir | grep '.*psl' | grep -v "allContigs"| awk '{split($1,a,":"); print a[1]":"a[2];}' | sort | uniq > $insertionsList

## 4.2 Prepare input file for insertion breakpoint analysis
# Output:
# - $outDir/paths2bkpAnalysis.txt
paths2bkpAnalysis=$outDir/paths2bkpAnalysis.txt
echo -n "" > $paths2bkpAnalysis
 
cat $insertionsList | while read insertionId; 
do 	
	contigPlusPath=${contigsDir}/${insertionId}:+.contigs.fa
	contigMinusPath=${contigsDir}/${insertionId}:-.contigs.fa
	blatPlusPath=${blatDir}/${insertionId}:+.psl
	blatMinusPath=${blatDir}/${insertionId}:-.psl
	
	printf ${insertionId}"\t"${contigPlusPath}","${contigMinusPath}"\t"${blatPlusPath}","${blatMinusPath}"\n" >> $paths2bkpAnalysis
done


######################
# 4) CLEANUP AND END #
######################

end=$(date +%s)
printHeader "TE-bkpAssembler for $sampleId completed in $(echo "($end-$start)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0

