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

while getopts ":i:f:s:k:o:h" opt "$@"; 
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
if [[ ! -e $fasta ]]; then log "The fasta with the TE insertion supporting reads does not exist. Mandatory argument -e" "ERROR" >&2; usageDoc; exit -1; fi
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
contigsDir=$outDir/Contigs
assemblyLogsDir=$contigsDir/logs

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

# 1) Produces for each TE insertion two fasta (one for read pairs supporting + cluster and another one for read pairs supporting - cluster)
############################################################################################################################################
# These fasta files will be used to assemble the 5' and 3' TE insertion breakpoint sequences 
#############################################################################################
# Each fasta will be named according to this convention: 
########################################################
# - ${sampleId}:${family}:${chr}_${beg}_${end}:${orientation}.txt
# where:
# 	- sampleId: identifier to name log, intermediate and output files
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

# 2) Assemble the 5' and 3' TE insertion breakpoints with velvet. 
#################################################################
# An independent assembly is performed for each insertion breakpoint
#####################################################################
# output: 
############

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
	run "mv ${contigsDir}/Log ${contigsDir}/logs/${bkpId}.log" "$ECHO"
	run "mv ${contigsDir}/contigs.fa ${contigsDir}/${bkpId}.contigs.fa" "$ECHO"
	log "4. Cleaning\n" $step
	run "rm ${contigsDir}/Sequences ${contigsDir}/Roadmaps ${contigsDir}/PreGraph ${contigsDir}/stats.txt ${contigsDir}/LastGraph ${contigsDir}/Graph2" "$ECHO"
done


######################
# 3) CLEANUP AND END #
######################

end=$(date +%s)
printHeader "TE-bkpAssembler for $sampleId completed in $(echo "($end-$start)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0

