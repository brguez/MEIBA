#!/bin/bash


<<authors
******************************************************************************
	
	TEIBA.sh

	Transposable Element Insertion Breakpoint Analyzer (TEIBA)
	
	Copyright (c) 2016 Bernardo Rodríguez-Martín
	
	Mobile Genomes & Disease Lab.
	Universidad de Vigo (Spain)

	Licenced under the GNU General Public License 3.0 license.
******************************************************************************
authors


# Function 1. Print basic usage information
############################################
function usageDoc
{
cat <<help
	

**** TEIBA version $version ****
Execute for one dataset (sample).
	
*** USAGE

	$0 -i <TraFic_insertions> -f <FASTA> -s <sample_identifier> [OPTIONS]

*** MANDATORY 
		
	-i	<TXT>			TraFiC TE somatic insertion calls for a given sample.		

	-f	<FASTA>			Fasta containing TE insertions supporting reads. 

	-g 	<FASTA>			Reference Genome in fasta format (RG). Please make sure you provide the same RG version you used to run TraFiC. 
					Also, make sure the same chromosome naming conventions are used.
	
	-d      <BED>			Database of repetitive sequences according to RepeatMasker in BED format. 

	
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

while getopts ":i:f:g:d:s:k:o:h" opt "$@"; 
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

      d)
	  if [ -n "$OPTARG" ];
	  then
              repeatsDb=$OPTARG
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

# TEIBA version 
version=0.1

# Enable extended pattern matching 
shopt -s extglob

# Exit if error
set -e -o pipefail

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
if [[ ! -e $repeatsDb ]]; then log "The RepeatMasker repeats database does not exist. Mandatory argument -d" "ERROR" >&2; usageDoc; exit -1; fi 

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
binDir=$rootDir/bin
srcDir=$rootDir/src
pyDir=$srcDir/python
bashDir=$srcDir/bash 

## references:
refDir=$rootDir/ref

## Output files directories
logsDir=$outDir/Logs
fastaDir=$outDir/Fasta
contigsDir=$outDir/Contigs
assemblyLogsDir=$logsDir/2_assembly
blatDir=$outDir/Blat
bkpAnalysisDir=$outDir/BkpAnalysis
annotDir=$outDir/Annot

# 5. Scripts/references
########################
# scripts
CLUSTERS2FASTA=$pyDir/clusters2fasta.py
ALIGN_CONTIGS=$bashDir/alignContigs2reference.sh
BKP_ANALYSIS=$pyDir/insertionBkpAnalysis.py
ANNOTATOR=$bashDir/variants_annotator.sh

# references
driverDb=$refDir/driverDb_COSMIC_CPG.20162706.tsv
consensusL1=$refDir/L1_consensus.fa 
consensusAlu=$refDir/Alu_consensus.fa
consensusSVA=$refDir/SVA_consensus.fa
consensusERVK=$refDir/ERVK_consensus.fa

## DISPLAY PROGRAM CONFIGURATION  
##################################
printf "\n"
header=" TEIBA CONFIGURATION FOR $sampleId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n" "***** MANDATORY ARGUMENTS *****"
printf "  %-34s %s\n" "input:" "$input"
printf "  %-34s %s\n" "fasta:" "$fasta"
printf "  %-34s %s\n" "genome:" "$genome"
printf "  %-34s %s\n" "repeats-db:" "$repeatsDb"
printf "  %-34s %s\n\n" "sampleId:" "$sampleId"
printf "  %-34s %s\n" "***** OPTIONAL ARGUMENTS *****"
printf "  %-34s %s\n" "K-mer length:" "$kmerLen"
printf "  %-34s %s\n\n" "outDir:" "$outDir"


##########
## START #
##########
header="Executing TEIBA $version for $sampleId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
start=$(date +%s)

#######################    	
# 0) PRELIMINARY STEPS #
#######################

## 0.1) Make directories
#######################
# Fastas for assembly
if [[ ! -d $fastaDir ]]; then mkdir $fastaDir; fi

# Assembled contigs of TE insertion breakpoints and assembly logs
if [[ ! -d $contigsDir ]]; then mkdir $contigsDir; fi

# Contigs alignment with Blat:
if [[ ! -d $blatDir ]]; then mkdir $blatDir; fi

# MEI breakpoint analysis directory
if [[ ! -d $bkpAnalysisDir ]]; then mkdir $bkpAnalysisDir; fi

# MEI annotation directory:
if [[ ! -d $annotDir ]]; then mkdir $annotDir; fi

# Log directories:
if [[ ! -d $logsDir ]]; then mkdir $logsDir; fi
if [[ ! -d $assemblyLogsDir ]]; then mkdir $assemblyLogsDir; fi


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
printHeader "Prepare fasta for assembly"  
log "Producing per TE insertion two fasta for insertion bkp assembly\n" $step
run "python $CLUSTERS2FASTA $input $fasta --outDir $fastaDir 1> $logsDir/1_clusters2fasta.out 2> $logsDir/1_clusters2fasta.err" "$ECHO"	
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
	run "velveth $contigsDir $kmerLen -fasta -short $fastaPath 1>> $logsDir/2_assembly.out 2>> $logsDir/2_assembly.err" "$ECHO"
	log "2. Breakpoint assembly with velvet\n" $step
	run "velvetg $contigsDir -exp_cov auto -cov_cutoff auto 1>> $logsDir/2_assembly.out 2>> $logsDir/2_assembly.err" "$ECHO"
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

ls $contigsDir | grep '.*fa' | while read bkpContigs;
do
	bkpContigsPath=${contigsDir}/${bkpContigs}	
	bkpId=${bkpContigs%.contigs.fa}
	category=`echo $bkpId | awk '{split($1,a,":"); print a[1];}'` 	

	# For each insertion, select the corresponding consensus TE sequence to align the contigs into
	case $category in   	
            
	    L1)
              	consensusTE=$consensusL1
	       	;;
      
            Alu)
                consensusTE=$consensusAlu
	       	;;

	    SVA)
		consensusTE=$consensusSVA	
		;;

 	    ERVK)
		consensusTE=$consensusERVK
		;;
	esac
	
	# Align contigs into the insertion target region and TE sequence
	log "** ${bkpId} breakpoint **\n" $step
	run "bash $ALIGN_CONTIGS $bkpContigsPath $bkpId $genome $consensusTE 1000 $blatDir 1>> $logsDir/3_blat.out 2>> $logsDir/3_blat.err" "$ECHO"
done 

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


# 4) TE insertion breakpoint analysis from assembled and aligned contigs
##########################################################################
# It will produce a single file containing the insertion breakpoint and 
########################################################################
# many pieces of information associated (orientation, TSD, length...)
# ####################################################################
# - $outDir/TEIBA.results.txt 

## 4.1 Make a list with the TE insertion ids:
# Output:
# - $outDir/insertions_list.txt
insertionsList=$bkpAnalysisDir/insertions_list.txt

ls $blatDir | grep '.*psl' | grep -v "allContigs"| awk '{split($1,a,":"); print a[1]":"a[2];}' | sort | uniq > $insertionsList

## 4.2 Prepare input file for insertion breakpoint analysis
# Output:
# - $outDir/paths2bkpAnalysis.txt
paths2bkpAnalysis=$bkpAnalysisDir/paths2bkpAnalysis.txt
echo -n "" > $paths2bkpAnalysis
 
cat $insertionsList | while read insertionId; 
do 	
	contigPlusPath=${contigsDir}/${insertionId}:+.contigs.fa
	contigMinusPath=${contigsDir}/${insertionId}:-.contigs.fa
	blatPlusPath=${blatDir}/${insertionId}:+.psl
	blatMinusPath=${blatDir}/${insertionId}:-.psl
	
	printf ${insertionId}"\t"${contigPlusPath}","${contigMinusPath}"\t"${blatPlusPath}","${blatMinusPath}"\n" >> $paths2bkpAnalysis
done

## 4.3 Perform breakpoint analysis
# Output:
# - $bkpAnalysisDir/$sampleId.vcf 
rawVCF=$bkpAnalysisDir/$sampleId.vcf 

if [ ! -s $rawVCF ]; 
then
	step="BKP-ANALYSIS"
	startTime=$(date +%s)
	printHeader "Performing MEI breakpoint analysis"
	log "Identifying insertion breakpoints, TSD, TE length, TE orientation and TE structure\n" $step  
	run "python $BKP_ANALYSIS $paths2bkpAnalysis $sampleId $genome --outDir $bkpAnalysisDir 1>> $logsDir/4_bkpAnalysis.out 2>> $logsDir/4_bkpAnalysis.err" "$ECHO"
	
	if [ ! -s $rawVCF ]; 
	then	
		log "Error performing breakpoint analysis\n" "ERROR" 
        	exit -1
	else
		endTime=$(date +%s)
		printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	fi
else
	printHeader "Output file already exists... skipping step"
fi


# 5) Annotate mobile elements insertions
##########################################
## Output:
# -  $outDir/$sampleId.annotated.vcf
annotVCF=$annotDir/$sampleId.annotated.vcf 

if [ ! -s $annotVCF ]; 
then
	step="ANNOTATION"
	startTime=$(date +%s)
	printHeader "Performing MEI breakpoint annotation"
	log "Annotating MEI\n" $step  
	run "bash $ANNOTATOR $rawVCF $repeatsDb $driverDb $sampleId $annotDir 1>> $logsDir/5_annotation.out 2>> $logsDir/5_annotation.err" "$ECHO"
	
	if [ ! -s $annotVCF ]; 
	then	
		log "Error performing annotation\n" "ERROR" 
        	exit -1
	else
		endTime=$(date +%s)
		printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	fi
else
	printHeader "Output file already exists... skipping step"
fi


######################################
# 6) MAKE OUTPUT GZ, CLEANUP AND END #
######################################

## Produce a compressed vcf as output
finalVCF=$outDir/$sampleId.vcf

cp $annotVCF $finalVCF

## Cleaning
# rm -r $fastaDir $contigsDir $blatDir $bkpAnalysisDir $annotDir

## End
end=$(date +%s)
printHeader "TEIBA for $sampleId completed in $(echo "($end-$start)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0

