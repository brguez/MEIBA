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

    $0 -i <insertions> -f <insertions_fasta> -g <genome> -d <driver_db> -s <sample_identifier> [OPTIONS]

*** MANDATORY

    -i     <TXT>              TraFiC MEI somatic insertion calls for a given sample.

    -f     <FASTA>            Fasta containing MEI insertions supporting reads.

    -g     <FASTA>            Reference Genome in fasta format (RG). Please make sure you provide the same RG version you used to run TraFiC.
                              Also, make sure the same chromosome naming conventions are used.

    -d     <BED>              Database of repetitive sequences according to RepeatMasker in BED format.


    -s     <STRING>           Sample id. Output file will be named accordingly.

*** [OPTIONS] can be:
* General:

    -k     <INTEGER>          K-mer length of the words being hashed for the assembly. Default: N=21.

    -e     <INTEGER>          Minimum MEI score to pass the filtering. Default: N=2.

    -o     <PATH>             Output directory. Default current working directory.

    -h                        Display usage information.


help
}

# Function 2. Parse user's input
################################
function getoptions {

while getopts ":i:f:g:d:s:k:e:o:h" opt "$@";
do
   case $opt in

      ## MANDATORY ARGUMENTS
      i)
      if [ -n "$OPTARG" ];
      then
              insertions=$OPTARG
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

      e)
      if [ -n "$OPTARG" ];
      then
              minScore=$OPTARG
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
            printf "[${label}] $string\n"
        else
            printf "$string\n"
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
version=0.4.0

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
  log "Path not accessible to the script" "ERROR"
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

## Check that input files are ok:
if [[ ! -s $insertions ]]; then log "The TraFiC MEI insertion calls file does not exist or is empty. Mandatory argument -i" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -s $fasta ]]; then log "The MEI insertion supporting reads fasta file does not exist or is empty. Mandatory argument -f" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -s $genome ]]; then log "The reference genome fasta file does not not exist or is empty. Mandatory argument -g" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -s $repeatsDb ]]; then log "The RepeatMasker repeats database does not exist or is empty. Mandatory argument -d" "ERROR" >&2; usageDoc; exit -1; fi
if [[ $sampleId == "" ]]; then log "Sample id does not provided. Mandatory argument -s" "ERROR" >&2; usageDoc; exit -1; fi

## Optional arguments
## ~~~~~~~~~~~~~~~~~~

# K-mer length for assembly
if [[ "$kmerLen" == "" ]];
then
    kmerLen='21';
fi

# minimum MEI score
if [[ "$minScore" == "" ]];
then
    minScore='2';
fi

# Output directory
if [[ "$outDir" == "" ]];
then
    outDir=${SGE_O_WORKDIR-$PWD};
else
    if [[ ! -e "$outDir" ]];
    then
        log "Your output directory does not exist. Option -o" "ERROR" >&2;
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
awkDir=$srcDir/awk

## references:
refDir=$rootDir/ref

## Output files directories
logsDir=$outDir/Logs
fastaDir=$outDir/Fasta
contigsDir=$outDir/Contigs
blatDir=$outDir/Blat
bkpAnalysisDir=$outDir/BkpAnalysis
annotDir=$outDir/Annot
filterDir=$outDir/Filter
srcRegDir=$outDir/SrcRegions

# 5. Scripts/references
########################
# scripts
CLUSTERS2FASTA=$pyDir/clusters2fasta.py
ALIGN_CONTIGS=$bashDir/alignContigs2reference.sh
ADD_INFO=$awkDir/addInfo2insertionList.awk
BKP_ANALYSIS=$pyDir/insertionBkpAnalysis.py
ANNOTATOR=$bashDir/variants_annotator.sh
FILTER=$pyDir/filterVCF.MEI.py

# references
driverDb=$refDir/cancerGenes_COSMIC_CPG.tsv
germlineMEIdb=$refDir/germline_MEI_1KGENOMES.bed
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
printf "  %-34s %s\n" "insertions:" "$insertions"
printf "  %-34s %s\n" "fasta:" "$fasta"
printf "  %-34s %s\n" "genome:" "$genome"
printf "  %-34s %s\n" "repeats-db:" "$repeatsDb"
printf "  %-34s %s\n\n" "sampleId:" "$sampleId"
printf "  %-34s %s\n" "***** OPTIONAL ARGUMENTS *****"
printf "  %-34s %s\n" "K-mer length:" "$kmerLen"
printf "  %-34s %s\n" "min-score:" "$minScore"
printf "  %-34s %s\n\n" "outDir:" "$outDir"


##########
## START #
##########
header="Executing TEIBA $version for $sampleId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
start=$(date +%s)

########################
# I) PRELIMINARY STEPS #
########################

# Create log directory
if [[ ! -d $logsDir ]]; then mkdir $logsDir; fi

# 0) For each orphan transduction (TD2) event, extract source region +/- $windowSize to use as target in BLAT alignment.
########################################################################################################################
# Each fasta will be named according to this convention:
########################################################
# - ${srcRegDir}/${family}:${type}:${chr}_${beg}_${end}.src.fa
# where:
#    - MEI: unique identifier
#    - family: MEI family (L1, ALU, SVA...)
#    - type: td0 (solo-insertion), td1 (partnered-transduccion) and td2 (orphan-transduction)
#    - chr: insertion chromosome
#    - beg: insertion beginning
#    - end: insertion end

## parameters
windowSize=1000
readLen=100

## Make source region fasta directory:
if [[ ! -d $srcRegDir ]]; then mkdir $srcRegDir; fi

## Execute the step
step="SOURCES2FASTA"
startTime=$(date +%s)
printHeader "Prepare fasta for assembly"
log "Extracting region downstream of source element for orphan transductions" $step

cat $insertions | while read chrP begP endP nReadsP famP readsP chrM begM endM nReadsM famM readsM tdType chrSrc begSrc endSrc strSrc begTd endTd lenRna lenTd; do
    if [[ $tdType == "TD2" ]]; then
        endP=`expr $endP + $readLen`
        # NOTE: insertion IDs used here correspond to format generated by $CLUSTERS2FASTA
        targetRegionFile=${famP}":"${tdType}":"${chrP}"_"${endP}"_"${begM}".src.fa"
        targetBeg=`expr $begTd - $windowSize`
        targetEnd=`expr $endTd + $windowSize`
        if [ "$targetBeg" -lt 0 ]; then targetBeg=0; fi # Set lower-bound to 0 (avoid negative coordinates)
        targetInterval=$chrSrc":"$targetBeg"-"$targetEnd
        echo "samtools faidx $genome $targetInterval > $srcRegDir/$targetRegionFile" >&1
        samtools faidx $genome $targetInterval > $srcRegDir/$targetRegionFile
    fi
done

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


##################
# II) CORE STEPS #
##################

# 1) Produces for each Mobile Element Insertion (MEI) two fasta (one for read pairs supporting + cluster and another one for read pairs
##########################################################################################################################################
# supporting - cluster) These fasta files will be used to assemble the 5' and 3' MEI insertion breakpoint sequences
###################################################################################################################
# Each fasta will be named according to this convention:
########################################################
# - ${fastaDir}/${family}:${type}:${chr}_${beg}_${end}:${orientation}.fa
# where:
#    - family: MEI family (L1, ALU, SVA...)
#    - type: td0 (solo-insertion), td1 (partnered-transduccion) and td2 (orphan-transduction)
#    - chr: insertion chromosome
#    - beg: insertion beginning
#    - end: insertion end
#    - orientation: cluster (+ or -)

## Make fasta for assembly directory:
if [[ ! -d $fastaDir ]]; then mkdir $fastaDir; fi

## Execute the step
step="CLUSTERS2FASTA"
startTime=$(date +%s)
printHeader "Prepare fasta for assembly"
log "Producing per MEI two fasta for insertion bkp assembly" $step
run "python $CLUSTERS2FASTA $insertions $fasta $genome --outDir $fastaDir 1> $logsDir/1_clusters2fasta.out 2> $logsDir/1_clusters2fasta.err" "$ECHO"
endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

# 2) Assemble the 5' and 3' MEI breakpoints with velvet. An independent assembly is performed for each insertion breakpoint
######################################################################################################################################
# It will produce a fasta containing the assembled contigs for each cluster insertion
######################################################################################
# Fasta naming convention:
###########################
# - ${contigsDir}/${family}:${type}:${chr}_${beg}_${end}:${orientation}.contigs.fa
# where:
#    - family: MEI family (L1, ALU, SVA...)
#    - type: td0 (solo-insertion), td1 (partnered-transduccion) and td2 (orphan-transduction)
#    - chr: insertion chromosome
#    - beg: insertion beginning
#    - end: insertion end
#    - orientation: cluster (+ or -)

## Make assembly directory:
if [[ ! -d $contigsDir ]]; then mkdir $contigsDir; fi

## Execute the step
step="BKP-ASSEMBLY"
startTime=$(date +%s)
printHeader "Assembling the 5' and 3' MEI breakpoints with velvet"

ls $fastaDir | grep '.*fa' | while read bkpFasta;
do
    bkpId=${bkpFasta%.fa}
    fastaPath=${fastaDir}/${bkpFasta}
    contigPath=${contigsDir}/${bkpId}".contigs.fa"

    log "** ${bkpId} breakpoint **" $step
    log "1. Preparing files for assembly" $step
    run "velveth $contigsDir $kmerLen -fasta -short $fastaPath 1>> $logsDir/2_assembly.out 2>> $logsDir/2_assembly.err" "$ECHO"
    log "2. Breakpoint assembly with velvet" $step
    run "velvetg $contigsDir -exp_cov auto -cov_cutoff auto 1>> $logsDir/2_assembly.out 2>> $logsDir/2_assembly.err" "$ECHO"
    log "3. Rename output files" $step
    run "mv ${contigsDir}/contigs.fa ${contigsDir}/${bkpId}.contigs.fa" "$ECHO"
    log "4. Second Cleaning" $step
#    run "rm $fastaPath ${contigsDir}/Log ${contigsDir}/Sequences ${contigsDir}/Roadmaps ${contigsDir}/PreGraph ${contigsDir}/stats.txt ${contigsDir}/LastGraph ${contigsDir}/Graph2" "$ECHO"

done

## Remove temporary fasta directory
#rm -r $fastaDir

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


# 3) Align the assembled bkp contigs into the reference genome with blat
#########################################################################
# It will produce a psl with the blat alignments for the assembled contigs
###########################################################################
# for each cluster insertion
############################
# Psl naming convention:
########################
# - ${blatPath}/${family}:${type}:${chr}_${beg}_${end}:${orientation}.psl
# where:
#    - family: MEI family (L1, ALU, SVA...)
#    - type: td0 (solo-insertion), td1 (partnered-transduccion) and td2 (orphan-transduction)
#    - chr: insertion chromosome
#    - beg: insertion beginning
#    - end: insertion end
#    - orientation: cluster (+ or -)

## Make blat directory:
if [[ ! -d $blatDir ]]; then mkdir $blatDir; fi

## Execute the step
step="BLAT"
startTime=$(date +%s)

printHeader "Aligning the assembled bkp contigs into the reference genome with blat"

ls $contigsDir | grep '.*fa' | while read bkpContig;
do
    bkpContigPath=${contigsDir}/${bkpContig}
    bkpId=${bkpContig%.contigs.fa}
    read family tdType <<<$(echo $bkpId | awk '{split($1, info, ":"); family=info[1]; tdType=info[2]; print family, tdType;}')

    # For each insertion, select the corresponding consensus MEI sequence to align the contigs into
    case $family in

        L1)
            consensusMEI=$consensusL1
        ;;

        Alu)
            consensusMEI=$consensusAlu
        ;;

        SVA)
            consensusMEI=$consensusSVA
        ;;

        ERVK)
            consensusMEI=$consensusERVK
        ;;
    esac

    # When processing an orphan transduction (TD2) event,
    #   use genomic DNA downstream of source element
    #   otherwise, use consensus MEI sequence
    # as BLAT target.
    if [[ $tdType == "TD2" ]]; then
        # remove +/- from breakpoint ID to get source region filename
        srcTarget=$srcRegDir/${bkpId/:[+-]/}".src.fa"
    else
        srcTarget=$consensusMEI
    fi

    # Align contigs into the insertion target region and MEI sequence
    log "** ${bkpId} breakpoint **" $step
    log "1. Blat alignment" $step
    run "bash $ALIGN_CONTIGS $bkpContigPath $bkpId $genome $srcTarget 1000 $blatDir 1>> $logsDir/3_blat.out 2>> $logsDir/3_blat.err" "$ECHO"
done

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


# 4) MEI breakpoint analysis from assembled and aligned contigs
#####################################################################
# It will produce a single file containing the insertion breakpoint and
########################################################################
# many pieces of information associated (orientation, TSD, length...)
# ####################################################################
# - $outDir/TEIBA.results.txt

## Make bkp analysis directory:
if [[ ! -d $bkpAnalysisDir ]]; then mkdir $bkpAnalysisDir; fi

### Execute the step
## 4.1 Make a list with the MEI ids:
# Output:
# - $outDir/insertions_list.txt
insertionList=$bkpAnalysisDir/insertionList.txt

ls $blatDir | grep '.*psl' | grep -v "allContigs"| awk '{split($1,a,":"); print a[1]":"a[2]":"a[3];}' | sort | uniq > $insertionList

## 4.2 For each insertion add the list of read pairs supporting the clusters, the source element and transduction info (only for transductions)
# Output:
# - $bkpAnalysisDir/insertionList_plusInfo.txt
insertionListInfo=$bkpAnalysisDir/insertionList_plusInfo.txt

awk -v OFS='\t' -v fileRef=$insertions -f $ADD_INFO $insertionList > $insertionListInfo

# Remove intermediate files:
# rm $insertionListInfo

## 4.3 Prepare input file for insertion breakpoint analysis
# Output:
# - $outDir/paths2bkpAnalysis.txt
paths2bkpAnalysis=$outDir/paths2bkpAnalysis.txt
echo -n "" > $paths2bkpAnalysis

cat $insertionListInfo | while read insertionId readPairsPlus readPairsMinus sourceElementInfo transductionInfo;
do
    contigPlusPath=${contigsDir}/${insertionId}:+.contigs.fa
    contigMinusPath=${contigsDir}/${insertionId}:-.contigs.fa
    blatPlusPath=${blatDir}/${insertionId}:+.psl
    blatMinusPath=${blatDir}/${insertionId}:-.psl

    printf ${insertionId}"\t"${contigPlusPath}","${contigMinusPath}"\t"${blatPlusPath}","${blatMinusPath}"\t"${readPairsPlus}"\t"${readPairsMinus}"\t"${sourceElementInfo}"\t"${transductionInfo}"\n" >> $paths2bkpAnalysis
done

# Remove intermediate files:
#rm $insertionListSupReads

## 4.3 Perform breakpoint analysis
# Output:
# - $bkpAnalysisDir/$sampleId.vcf
rawVCF=$bkpAnalysisDir/$sampleId.vcf

if [ ! -s $rawVCF ];
then
    step="BKP-ANALYSIS"
    startTime=$(date +%s)
    printHeader "Performing MEI breakpoint analysis"
    log "Identifying insertion breakpoints, TSD, MEI length, orientation and structure" $step
    run "python $BKP_ANALYSIS $paths2bkpAnalysis $sampleId $genome --outDir $bkpAnalysisDir 1>> $logsDir/4_bkpAnalysis.out 2>> $logsDir/4_bkpAnalysis.err" "$ECHO"

    if [ ! -s $rawVCF ];
    then
        log "Error performing breakpoint analysis" "ERROR"
            exit -1
    else
        endTime=$(date +%s)
        printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    fi
else
    printHeader "Output file already exists... skipping step"
fi

## Remove temporary contigs and blat directories
#rm $paths2bkpAnalysis
#rm -r $contigsDir $blatDir


# 5) Annotate MEI
###################
## Output:
# -  $annotDir/$sampleId.annotated.vcf
annotVCF=$annotDir/$sampleId.annotated.vcf

## Make MEI annotation directory:
if [[ ! -d $annotDir ]]; then mkdir $annotDir; fi

### Execute the step
if [ ! -s $annotVCF ];
then
    step="ANNOTATION"
    startTime=$(date +%s)
    printHeader "Performing MEI breakpoint annotation"
    log "Annotating MEI" $step
    run "bash $ANNOTATOR $rawVCF $repeatsDb $driverDb $germlineMEIdb $sampleId $annotDir 1>> $logsDir/5_annotation.out 2>> $logsDir/5_annotation.err" "$ECHO"

    if [ ! -s $annotVCF ];
    then
        log "Error performing annotation" "ERROR"
            exit -1
    else
        endTime=$(date +%s)
        printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    fi
else
    printHeader "Output file already exists... skipping step"
fi

## Remove temporary bkp analysis directory
#rm -r $bkpAnalysisDir


# 6) Filter MEI
#################
## Output:
# -  $filterDir/$sampleId.filtered.vcf
filteredVCF=$filterDir/$sampleId.filtered.vcf

## Make MEI filtering directory:
if [[ ! -d $filterDir ]]; then mkdir $filterDir; fi

### Execute the step
# NOTE: for germline variants use a minimum score of 5...

if [ ! -s $filteredVCF ];
then
    step="FILTER"
    startTime=$(date +%s)
    printHeader "Performing MEI filtering"
    log "Filtering MEI" $step
    run "$FILTER $annotVCF $sampleId --min-score $minScore --min-score-ERVK 5 --max-divergence 300 --outDir $filterDir 1>> $logsDir/6_filter.out 2>> $logsDir/6_filter.err" "$ECHO"

    if [ ! -s $filteredVCF ];
    then
        log "Error performing filtering" "ERROR"
            exit -1
    else
        endTime=$(date +%s)
        printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    fi
else
    printHeader "Output file already exists... skipping step"
fi

## Remove temporary annotation directory
#rm -r $annotDir

#############################
# 7) MAKE OUTPUT VCF AND END #
#############################

## Produce output VCF
finalVCF=$outDir/$sampleId.vcf

cp $filteredVCF $finalVCF

## Remove temporary filter directory
#rm -r $filterDir

## End
end=$(date +%s)
printHeader "TEIBA for $sampleId completed in $(echo "($end-$start)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0
