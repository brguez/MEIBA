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

##### TEIBA

## TEIBA is aimed to characterize at base pair resolution a set of mobile element insertions provided as input throught local assembly and breakpoint analysis.
# It currently handles four types of insertion events:

# - TD0: solo L1, Alu, SVA and ERVK insertions
# - TD1: L1 partnered transductions
# - TD2: L1 orphan transductions
# - PSD: processed pseudogene insertions

##### DESCRIPTION OF THE FIELDS IN THE INSERTIONS INPUT FILE FORMAT:

### Generic fields (TD0, TD1, TD2 and PSD)
# 1.  chrom_plus_cluster
# 2.  beg_plus_cluster
# 3.  end_plus_cluster
# 4.  nbReads_plus_cluster
# 5.  class_plus_cluster, "NA" for PSD
# 6.  readList_plus_cluster
# 7.  chrom_minus_cluster
# 8.  beg_minus_cluster
# 9.  end_minus_cluster
# 10. nbReads_minus_cluster
# 11. class_minus_cluster, "NA" for PSD
# 12. readList_minus_cluster
# 13. insertion_type (TD0, TD1, TD2 or PSD)

### L1 transductions specific fields (TD1 and TD2). "NA" for TD0 and PSD
# 14. cytobandId	
# 15. sourceType (GERMLINE or SOMATIC)
# 16. chrom_source_element
# 17. beg_source_element
# 18. end_source_element
# 19. orientation_source_element
# 20. transduction_beg
# 21. transduction_end
# 22. transduction_rna_length
# 23. transduction_length

### Processed pseudogene specific fields (PSD). "NA" for TD0, TD1 and TD2
# 24. psd_gene
# 25. chrom_exonA_cluster
# 26. beg_exonA_cluster
# 27. end_exonA_cluster
# 28. chrom_exonB_cluster
# 29. beg_exonB_cluster
# 30. end_exonB_cluster

### L-mediated genomic rearrangement specific fields. "NA" for standard L1, Alu, SVA, ERVK and PSD insertions
# 31. gr_type (DEL, DUP or TRANS)

# Function 1. Print usage information
#######################################
function usageDoc
{
cat <<help

**** TEIBA version $version ****

Execute TEIBA on one dataset (sample).

*** USAGE
    $0 -i <insertions> -f <insertions_fasta> -g <genome> -d <driver_db> --sample-id <sample_identifier> --file-name <output_fileName>[OPTIONS]

*** MANDATORY
    -i|--insertions     <TSV>              TraFiC mobile element insertion (MEI) candidate calls for a given sample.
    -b|--bam            <BAM>              BAM file containing MEI insertions supporting reads.
    -g|--genome         <FASTA>            Reference Genome in fasta format (RG). Make sure the same RG version is used than for running TraFiC.
                                           Also, make sure the same chromosome naming conventions are used.
    -d|--repeats-db     <BED>              Database of repetitive sequences according to RepeatMasker in BED format.
    --sample-id	        <STRING>           Sample identifier to be incorporated in the SL field of the output VCF. 
                                           In PCAWG is used the normal_wgs_aliquot_id or tumor_wgs_aliquot_id.
    --file-name	        <STRING>           Output VCF name. In PCAWG we use the submitted_donor_id.

*** [OPTIONS] can be:
* General:
    -o|--output-dir     <PATH>             Output directory. Default current working directory.
    --tmp-dir		<PATH>		   Temporary directory. Default /tmp.
    --no-annotation	                   Skip MEI annotation step.
    --no-cleanup	                   Keep intermediate files.
    -h|--help			           Display usage information


* Filters:
    --filters           <(FILTER_1)>, ... ,<(FILTER_N)>	List of filters to be applied out of 6 possible filtering criteria: 
                                                        NONE, SCORE, REP, DUP, FPSOURCE and GERMLINE. 'NONE' will disable any filtering.
                                                        Default='SCORE,DUP,FPSOURCE'
    --score-L1-TD0      <INTEGER>                       Minimum assembly score for solo L1 insertions. Default 2.
    --score-L1-TD1      <INTEGER>                       Minimum assembly score for L1 partnered transductions. Default 2.
    --score-L1-TD2      <INTEGER>                       Minimum assembly score for L1 orphan transductions. Default 2.
    --score-Alu         <INTEGER>                       Minimum assembly score for Alu insertions. Default 2.
    --score-SVA         <INTEGER>                       Minimum assembly score for SVA insertions. Default 2.
    --score-ERVK        <INTEGER>                       Minimum assembly score for ERVK insertions. Default 2.
    --score-PSD         <INTEGER>                       Minimum assembly score for processed-pseudogene (PSD) insertions. Default 2.

* Files:
    --germline-VCF      <VCF>                           VCF with germline MEI calls for a given donor. If provided, input insertions are considered to be somatic.
                                                        Necesary for GERMLINE filtering.
help
}


# Function 2. Parse user's input
################################
function getoptions {

ARGS=`getopt -o "i:b:g:d:o:h" -l "insertions:,bam:,genome:,repeats-db:,sample-id:,file-name:,output-dir:,tmp-dir:,no-annotation,no-cleanup,help,filters:,score-L1-TD0:,score-L1-TD1:,score-L1-TD2:,score-Alu:,score-SVA:,score-ERVK:,score-PSD:,score-GR:,germline-VCF:" \
      -n "$0" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi

# A little magic
eval set -- "$ARGS"

while true;
do

    case "$1" in

        ## MANDATORY ARGUMENTS
        -i|--insertions)
            if [ -n "$2" ];
            then
                insertions=$2
            fi
            shift 2;;

        -b|--bam)
            if [ -n "$2" ];
            then
                bam=$2
            fi
            shift 2;;

        -g|--genome)
            if [ -n "$2" ];
            then
                genome=$2
            fi
            shift 2;;

        -d|--repeats-db)
            if [ -n "$2" ];
            then
                repeatsDb=$2
            fi
            shift 2;;

        --sample-id)
            if [ -n "$2" ];
            then
                sampleId=$2
            fi
            shift 2;;

        --file-name)
            if [ -n "$2" ];
            then
                fileName=$2
            fi
            shift 2;;

        ## OPTIONS
        # General:
        -o|--output-dir)
            if [ -n "$2" ];
            then
                outDir=$2
            fi
            shift 2;;

        --tmp-dir)
      	  if [ -n $2 ];
      	  then
              TMPDIR=$2
      	  fi
      	  shift 2;;

        --no-annotation)
            annot="FALSE";
            shift;;

        --no-cleanup)
            cleanup="FALSE";
            shift;;

        -h|--help)
            usagedoc;
            exit 1
            shift;;

        # Filters:
        --filters)
            if [ -n "$2" ];
            then
                filterList=$2
            fi
            shift 2;;

        --score-L1-TD0)
            if [ -n "$2" ];
            then
                scoreL1_TD0=$2
            fi
            shift 2;;

        --score-L1-TD1)
            if [ -n "$2" ];
            then
                scoreL1_TD1=$2
            fi
            shift 2;;

        --score-L1-TD2)
            if [ -n "$2" ];
            then
                scoreL1_TD2=$2
            fi
            shift 2;;

        --score-Alu)
            if [ -n "$2" ];
            then
                scoreAlu=$2
            fi
            shift 2;;

        --score-SVA)
            if [ -n "$2" ];
            then
                scoreSVA=$2
            fi
            shift 2;;

        --score-ERVK)
            if [ -n "$2" ];
            then
                scoreERVK=$2
            fi
            shift 2;;

        --score-PSD)
            if [ -n "$2" ];
            then
                scorePSD=$2
            fi
            shift 2;;

        --score-GR)
            if [ -n "$2" ];
            then
                scoreGR=$2
            fi
            shift 2;;

        # Files:
        --germline-VCF)
            if [ -n "$2" ];
            then
                germlineVCF=$2
            fi
            shift 2;;

        # Exit loop when finish
        --)
            shift
            break;;
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

# Function 7. Extract a region from a reference sequence
#########################################################
function extractRegion {
    ref=$1     # reference FASTA file name
    seq=$2     # target sequence id
    beg=$3     # left-most coordinate
    end=$4     # right-most coordinate
    extend=$5  # number of bp to extend the target region by
    tdType=$6  # insertion type
    outPath=$7 # output file path

    targetBeg=`expr $beg - $extend`
    targetEnd=`expr $end + $extend`
    if [ "$targetBeg" -lt 0 ]; then targetBeg=0; fi # Set lower-bound to 0 (avoid negative coordinates)
    targetInterval=$seq":"$targetBeg"-"$targetEnd
    echo "samtools faidx $ref $targetInterval | sed -E s/>[0-9]+:/>$tdType:/g > $outPath" >&1
    samtools faidx $ref $targetInterval | sed -E "s/>[0-9]+:/>$tdType:/g" > $outPath
}

# Function 8. If cleanup is enabled, make sure the pipeline removes intermediate files even if it fails
########################################################################################################
function cleanupFunc {

    ## Cleanup enabled
    if [[ "$cleanup" == "TRUE" ]]; then
        log "Finish execution. Performing cleanup before exit" "FINISH"

        if [[ -d $srcRegDir ]]; then rm -r $srcRegDir; fi
        if [[ -d $clippedDir ]]; then rm -r $clippedDir; fi
        if [[ -d $blatDir  ]]; then rm -r $blatDir; fi
        if [[ -d $bkpAnalysisDir  ]]; then rm -r $bkpAnalysisDir; fi
        if [[ -d $annotDir  ]]; then rm -r $annotDir; fi
        if [[ -d $filterDir  ]]; then rm -r $filterDir; fi
    else
        log "Finish execution. Don't remove intermediate files as cleanup not enabled" "FINISH"
    fi
}


# SETTING UP THE ENVIRONMENT
############################

# TEIBA version
version=0.7.0

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

#### Mandatory arguments
##   ~~~~~~~~~~~~~~~~~~~

## Check that everything is ok:
if [[ ! -s $insertions ]]; then log "TraFiC MEI calls file does not exist or is empty. Mandatory argument -i|--insertions" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -s $bam ]]; then log "MEI supporting reads BAM file does not exist or is empty. Mandatory argument -b|--bam" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -s $genome ]]; then log "The reference genome fasta file does not not exist or is empty. Mandatory argument -g|--genome" "ERROR" >&2; usageDoc; exit -1; fi
if [[ ! -s $repeatsDb ]]; then log "The RepeatMasker repeats database does not exist or is empty. Mandatory argument -d|--repeats-db" "ERROR" >&2; usageDoc; exit -1; fi
if [[ $sampleId == "" ]]; then log "Sample id does not provided. Mandatory argument --sample-id" "ERROR" >&2; usageDoc; exit -1; fi
if [[ $fileName == "" ]]; then log "Output file name does not provided. Mandatory argument --file-name" "ERROR" >&2; usageDoc; exit -1; fi


#### Optional arguments
##   ~~~~~~~~~~~~~~~~~~

#### General
## Output directory
if [[ "$outDir" == "" ]];
then
	outDir=${SGE_O_WORKDIR-$PWD};
else
	if [[ ! -e "$outDir" ]];
	then
		log "Your output directory does not exist. Option -o|--output-dir\n" "ERROR" >&2;
		usageDoc;
		exit -1;
	fi
fi

# Temporary directory
if [[ "$TMPDIR" == "" ]];
then
	TMPDIR='/tmp';
else
	if [[ ! -e "$TMPDIR" ]];
	then
		log "Your temporary directory does not exist. Option --tmp-dir\n" "ERROR" >&2;
		usagedoc;
		exit -1;
	fi
fi

## Annotation
if [[ "$annot" != "FALSE" ]];
then
    annot="TRUE";
fi

## Clean up
if [[ "$cleanup" != "FALSE" ]];
then
    cleanup='TRUE';
fi


#### Filters:

## List of filters to be applied
if [[ "$filterList" == "" ]];
then
	filterList='SCORE,DUP,FPSOURCE';
fi

## Minimum assembly score for solo L1 insertions
if [[ "$scoreL1_TD0" == "" ]];
then
    scoreL1_TD0=2;
else
	if [[ ! "$scoreL1_TD0" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for solo L1 insertions. Option --score-L1-TD0\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi

## Minimum assembly score for L1 partnered transductions
if [[ "$scoreL1_TD1" == "" ]];
then
    scoreL1_TD1=2;
else
	if [[ ! "$scoreL1_TD1" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for L1 partnered transductions. Option --score-L1-TD1\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi

## Minimum assembly score for L1 orphan transductions
if [[ "$scoreL1_TD2" == "" ]];
then
    scoreL1_TD2=2;
else
	if [[ ! "$scoreL1_TD2" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for L1 orphan transduction. Option --score-L1-TD2\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi

## Minimum assembly score for Alu insertions
if [[ "$scoreAlu" == "" ]];
then
    scoreAlu=2;
else
	if [[ ! "$scoreAlu" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for Alu insertions. Option --score-Alu\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi

## Minimum assembly score for SVA insertions
if [[ "$scoreSVA" == "" ]];
then
    scoreSVA=2;
else
	if [[ ! "$scoreSVA" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for SVA insertions. Option --score-SVA\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi

## Minimum assembly score for ERVK insertions
if [[ "$scoreERVK" == "" ]];
then
    scoreERVK=2;
else
	if [[ ! "$scoreERVK" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for ERVK insertions. Option --score-ERVK\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi

## Minimum assembly score for pseudogene insertions
if [[ "$scorePSD" == "" ]];
then
    scorePSD=2;
else
	if [[ ! "$scorePSD" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for PSD insertions. Option --score-PSD\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi

## Minimum assembly score for L1-mediated rearrangements
if [[ "$scoreGR" == "" ]];
then
    scoreGR=2;
else
	if [[ ! "$scoreGR" =~ ^[0-9]+$ ]];
    then
	    log "Please specify a proper minimum assembly score for L1-mediated rearrangements. Option --score-GR\n" "ERROR" >&2;
        usageDoc;
        exit -1;
    fi
fi



#### Files:

## VCF with germline MEI calls for filtering out GERMLINE insertions miscalled as SOMATIC
if [[ "$germlineVCF" == "" ]]
then
    germlineVCF="NOT_PROVIDED";
fi


# 4. Directories
################
## binaries and scripts
srcDir=$rootDir/src
pyDir=$srcDir/python
bashDir=$srcDir/bash
awkDir=$srcDir/awk


## references:
refDir=$rootDir/ref

## Output files directories
logsDir=$outDir/Logs
clippedDir=$outDir/Clipped
blatDir=$outDir/Blat
bkpAnalysisDir=$outDir/BkpAnalysis
annotDir=$outDir/Annot
filterDir=$outDir/Filter
srcRegDir=$outDir/SrcRegions

# The temporary directory will be exported as an environmental variable since it will
# be used by every TEIBA's scripts
export TMPDIR=$TMPDIR

## make sure the pipeline removes intermediate files even if it fails
trap cleanupFunc EXIT

# 5. Scripts/references
########################
# scripts
EMPTYVCF=$pyDir/makeEmptyVCF.py
CLIPPED=$pyDir/clusterClippedReads.py
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

## DISPLAY PIPELINE CONFIGURATION
##################################

printf "\n"
header="TEIBA CONFIGURATION FOR $sampleId"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n\n" "TEIBA Version $version"
printf "  %-34s %s\n" "***** MANDATORY ARGUMENTS *****"
printf "  %-34s %s\n" "insertions:" "$insertions"
printf "  %-34s %s\n" "bam:" "$bam"
printf "  %-34s %s\n" "genome:" "$genome"
printf "  %-34s %s\n" "repeatsDb:" "$repeatsDb"
printf "  %-34s %s\n" "sample-id:" "$sampleId"
printf "  %-34s %s\n\n" "file-name:" "$fileName"

printf "  %-34s %s\n" "***** OPTIONAL ARGUMENTS *****"
printf "  %-34s %s\n" "*** General ***"
printf "  %-34s %s\n" "output-dir:" "$outDir"
printf "  %-34s %s\n" "tmp-dir:" "$TMPDIR"
printf "  %-34s %s\n" "MEI-annotation:" "$annot"
printf "  %-34s %s\n\n" "cleanup:" "$cleanup"

printf "  %-34s %s\n" "*** Filters ***"
printf "  %-34s %s\n" "filters:" "$filterList"
printf "  %-34s %s\n" "score-L1-TD0:" "$scoreL1_TD0"
printf "  %-34s %s\n" "score-L1-TD1:" "$scoreL1_TD1"
printf "  %-34s %s\n" "score-L1-TD2:" "$scoreL1_TD2"
printf "  %-34s %s\n" "score-Alu:" "$scoreAlu"
printf "  %-34s %s\n" "score-SVA:" "$scoreSVA"
printf "  %-34s %s\n" "score-ERVK:" "$scoreERVK"
printf "  %-34s %s\n" "score-PSD:" "$scorePSD"
printf "  %-34s %s\n\n" "score-GR:" "$scoreGR"

printf "  %-34s %s\n" "*** Files ***"
printf "  %-34s %s\n\n" "germline-VCF:" "$germlineVCF"


##########
## START #
##########
header="Executing TEIBA $version for $sampleId sample"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
start=$(date +%s)

########################
# I) PRELIMINARY STEPS #
########################

# Create log directory
if [[ ! -d $logsDir ]]; then mkdir $logsDir; fi


# 1) Check the number of insertions in the input file,
#######################################################
# stop execution producing an empty VCF if 0 insertions:
##########################################################

nbLines=`cat $insertions | wc -l`
nbInsertions=$(( nbLines - 1 )) # Substract one in order not to count the header.

printHeader "$sampleId sample has $nbInsertions candidate retrotransposition events"

if [[ $nbInsertions == 0 ]]
then

    log "Donor with 0 retrotransposition events. Generate empty VCF and stop execution\n" "INFO"

    # Print empty VCF only with header
    run "python $EMPTYVCF $fileName -o $outDir 1> $logsDir/1_emptyVCF.out 2> $logsDir/1_emptyVCF.err" "$ECHO"

    ## End
    end=$(date +%s)
    printHeader "TEIBA for $sampleId completed in $(echo "($end-$start)/60" | bc -l | xargs printf "%.2f\n") min "

    exit 0

fi

# 2) For each orphan transduction (TD2) event, extract source region +/- $windowSize to use as target in BLAT alignment.
########################################################################################################################
# Each fasta will be named according to this convention:
########################################################
# - ${srcRegDir}/${family}:${type}:${chr}_${beg}_${end}.src.fa
# where:
#    - MEI: unique identifier
#    - family: MEI family (L1, ALU, SVA, ERVK) and "NA" for pseudogenes
#    - type: TD0 (solo-insertion), TD1 (partnered-transduccion), TD2 (orphan-transduction) and PSD (processed-pseudogenes)
#    - chr: insertion chromosome
#    - beg: insertion beginning
#    - end: insertion end

## parameters
windowSize=1000

## Make source region fasta directory:
if [[ ! -d $srcRegDir ]]; then mkdir $srcRegDir; fi

## Execute the step
step="SOURCES2FASTA"
startTime=$(date +%s)
printHeader "Prepare fasta for assembly"
log "Extracting region downstream of source element for orphan transductions" $step

# extract source gene region
# NOTE: insertion IDs used here correspond to format generated by $CLUSTERS2FASTA
cat $insertions | while read chrP begP endP nReadsP famP readsP chrM begM endM nReadsM famM readsM tdType cytobandId sourceType chrSrc begSrc endSrc strSrc begTd endTd lenRna lenTd psdGene chrExA begExA endExA chrExB begExB endExB grType; do
    if [[ $tdType == "TD2" ]]; then
        targetRegionFile=${famP}":"${tdType}":"${chrP}"_"${endP}"_"${begM}".src.fa"
        extractRegion $genome $chrSrc $begTd $endTd $windowSize $tdType $srcRegDir/$targetRegionFile
    elif [[ $tdType == "PSD" ]]; then
        # sanity check: are both source exons on the same chromosome?
        if [[ "$chrExA" -ne "$chrExB" ]]; then
            echo "[ERROR] PSD event has conflicting source chromsomes. (check input file)" >&2
            continue
        fi
        [[ $begExA -lt $begExB ]] && begEx=$begExA || begEx=$begExB
        [[ $endExA -gt $endExB ]] && endEx=$endExA || endEx=$endExB
        targetRegionFile=${famP}":"${tdType}":"${chrP}"_"${endP}"_"${begM}".src.fa"
        extractRegion $genome $chrExA $begEx $endEx $windowSize $srcRegDir/$targetRegionFile
    fi
done

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


##################
# II) CORE STEPS #
##################

# 1) Cluster clipped reads supporting the same breakpoint. 
###########################################################
# Generate consensus sequence for each cluster
##############################################
# - ${fastaDir}/${family}:${type}:${chr}_${beg}_${end}.fa
# where:
#    - family: MEI family (L1, ALU, SVA, ERVK) and "NA" for pseudogenes
#    - type: TD0 (solo-insertion), TD1 (partnered-transduccion), TD2 (orphan-transduction) and PSD (processed-pseudogenes)
#    - chr: insertion chromosome
#    - beg: insertion beginning
#    - end: insertion end

## Make fasta for assembly directory:
if [[ ! -d $clippedDir ]]; then mkdir $clippedDir; fi

## Execute the step
step="CLUSTER-CLIPPED"
startTime=$(date +%s)
printHeader "Cluster clipped"
log "Cluster clipped" $step
run "python $CLIPPED $insertions $bam -o $outDir --outDir $clippedDir 1> $logsDir/1_clipped.out 2> $logsDir/1_clipped.err" "$ECHO"
endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

# 2) Align the consensus sequence for each breakpoint cluster into the reference genome with blat
##################################################################################################
# It will produce a psl with the blat alignments for the assembled contigs
###########################################################################
# for each cluster insertion
############################
# Psl naming convention:
########################
# - ${blatPath}/${family}:${type}:${chr}_${beg}_${end}.psl
# where:
#    - family: MEI family (L1, ALU, SVA, ERVK) and "NA" for pseudogenes
#    - type: TD0 (solo-insertion), TD1 (partnered-transduccion), TD2 (orphan-transduction) and PSD (processed-pseudogenes)
#    - chr: insertion chromosome
#    - beg: insertion beginning
#    - end: insertion end

## Make blat directory:
if [[ ! -d $blatDir ]]; then mkdir $blatDir; fi

## Execute the step
step="BLAT"
startTime=$(date +%s)

printHeader "Aligning the assembled bkp contigs into the reference genome with blat"

ls $clippedDir | grep '.*fa' | grep -v 'allReadPairs' | while read bkpContig;
do
    bkpContigPath=${clippedDir}/${bkpContig}
    bkpId=${bkpContig%.fa}
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
    if [[ $tdType == "TD2" ]] || [[ $tdType == "PSD" ]]; then
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

## Remove temporary src region directory
if [[ "$cleanup" == "TRUE" ]]; then rm -r $srcRegDir; fi

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

# 3) MEI breakpoint analysis of aligned consensus sequence
###########################################################
# It will produce a VCF with the MEI and many pieces of 
########################################################
# associated information (orientation, TSD, length...)
# ####################################################

## Make bkp analysis directory:
if [[ ! -d $bkpAnalysisDir ]]; then mkdir $bkpAnalysisDir; fi

### Execute the step
## 3.1 Make a list with the MEI ids:
# Output:
# - $outDir/insertions_list.txt
insertionList=$bkpAnalysisDir/insertionList.txt

ls $clippedDir | grep '.*fa' | grep -v 'allReadPairs' | awk '{split($1,a,".fa"); print a[1];}' > $insertionList

## 3.2 For each insertion add the list of read pairs supporting the clusters, the source element, transduction and rearrangement info 
# Output:
# - $bkpAnalysisDir/insertionList_plusInfo.txt
insertionListInfo=$bkpAnalysisDir/insertionList_plusInfo.txt

awk -v OFS='\t' -v fileRef=$insertions -f $ADD_INFO $insertionList > $insertionListInfo

# Remove intermediate files:
if [[ "$cleanup" == "TRUE" ]]; then rm $insertionList ; fi

## 3.3 Prepare input file for insertion breakpoint analysis
# Output:
# - $outDir/paths2bkpAnalysis.txt
paths2bkpAnalysis=$outDir/paths2bkpAnalysis.txt
echo -n "" > $paths2bkpAnalysis

cat $insertionListInfo | while read insertionId readPairsPlus readPairsMinus sourceElementInfo transductionInfo pseudogeneInfo rgInfo;
do
    ## Contigs
    contigPath=${clippedDir}/${insertionId}.fa
    if [[ ! -f $contigPath ]]; then contigPath="NA"; fi
        
    ## Blat
    blatPath=${blatDir}/${insertionId}.psl
    if [[ ! -f $blatPath ]]; then blatPath="NA"; fi    

    printf ${insertionId}"\t"${contigPath}"\t"${blatPath}"\t"${readPairsPlus}"\t"${readPairsMinus}"\t"${sourceElementInfo}"\t"${transductionInfo}"\t"${pseudogeneInfo}"\t"${rgInfo}"\n" >> $paths2bkpAnalysis
done

# Remove intermediate files:
if [[ "$cleanup" == "TRUE" ]]; then rm $insertionListInfo ; fi

## 3.4 Perform breakpoint analysis
# Output:
# - $bkpAnalysisDir/$fileName.vcf
rawVCF=$bkpAnalysisDir/$fileName.vcf

printHeader "Performing MEI breakpoint analysis"

if [ ! -s $rawVCF ];
then
    step="BKP-ANALYSIS"
    startTime=$(date +%s)

    log "Identifying insertion breakpoints, TSD, MEI length, orientation and structure" $step
    run "python $BKP_ANALYSIS $paths2bkpAnalysis $sampleId $fileName $genome --outDir $bkpAnalysisDir 1>> $logsDir/4_bkpAnalysis.out 2>> $logsDir/4_bkpAnalysis.err" "$ECHO"

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
if [[ "$cleanup" == "TRUE" ]]; then rm $paths2bkpAnalysis ; fi
if [[ "$cleanup" == "TRUE" ]]; then rm -r $clippedDir $blatDir ; fi

# 4) Annotate MEI
###################
## Output:
# -  $annotDir/$fileName.annotated.vcf
annotVCF=$annotDir/$fileName.annotated.vcf

## Make MEI annotation directory:
if [[ ! -d $annotDir ]]; then mkdir $annotDir; fi

## If annot disabled copy the previously generated VCF as annot output VCF. 
# Then, annotation step will be skipped as file already exits 
if [[ "$annot" == "FALSE" ]]; then cp $rawVCF $annotVCF; fi

printHeader "Performing MEI breakpoint annotation"

### Execute the step
if [ ! -s $annotVCF ];
then
    step="ANNOTATION"
    startTime=$(date +%s)

    log "Annotating MEI" $step
    run "bash $ANNOTATOR $rawVCF $repeatsDb $driverDb $germlineMEIdb $fileName $annotDir 1>> $logsDir/5_annotation.out 2>> $logsDir/5_annotation.err" "$ECHO"

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
if [[ "$cleanup" == "TRUE" ]]; then rm -r $bkpAnalysisDir ; fi


# 5) Filter MEI
#################
## Output:
# -  $filterDir/$fileName.filtered.vcf
filteredVCF=$filterDir/$fileName.filtered.vcf

## Make MEI filtering directory:
if [[ ! -d $filterDir ]]; then mkdir $filterDir; fi

### Execute the step
# NOTE: for germline variants use a minimum score of 5...

printHeader "Performing MEI filtering"

if [ ! -s $filteredVCF ];
then
    step="FILTER"
    startTime=$(date +%s)
    log "Filtering MEI" $step

    if [[ "$germlineVCF" == "NOT_PROVIDED" ]]
    then
        command="$FILTER $annotVCF $fileName $filterList --score-L1-TD0 $scoreL1_TD0 --score-L1-TD1 $scoreL1_TD1 --score-L1-TD2 $scoreL1_TD2 --score-Alu $scoreAlu --score-SVA $scoreSVA --score-ERVK $scoreERVK --score-PSD $scorePSD --outDir $filterDir 1>> $logsDir/6_filter.out 2>> $logsDir/6_filter.err"
    else
        command="$FILTER $annotVCF $fileName $filterList --score-L1-TD0 $scoreL1_TD0 --score-L1-TD1 $scoreL1_TD1 --score-L1-TD2 $scoreL1_TD2 --score-Alu $scoreAlu --score-SVA $scoreSVA --score-ERVK $scoreERVK --score-PSD $scorePSD --germline-VCF $germlineVCF --outDir $filterDir 1>> $logsDir/6_filter.out 2>> $logsDir/6_filter.err"
    fi

    run "$command" "$ECHO"

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
if [[ "$cleanup" == "TRUE" ]]; then rm -r $annotDir ; fi


#############################
# 6) MAKE OUTPUT VCF AND END #
#############################

## Produce output VCF
finalVCF=$outDir/$fileName.vcf

cp $filteredVCF $finalVCF

## Remove temporary filter directory
if [[ "$cleanup" == "TRUE" ]]; then rm -r $filterDir ; fi

## End
end=$(date +%s)
printHeader "TEIBA for $sampleId sample completed in $(echo "($end-$start)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0
