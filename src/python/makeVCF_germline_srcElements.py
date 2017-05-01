#!/usr/bin/env python
#coding: utf-8


#### FUNCTIONS ####
def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"


def subHeader(string):
    """
        Display  subheader
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, "**", string, "**"


def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


def write_header(outFilePath):
    ## Define variables
    date = time.strftime("%Y%m%d")

    context = {
    "date": date,
    "source": "TraFiCv2.0",
    "reference": "hs37d5",
     }

    ## Header template
    template = """##fileformat=VCFv4.2
##fileDate={date}
##source={source}
##reference={reference}
##contig=<ID=1,assembly=GRCh37,length=249250621,species=human>
##contig=<ID=2,assembly=GRCh37,length=243199373,species=human>
##contig=<ID=3,assembly=GRCh37,length=198022430,species=human>
##contig=<ID=4,assembly=GRCh37,length=191154276,species=human>
##contig=<ID=5,assembly=GRCh37,length=180915260,species=human>
##contig=<ID=6,assembly=GRCh37,length=171115067,species=human>
##contig=<ID=7,assembly=GRCh37,length=159138663,species=human>
##contig=<ID=8,assembly=GRCh37,length=146364022,species=human>
##contig=<ID=9,assembly=GRCh37,length=141213431,species=human>
##contig=<ID=10,assembly=GRCh37,length=135534747,species=human>
##contig=<ID=11,assembly=GRCh37,length=135006516,species=human>
##contig=<ID=12,assembly=GRCh37,length=133851895,species=human>
##contig=<ID=13,assembly=GRCh37,length=115169878,species=human>
##contig=<ID=14,assembly=GRCh37,length=107349540,species=human>
##contig=<ID=15,assembly=GRCh37,length=102531392,species=human>
##contig=<ID=16,assembly=GRCh37,length=90354753,species=human>
##contig=<ID=17,assembly=GRCh37,length=81195210,species=human>
##contig=<ID=18,assembly=GRCh37,length=78077248,species=human>
##contig=<ID=19,assembly=GRCh37,length=59128983,species=human>
##contig=<ID=20,assembly=GRCh37,length=63025520,species=human>
##contig=<ID=21,assembly=GRCh37,length=48129895,species=human>
##contig=<ID=22,assembly=GRCh37,length=51304566,species=human>
##contig=<ID=hs37d5,assembly=GRCh37,length=35477943,species=human>
##contig=<ID=GL000191.1,assembly=GRCh37,length=106433,species=human>
##contig=<ID=GL000192.1,assembly=GRCh37,length=547496,species=human>
##contig=<ID=GL000193.1,assembly=GRCh37,length=189789,species=human>
##contig=<ID=GL000194.1,assembly=GRCh37,length=191469,species=human>
##contig=<ID=GL000195.1,assembly=GRCh37,length=182896,species=human>
##contig=<ID=GL000196.1,assembly=GRCh37,length=38914,species=human>
##contig=<ID=GL000197.1,assembly=GRCh37,length=37175,species=human>
##contig=<ID=GL000198.1,assembly=GRCh37,length=90085,species=human>
##contig=<ID=GL000199.1,assembly=GRCh37,length=169874,species=human>
##contig=<ID=GL000200.1,assembly=GRCh37,length=187035,species=human>
##contig=<ID=GL000201.1,assembly=GRCh37,length=36148,species=human>
##contig=<ID=GL000202.1,assembly=GRCh37,length=40103,species=human>
##contig=<ID=GL000203.1,assembly=GRCh37,length=37498,species=human>
##contig=<ID=GL000204.1,assembly=GRCh37,length=81310,species=human>
##contig=<ID=GL000205.1,assembly=GRCh37,length=174588,species=human>
##contig=<ID=GL000206.1,assembly=GRCh37,length=41001,species=human>
##contig=<ID=GL000207.1,assembly=GRCh37,length=4262,species=human>
##contig=<ID=GL000208.1,assembly=GRCh37,length=92689,species=human>
##contig=<ID=GL000209.1,assembly=GRCh37,length=159169,species=human>
##contig=<ID=GL000210.1,assembly=GRCh37,length=27682,species=human>
##contig=<ID=GL000211.1,assembly=GRCh37,length=166566,species=human>
##contig=<ID=GL000212.1,assembly=GRCh37,length=186858,species=human>
##contig=<ID=GL000213.1,assembly=GRCh37,length=164239,species=human>
##contig=<ID=GL000214.1,assembly=GRCh37,length=137718,species=human>
##contig=<ID=GL000215.1,assembly=GRCh37,length=172545,species=human>
##contig=<ID=GL000216.1,assembly=GRCh37,length=172294,species=human>
##contig=<ID=GL000217.1,assembly=GRCh37,length=172149,species=human>
##contig=<ID=GL000218.1,assembly=GRCh37,length=161147,species=human>
##contig=<ID=GL000219.1,assembly=GRCh37,length=179198,species=human>
##contig=<ID=GL000220.1,assembly=GRCh37,length=161802,species=human>
##contig=<ID=GL000221.1,assembly=GRCh37,length=155397,species=human>
##contig=<ID=GL000222.1,assembly=GRCh37,length=186861,species=human>
##contig=<ID=GL000223.1,assembly=GRCh37,length=180455,species=human>
##contig=<ID=GL000224.1,assembly=GRCh37,length=179693,species=human>
##contig=<ID=GL000225.1,assembly=GRCh37,length=211173,species=human>
##contig=<ID=GL000226.1,assembly=GRCh37,length=15008,species=human>
##contig=<ID=GL000227.1,assembly=GRCh37,length=128374,species=human>
##contig=<ID=GL000228.1,assembly=GRCh37,length=129120,species=human>
##contig=<ID=GL000229.1,assembly=GRCh37,length=19913,species=human>
##contig=<ID=GL000230.1,assembly=GRCh37,length=43691,species=human>
##contig=<ID=GL000231.1,assembly=GRCh37,length=27386,species=human>
##contig=<ID=GL000232.1,assembly=GRCh37,length=40652,species=human>
##contig=<ID=GL000233.1,assembly=GRCh37,length=45941,species=human>
##contig=<ID=GL000234.1,assembly=GRCh37,length=40531,species=human>
##contig=<ID=GL000235.1,assembly=GRCh37,length=34474,species=human>
##contig=<ID=GL000236.1,assembly=GRCh37,length=41934,species=human>
##contig=<ID=GL000237.1,assembly=GRCh37,length=45867,species=human>
##contig=<ID=GL000238.1,assembly=GRCh37,length=39939,species=human>
##contig=<ID=GL000239.1,assembly=GRCh37,length=33824,species=human>
##contig=<ID=GL000240.1,assembly=GRCh37,length=41933,species=human>
##contig=<ID=GL000241.1,assembly=GRCh37,length=42152,species=human>
##contig=<ID=GL000242.1,assembly=GRCh37,length=43523,species=human>
##contig=<ID=GL000243.1,assembly=GRCh37,length=43341,species=human>
##contig=<ID=GL000244.1,assembly=GRCh37,length=39929,species=human>
##contig=<ID=GL000245.1,assembly=GRCh37,length=36651,species=human>
##contig=<ID=GL000246.1,assembly=GRCh37,length=38154,species=human>
##contig=<ID=GL000247.1,assembly=GRCh37,length=36422,species=human>
##contig=<ID=GL000248.1,assembly=GRCh37,length=39786,species=human>
##contig=<ID=GL000249.1,assembly=GRCh37,length=38502,species=human>
##contig=<ID=MT,assembly=GRCh37,length=16569,species=human>
##contig=<ID=NC_007605,assembly=GRCh37,length=171823,species=human>
##contig=<ID=X,assembly=GRCh37,length=155270560,species=human>
##contig=<ID=Y,assembly=GRCh37,length=59373566,species=human>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant. (All sequence is on the plus strand and in the forward direction).">
##INFO=<ID=CLASS,Number=1,Type=String,Description="Mobile element class (L1, ALU, SVA or ERVK)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Insertion type (TD0: solo, TD1: partnered-3'transduction, TD2: orphan-3'transduction, PSD: processed-pseudogene)>
##INFO=<ID=SCORE,Number=1,Type=Integer,Description="Insertion score (5: 5' and 3' breakpoints (bkp) assembled, 4: 3'bkp assembled, 3: 5'bkp assembled, 2: no bkp assembled, 1: inconsistent (contradictory orientation, bkp or TSD))">
##INFO=<ID=MANUAL,Number=0,Type=Flag,Description="MEI manually verified and curated through BAM inspection (Only used for PSD)">
##INFO=<ID=BKPB,Number=1,Type=Integer,Description="MEI right-most breakpoint position (bkp B). Left-most breakpoint position (bkp A) represented in the POS field">
##INFO=<ID=CIPOS,Number=1,Type=Integer,Description="Confidence interval around insertion breakpoints">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion DNA strand (+ or -)">
##INFO=<ID=STRUCT,Number=1,Type=String,Description="Mobile element structure (INV: 5'inverted, DEL: 5'deleted, FULL: full-length)">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Mobile element length">
##INFO=<ID=TSLEN,Number=1,Type=Integer,Description="Target site duplication (+_value) or deletion (-_value) length">
##INFO=<ID=TSSEQ,Number=1,Type=String,Description="Target site duplication sequence">
##INFO=<ID=POLYA,Number=1,Type=String,Description="Poly-A sequence">
##INFO=<ID=SRCID,Number=1,Type=String,Description="Source element cytoband identifier. Only for gemline source elements">
##INFO=<ID=SRCTYPE,Number=1,Type=String,Description="Source element type (GERMLINE or SOMATIC)">
##INFO=<ID=SRC,Number=1,Type=String,Description="Coordinates of the source element that mediated the transduction in the format: chrom_beg_end">
##INFO=<ID=TDC,Number=1,Type=String,Description="Begin and end coordinates of the integrated transduced or pseudogene sequence in the format: chrom_beg_end">
##INFO=<ID=TDLEN,Number=1,Type=Integer,Description="Transduced region length">
##INFO=<ID=TDLENR,Number=1,Type=Integer,Description="Transduced region length at RNA level">
##INFO=<ID=SRCGENE,Number=1,Type=String,Description="Source gene of the processed pseudogene insertion">
##INFO=<ID=GERMDB,Number=1,Type=String,Description="MEI already reported as germinal in a database (1KGENOMES: 1000 genomes project (source_papers_doi: 10.1038/nature15394 and 10.1073/pnas.1602336113), TRAFIC: TraFic in-house database)">
##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region where the mobile element is inserted (exonic, splicing, ncRNA, UTR5, UTR3, intronic, upstream, downstream, intergenic)">
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO gene symbol">
##INFO=<ID=ROLE,Number=1,Type=String,Description="Role in cancer (oncogene, TSG: tumor suppressor gene, oncogene/TSG: both roles)">
##INFO=<ID=COSMIC,Number=0,Type=Flag,Description="Reported as cancer driver in COSMIC cancer gene census database">
##INFO=<ID=CPG,Number=0,Type=Flag,Description="Reported as cancer predisposition gene in 10.1038/nature12981 (DOI).">
##INFO=<ID=REP,Number=1,Type=String,Description="Repetitive element overlapping the insertion breakpoint">
##INFO=<ID=DIV,Number=1,Type=Integer,Description="Millidivergence of the overlapping repetitive element with respect a consensus sequence">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Assembled contig sequence spanning 1st bkp (lowest genomic position)">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Assembled contig sequence spanning 2nd bkp (highest genomic position)">
##INFO=<ID=RP,Number=.,Type=String,Description="Reads from the tumour sample and positive cluster that support this insertion">
##INFO=<ID=RN,Number=.,Type=String,Description="Reads from the tumour sample and negative cluster that support this insertion">
##FILTER=<ID=SCORE,Description="Insertion with an score < threshold">
##FILTER=<ID=REP,Description="Insertion overlapping a satellite region or a repetitive element of the same class">
##FILTER=<ID=DUP,Description="Duplicated MEI call">
##FILTER=<ID=GERMLINE,Description="Germline MEI miscalled as somatic">
##FILTER=<ID=TD,Description="L1 transduction incorrectly identified as a processed pseudogene insertion">
##FORMAT=<ID=RCP,Number=1,Type=Integer,Description="Count of positive cluster supporting reads">
##FORMAT=<ID=RCN,Number=1,Type=Integer,Description="Count of negative cluster supporting reads">
##FORMAT=<ID=SL,Number=1,Type=String,Description="List of samples where the variant was found (specially relevant for multi-tumor donors)">
##FORMAT=<ID=REPR,Number=1,Type=String,Description="Sample selected as representative among all the samples where the variant was found (specially relevant for multi-sample VCF).">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting the variant in this sample">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads covering variant location in this sample">
"""

    ## Replace variables into the template and print header into the output file
    with  open(outFilePath,'w') as outFile:
        outFile.write(template.format(**context))


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time

## Get user's input ##
parser = argparse.ArgumentParser(description= "Make VCF file containing annotated source elements")
parser.add_argument('inputBED', help='BED file with annotated source elements')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputBED = args.inputBED
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputBED: ", inputBED
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Read input BED file and create a VCF object with source element information
#####################################################################################
header("1. Read input BED file and create a VCF object with source element information")

# Initialize VCF object
VCFObj = formats.VCF()

inputBED = open(inputBED, 'r')

# For each source element in the bed
for line in inputBED:

    # Discard header
    if not line.startswith("#"):
                 
        line = line.rstrip('\n')
        line = line.rstrip('\r')
        line = line.split('\t')

        ### 1.  Get input values       
        chrom = line[0]	
        bkpA = int(line[1]) 
        bkpB = line[2]  
        category = line[3]	
        orientation = line[4]
        refBool = line[5]	       
        
        ### 2. Define vcf fields
        CHROM = chrom 	
        POS	= bkpA
        ID = '.' 
        
        ## Reference and alternative field
        # a) MEI in the reference genome
        if (refBool == "1"):
            REF	= '<MEI>'
            ALT = '.'
    
        # b) MEI not in the reference
        else:
            REF	= '.'
            ALT = '<MEI>'

        ## Source element orientation
        # a) Plus orientation
        if (orientation == "plus"):
            orientation = "+"
    
        # b) Minus orientation
        else:
            orientation = "-"

        QUAL = '.'    
        FILTER = '.' 
        INFO	 = 'SVTYPE=<MEI>;' + 'CLASS=' + category + ';TYPE=TD0;BKPB=' + str(bkpB) + ';CIPOS=0;STRAND=' + orientation 
        FORMAT = 'RCP:RCN' 
        GT = '.:.'

        ### 3. Create VCFline object and add it to the VCF
        lineList = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, GT]
        
        VCFlineObj = formats.VCFline(lineList)
        VCFlineObj.info = VCFlineObj.make_info()        
        VCFObj.addLine(VCFlineObj)      
    

## 3. Make output VCF
header("2. Write output VCF")

outFilePath = outDir + '/' + fileName + '.vcf'

# 3.1 Write header
write_header(outFilePath)

# 3.2 Write variants
VCFObj.write_variants(outFilePath)

#### END
header("Finished")
