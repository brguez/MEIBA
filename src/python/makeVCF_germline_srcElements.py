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
##INFO=<ID=CLASS,Number=1,Type=String,Description="L1 element class (L1, L1HS, L1PA2 or L1PA3)">
##INFO=<ID=BKPB,Number=1,Type=Integer,Description="MEI right-most breakpoint position (bkp B). Left-most breakpoint position (bkp A) represented in the POS field">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion DNA strand (+ or -)">
##INFO=<ID=SRCID,Number=1,Type=String,Description="Source element cytoband identifier. Only for gemline source elements">
##INFO=<ID=POLYMORPHIC,Number=0,Type=Flag,Description="Polymorphic source element.">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Novel source element. Not already reported in Tubio et. al. Science (2014)">
##FORMAT=<ID=RCP,Number=1,Type=Integer,Description="Count of positive cluster supporting reads">
##FORMAT=<ID=RCN,Number=1,Type=Integer,Description="Count of negative cluster supporting reads">
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
parser.add_argument('metadata', help='Source elements metadata')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
metadata = args.metadata
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "metadata: ", metadata
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read metadata file and create a VCF object with source element information
#####################################################################################
header("1. Read metadata file and create a VCF object with source element information")

# Initialize VCF object
VCFObj = formats.VCF()

metadata = open(metadata, 'r')

# For each source element in the bed
for line in metadata:


    # Discard header
    if not line.startswith("#"):
                 
        line = line.rstrip('\n')
        line = line.rstrip('\r')
        line = line.split('\t')

        ### 1.  Get input values  
        #cytobandId	chrom	bkpA	    bkpB	    class	orientation	ref	novel	polymorphic
     
        cytobandId = line[0]
        chrom = 	line[1]
        bkpA = int(line[2]) 
        bkpB = line[3]  
        category = line[4]	
        orientation = line[5]
        refBool = int(line[6])	       
        novelBool = int(line[7])
        polymorphicBool = int(line[8])

        ### 2. Define vcf fields
        CHROM = chrom 	
        POS	= bkpA
        ID = '.' 
        
        ## Reference and alternative field
        # a) MEI in the reference genome
        if (refBool == 1):
            REF	= '<MEI>'
            ALT = '.'
    
        # b) MEI not in the reference
        else:
            REF	= '.'
            ALT = '<MEI>'

        ## Orientation
        # a) Plus
        if (orientation == 'plus'):
            orientation = '+'

        # b) Minus
        else:
            orientation = '-'
        
        ## Polymorphic flag
        # a) Polymorphic element
        if (polymorphicBool == 1):
            polymorphicFlag = ';POLYMORPHIC'

        # b) Fixed element
        else:
            polymorphicFlag = ''

        ## Novel flag
        # a) Novel element
        if (novelBool == 1):
            novelFlag = ';NOVEL'

        # b) Known element
        else:
            novelFlag = ''

        QUAL = '.'    
        FILTER = '.' 
        INFO	 = 'SVTYPE=<MEI>;' + 'CLASS=' + category + ';BKPB=' + str(bkpB) + ';STRAND=' + orientation + ';SRCID=' + cytobandId + polymorphicFlag + novelFlag
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

# 3.2 Sort variants
VCFObj.lineList = VCFObj.sort()

# 3.3 Write variants
VCFObj.write_variants(outFilePath)

#### END
header("Finished")
