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

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time

## Get user's input ##
parser = argparse.ArgumentParser(description= "Select common MEI from genotyped 1000 genomes MEI and PCAWG not reduntand MEI VCFs. Produce 2 VCF as output: 1) Common not redundant MEI one-sample VCF and 2) common MEI multi-sample VCF with 1KGP genotyping results and PCAWG MEI features")
parser.add_argument('VCF1KGP', help='Multi-sample VCF file containing genotyped MEI from 1000 genomes project')
parser.add_argument('VCFPCAWG', help='Multi-sample VCF file containing genotyped MEI from PCAWG project')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF1KGP = args.VCF1KGP
VCFPCAWG =  args.VCFPCAWG
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF1KGP: ", VCF1KGP
print "VCFPCAWG: ", VCFPCAWG
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print



## Start ## 

#### 1. Read input VCFs and generate VCF objects
#############################################################
header("1. Process input VCFs ")

## 1000 genomes multi-sample VCF
VCFObj1KGP = formats.VCF()
donorIdList1KGP = VCFObj1KGP.read_VCF_multiSample(VCF1KGP)

## PCAWG VCF
VCFObjPCAWG = formats.VCF()
VCFObjPCAWG.read_VCF(VCFPCAWG)


#### 2. Select common MEI from 1000 genomes and PCAWG
############################################################
header("2. Select common MEI from 1000 genomes and PCAWG")

outVCFObj = formats.VCF()

counter = 1

## For each PCAWG MEI
for MEIObjPCAWG in VCFObjPCAWG.lineList:

    msg = "Assessing if there is 1KGP MEI overlapping PCAWG MEI: " + str(counter)
    info(msg)

    print "PCAWG: ", MEIObjPCAWG.chrom, MEIObjPCAWG.pos, MEIObjPCAWG.infoDict['CLASS']

    # For each 1KGP MEI
    for MEIObj1KGP  in VCFObj1KGP.lineList:

        MEItype = MEIObj1KGP.alt.split(":")[2][:-1]

        # Use our MEI type naming convention (1000 genomes uses a different one..)
        if (MEItype == "LINE1"):
            MEItype = "L1"

        elif (MEItype == "ALU"):
            MEItype = "Alu"

        else:
            MEItype = "SVA"

        # Current 1KGP MEI in the same chromosome and of the same type as PCAWG MEI
        if (str(MEIObj1KGP.chrom) == str(MEIObjPCAWG.chrom)) and (MEItype == MEIObjPCAWG.infoDict['CLASS']):

            # Define range +-20 bp to search for overlap between PCAWG and 1KGP MEI breakpoints
            beg = MEIObjPCAWG.pos - 20
            end = MEIObjPCAWG.pos + 20

            # Current 1KGP MEI overlapping current PCAWG breakpoint range
            if (int(MEIObj1KGP.pos) >= beg) and (int(MEIObj1KGP.pos) <= end):

                print "overlapping_1KGP: ", MEIObj1KGP.chrom, MEIObj1KGP.pos, MEItype

                ## Create MEI object containing PCAWG features (bkp, lenght..) but 1KGP genotyping
                # PCAWG features:
                featureList = [ MEIObjPCAWG.chrom, MEIObjPCAWG.pos, MEIObjPCAWG.id, MEIObjPCAWG.ref, MEIObjPCAWG.alt, MEIObjPCAWG.qual, MEIObjPCAWG.filter, MEIObjPCAWG.info, "GT", MEIObjPCAWG.genotype ]
                MEIObj = formats.VCFline(featureList)

                # 1KGP genotyping
                MEIObj.genotypesDict =  MEIObj1KGP.genotypesDict

                ## Add MEI object to ouput file:
                outVCFObj.addLine(MEIObj)

                # Stop iterating when match found
                break

            # Current MEI with a breakpoint with a position higher than end range -> stop iterating
            if (int(MEIObj1KGP.pos) >= end):
                break
    counter+=1

#### 3. Make multi-sample VCF file with common MEI as ouput
##########################################################
header("3. Make multi-sample VCF file with common MEI as ouput")

fileName = 'PCAWG_1KGP_commonMEI_1KGPgenotyped.vcf'
outFilePath = outDir + '/' + fileName

# 3.1 Write header
outVCFObj.header = VCFObjPCAWG.header
outVCFObj.write_header(outFilePath)

# 3.2 Write variants
outVCFObj.write_variants_multiSample(donorIdList1KGP, outFilePath)


#### 4. Make not redundant VCF file with common MEI as ouput
##########################################################
header("4. Make not redundant VCF file with common MEI as ouput")

fileName = 'PCAWG_1KGP_commonMEI_notRedundant.vcf'
outFilePath = outDir + '/' + fileName

# 3.1 Write header
outVCFObj.write_header(outFilePath)

# 3.2 Write variants
outVCFObj.write_variants(outFilePath)


#### END
header("Finished")
