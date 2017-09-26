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
import scipy.stats as stats
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputVCF', help='Multisample VCF with genotyped source elemetns')
parser.add_argument('metadata', help='Text file with the project and ancestry code per donor Id')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
metadata = args.metadata
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "metadata: ", metadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print



## Start ## 

#### 1. Read metadata file
##########################
# Initialize a dictionary with the following structure:
# - dict: key(ancestryCode) -> donorIdList: [donor1, donor2, donor3]

header("1. Read metadata file")
metadataFile = open(metadata, 'r')

targetAncestriesList = ["AFR", "ASN", "EUR"]
#targetAncestriesList = ["EUR", "AFR"]

# Initialize
ancestryDonorIdListDict = {}

for line in metadataFile:

    line = line.rstrip('\r\n')

    # Skip header
    if not line.startswith("#"):
        line = line.split('\t')

        donorId = line[0]
        status = line[3]
        ancestryCode = line[4]
         
        # Select TraFiC whitelisted donors belonging to AFR, ASN or EUR ancestries
        if (status == "Whitelist") and (ancestryCode in targetAncestriesList):
             
            if ancestryCode not in ancestryDonorIdListDict:

                # Create list
                ancestryDonorIdListDict[ancestryCode] = []

            # Add donor to the list 
            ancestryDonorIdListDict[ancestryCode].append(donorId)
            
### Generate a text file per ancestry containing the donor ids:
targetDonorsList = []

for ancestryCode, donorIdList in ancestryDonorIdListDict.iteritems():
    
    targetDonorsList = targetDonorsList + donorIdList

    ## Open output file
    outFilePath = outDir + '/' + ancestryCode + '_donorIdList.tsv'
    outFile = open(outFilePath, 'w')

    ## Write each donorId in the output file. One id per row
    for donorId in donorIdList:
        row = donorId + '\n'
        outFile.write(row)

#### 2. Read input multi-sample VCF and generate a VCF object
###############################################################
header("2. Read input multi-sample VCF and generate a VCF object")

VCFObj = formats.VCF()
VCFObj4Fst = formats.VCF()

VCFObj.read_VCF_multiSample(inputVCF)


#### 3. Select target donors and source elements
##################################################
header("3. Select target donors and source elements")

## target source elements are rare elements with a MAF < 1%
targetSourceList = ["1p35.2", "1q23.3", "2q21.3", "3p24.1", "3q26.1", "5q13.1", "7p12.3", "7q31.2", "8p23.1f", "9q22.33", "10q25.1", "11p11.2", "11q14.2", "21q21.1"]

## For each MEI:
for MEIObj in VCFObj.lineList:
    
    # Select rare source elements absent in the reference genome:
    #if (MEIObj.infoDict["SRCID"] in targetSourceList) and (MEIObj.alt == "<MEI>"):
    if (MEIObj.infoDict["SRCID"] in targetSourceList):
        print "MEIObj: ", MEIObj, MEIObj.pos, MEIObj.infoDict["SRCID"]

        MEIObj4Fst = MEIObj

        # Select target donors (whitelisted donors from AFR, ASN or EUR ancestries)
        MEIObj4Fst.genotypesDict = { targetDonorId: MEIObj.genotypesDict[targetDonorId] for targetDonorId in targetDonorsList }
        
        # Add MEI object to the VCF object:
        VCFObj4Fst.addLine(MEIObj4Fst)


#### 4. Make output multisample VCF file
##########################################
header("4. Make output multisample VCF file")

fileName = 'rare_germline_source_elements.genotyped.4fst.vcf'
outFilePath = outDir + '/' + fileName

# 5.1 Write header
VCFObj4Fst.write_header(outFilePath)

# 5.2 Write variants
VCFObj4Fst.write_variants_multiSample(targetDonorsList, outFilePath)

####
header("Finished")
