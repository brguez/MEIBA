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
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('inputVCF', help='Multi-sample VCF file containing genotyped MEI from 1000 genomes project')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF =  args.inputVCF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"


## Start ## 

#### 1. Read input VCF and generate VCF object
###############################################
header("1. Process input VCF")

## Read  VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)


##### 2. Produce output table with the MAF in Europeans for each variant
##########################################################################
header("2. Produce output table with the MAF in Europeans for each variant")

# Open output file
outFilePath = outDir + '/MEI_MAF_1KGP.EUR.bed'
outFile = open(outFilePath, 'w')

# Write header:
row = '#chrom' + "\t" + 'beg' +  "\t" + 'end' +  "\t" + 'Class' +  "\t" + 'EUR_MAF' + "\n"
outFile.write(row)

# For each MEI in the VCF
for MEIObj in VCFObj.lineList: 

    # a) Alu
    if (MEIObj.infoDict["SVTYPE"] == "ALU"):
        MEIclass = "Alu"

    # b) LINE 1    
    elif (MEIObj.infoDict["SVTYPE"] == "LINE1"):
        MEIclass = "L1"

    # c) SVA
    else:
        MEIclass = "SVA"

    EURMAF = float(MEIObj.infoDict["EUR_AF"]) * 100

    # Write data into the output file    
    row = MEIObj.chrom + "\t" + str(MEIObj.pos) +  "\t" + str(MEIObj.pos) +  "\t" + MEIclass +  "\t" + str(EURMAF) + "\n"

    outFile.write(row)

#### END
header("Finished")


