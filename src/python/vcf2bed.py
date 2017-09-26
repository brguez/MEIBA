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
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import seaborn as sns
import scipy
import formats

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('vcf', help='Multisample VCF containing genotyped MEI')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.vcf
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Read input VCF and generate a VCF object
#################################################
header("1. Process VCF as input")

VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

#### 2. Write each MEI as a BED data line
##########################################
header("2. Write each MEI as a BED data line")

# Open output file
outFilePath = outDir + '/' + fileName + '.bed'
outFile = open(outFilePath, 'a')

# Write header:
row = '#chrom' + "\t" + 'beg' +  "\t" + 'end' + "\t" + 'iClass' + "\t" + 'iType' + "\t" + 'srcId' + "\n"
outFile.write(row)

## For each MEI
for MEIObj in VCFObj.lineList:

    ## Select only those MEI that passes all the filters
    if (MEIObj.filter == "PASS"):

        chrom = MEIObj.chrom
        beg = str(MEIObj.pos - 1000)
        end = str(MEIObj.pos + 1000)
        iClass = MEIObj.infoDict["CLASS"]
        iType = MEIObj.infoDict["TYPE"]
        srcId = MEIObj.infoDict['SRCID'] if 'SRCID' in MEIObj.infoDict else 'NA'
    
        row = chrom + '\t' + beg + '\t' + end + '\t' + iClass + '\t' + iType + '\t' + srcId + '\n'
        outFile.write(row)
       

#### End
header("FINISH!!")


