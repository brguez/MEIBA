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
import formats

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('vcf', help='VCF containing MEI called with TraFiC-mem')
parser.add_argument('sampleId', help='Sample identifier. Output file will be named accordingly')
parser.add_argument('--overhang', default=0, dest='overhang', type=int, help='Extend beg and end interval with +- X base pairs. Default: 0')

parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.vcf
sampleId = args.sampleId
overhang = args.overhang
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "sampleId: ", sampleId
print "overhang: ", overhang
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read input VCF and generate a VCF object
###################################################
header("1. Read input VCF and generate a VCF object")

VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)


#### 2. Select MEI, extract relevant information and produce bed
######################################################################

header("2. Select MEI, extract relevant information and produce bed")

## Write header
outFilePath = outDir + '/MEI_' + sampleId + '.bed'
outFile = open(outFilePath, 'w')


fieldsList = ['#chrom', 'beg', 'end', 'family', 'length']
row = '\t'.join(fieldsList) + '\n' 
outFile.write(row)

## For each SV in the VCF
for VCFline in VCFObj.lineList:

    ## Select MEI that passes the filters
    if (VCFline.filter == 'PASS'):
        
		## Define beginning and end coordinates       
		beg = str(VCFline.pos - overhang)
		end = str(int(VCFline.infoDict['BKPB']) + overhang) if 'BKPB' in VCFline.infoDict else str(VCFline.pos + overhang)
		family = VCFline.infoDict['CLASS']
		length = VCFline.infoDict['LEN'] if 'LEN' in VCFline.infoDict else 'NA'
		
        ## Write into output file
		fieldsList = [VCFline.chrom, beg, end, family, length]

		row = '\t'.join(fieldsList) + '\n' 

		outFile.write(row)

        
#### End
header("FINISH!!")

