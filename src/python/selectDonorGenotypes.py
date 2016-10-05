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
parser = argparse.ArgumentParser(description= "")
parser.add_argument('inputVCF', help='Multi-sample VCF file containing genotyped MEI')
parser.add_argument('targetDonors', help='')
parser.add_argument('sampleId', help='Identifier to name output file.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
targetDonors =  args.targetDonors
sampleId =  args.sampleId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "targetDonors: ", targetDonors
print "sampleId: ", sampleId
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 

## Start ## 

#### 1. Read input VCF and generate VCF object
###############################################
header("1. Process input VCFs ")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Read target donors txt and make list with target donor ids
##################################################################
header("2. Read target donors txt")

targetDonorsList = [line.rstrip('\n') for line in open(targetDonors)]

#### 3. Make multi-sample VCF file with selected donors as ouput
################################################################
header("3. Make multi-sample VCF file with selected donors as ouput")

outFilePath = outDir + '/' + sampleId + '.vcf'

# 3.1 Write header
VCFObj.write_header(outFilePath)

# 3.2 Write variants
VCFObj.write_variants_multiSample(targetDonorsList, outFilePath)

#### END
header("Finished")
