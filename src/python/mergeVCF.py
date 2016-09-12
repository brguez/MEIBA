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
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCFPaths', help='text file with the path to the VCF files will be merged')
parser.add_argument('sampleId', help='Identifier to name output file.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCFPaths = args.VCFPaths
sampleId =  args.sampleId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCFPaths
print "sampleId: ", sampleId
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 

## Start ## 

#### 1. Create VCF object and read input VCF
header("1. Process input VCFs")
paths = open(VCFPaths, 'r')

completeVCFObj = formats.VCF()

#VCFObj.read_VCF(VCF)

for VCFfile in paths:
			
	VCFfile = VCFfile.rstrip('\n')
	VCFObj = formats.VCF()
	donorIdList = VCFObj.read_VCF_multiSample(VCFfile)

	for lineObj in VCFObj.lineList:
		completeVCFObj.addLine(lineObj)

	if completeVCFObj.header == "":
		completeVCFObj.header = VCFObj.header				

completeVCFObj.lineList = completeVCFObj.sort() 

#### . 
header("2. Split into N VCFs with N the number of chromosomes")

outFilePath = outDir + '/' + sampleId + ".vcf"

# .1 Write header
completeVCFObj.write_header(outFilePath)

# .2 Write variants
completeVCFObj.write_variants_multiSample(donorIdList, outFilePath)


header("Finished")
