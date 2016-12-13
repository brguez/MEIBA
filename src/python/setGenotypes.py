#!/usr/bin/env python
#coding: utf-8

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


def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string



    ## Iterate over donor id/BAM files, genotyping one donor at a time
    for line in donorIdBamPathList:

        # Extract donorId and BAM file path
        donorId = line[0]
        BAMPath = line[1]

        ## Genotype each MEI polymorphism for the current donor
        for VCFlineObj in VCFObj.lineList:

            genotypeField = genotype + ':' + str(nbReadsALT) + ':' + str(totalNbReads)
            VCFlineObj.genotypesDict[donorId] = genotypeField



#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time

# Global variables:
global debugBool ## debug logging mode. Boolean.

## Get user's input ##
parser = argparse.ArgumentParser(description= "Set a genotype for all the germline mobile element insertions (MEI) across a group of donors. Produce a multi-sample VCF as output")
parser.add_argument('VCF', help='VCF with the MEI')
parser.add_argument('donorIds', help='Tab separated text file with one column: 1) donor ids. One donor per row.')
parser.add_argument('outFileName', help='Identifier to name the output file.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
donorIds = args.donorIds
outFileName =  args.outFileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "donorIds: ", donorIds
print "outFileName: ", outFileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Create VCF object and read input VCF

header("1. Process input VCF with polymorphic MEI for genotyping")

VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)


#### 2. Add donors genotype
header("2. For each MEI set donors genotype")

# For each MEI 
for VCFlineObj in VCFObj.lineList:

    # For each donor set MEI genotype
    for line in open(donorIds):
        donorId = line.rstrip('\n')
        genotype = "0/0:.:."
    
        print "donorId: ", donorId, genotype
        VCFlineObj.genotypesDict[donorId] = genotype
          
#### 3. Make output multisample VCF file
header("3. Produce multi-sample VCF file as ouput")

fileName = outFileName + '.vcf'
outFilePath = outDir + '/' + fileName

# 3.1 Write header
VCFObj.write_header(outFilePath)

# 3.2 Write variants
donorIdList = [line.rstrip('\n') for line in open(donorIds)]

print donorIdList
VCFObj.write_variants_multiSample(donorIdList, outFilePath)

header("Finished")        
