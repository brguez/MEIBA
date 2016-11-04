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


#### FUNCTIONS ####
def splitVCF(VCFObj, outDir):
    '''
    '''

    ## 1. Make list of chromosome ids with the same order as in the original VCF
    chromList = [VCFlineObj.chrom for VCFlineObj in VCFObj.lineList]
    chromList = sorted(set(chromList), key=lambda x: chromList.index(x))

    ## 2. Generate a dictionary with the following format:
    # key(chromId) -> value(list of variant objects in a given chromosome)
    VCFlinesDict = {}

    for lineObj in VCFObj.lineList:

        # Iniciate key
        if not lineObj.chrom in VCFlinesDict:
            VCFlinesDict[lineObj.chrom] = []

        VCFlinesDict[lineObj.chrom].append(lineObj)

    ## 3. Generate a dictionary with the following format:
    # key(chromId) -> value(VCF object containing variants for a given chromosome)
    VCFDict = {}

    for chrom in chromList:
        VCFChrObj = formats.VCF()
        VCFChrObj.lineList =  VCFlinesDict[chrom]
        VCFDict[chrom] = VCFChrObj

    ## 4. Write a VCF per chromosome containing those variants in the chrom.

    # Iterate over the chromosomes generating a VCF at each time
    for chrom in chromList:
        outFilePath = outDir + '/' + chrom + '.vcf'

        # Write meta-information
        VCFObj.write_header(outFilePath)

        # Write header
        row = '#CHROM' + "\t" + 'POS' + "\t" + 'ID' + "\t" + 'REF' + "\t" + 'ALT' + "\t" + 'QUAL' + "\t" + 'FILTER' + "\t" + 'INFO' + "\t" + 'FORMAT' + "\t" + 'CONSENSUS' + "\n"
        outFile = open(outFilePath, 'a')
        outFile.write(row)
        outFile.close()

        # Write variants
        VCFDict[chrom].write_variants(outFilePath)

    return  VCFDict


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='input VCF')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 


#### 1. Create VCF object and read input VCF

header("1. Process input VCF")

VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)

#### 2.
header("2. Split into N VCFs with N the number of chromosomes")
splitVCF(VCFObj, outDir)


header("Finished")
