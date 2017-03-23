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
def splitVCF(inputVCFObj, nbChunks, outDir):
    '''
    '''

    ## 1. Split variants into X evenly sized chunks.
    print "* nb.insertions: ", len(inputVCFObj.lineList)

    chunkSize = int(round(len(inputVCFObj.lineList) / float(nbChunks)))

    chunkList = [inputVCFObj.lineList[i:i + chunkSize] for i in xrange(0, len(inputVCFObj.lineList), chunkSize)]

    print "* chunkSize: ", chunkSize
    print "* Number_chunks: ", len(chunkList)

    chunkId = 1

    print "* Create VCFs: ", len(chunkList)

    ## 2. For each chunk generate a VCF file:
    for chunk in chunkList:
        
        print "     VCF_id: ", chunkId

        ## Create VCF object only containing those variant in the chunk
        VCFObj = formats.VCF()
        VCFObj.header = inputVCFObj.header
        VCFObj.lineList = chunk
        
        ## Write VCF file
        outFilePath = outDir + '/chunk_' + str(chunkId) + '.vcf'
        
        print "     VCF_path: ", outFilePath

        # Write meta-information
        VCFObj.write_header(outFilePath)

        # Write variants
        VCFObj.write_variants(outFilePath)

        ## Update chunk identifier
        chunkId  += 1        




#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time

## Get user's input ##
parser = argparse.ArgumentParser(description= "Split variants in the VCF into X evenly sized chunks. X is an optional parameter")
parser.add_argument('VCF', help='input VCF')
parser.add_argument('--nb-clunks', default=1, dest='nbChunks', type=int, help='Number of chunks. This is an approximation. The final number can be deviated. Default: 1' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='Output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
nbChunks = args.nbChunks
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "nb-clunks: ", nbChunks
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Create VCF object and read input VCF

header("1. Process input VCF")

VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)

#### 2. Split into N VCFs with N the number of chromosomes
header("2. Split input VCF into N VCFs ")
splitVCF(VCFObj, nbChunks, outDir)


header("Finished")
