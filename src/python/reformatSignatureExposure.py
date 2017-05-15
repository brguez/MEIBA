#!/usr/bin/env python
#coding: utf-8

def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"

def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
from operator import itemgetter, attrgetter, methodcaller

## Get user's input ##
parser = argparse.ArgumentParser(description= "Reformat signature exposure file")
parser.add_argument('signatures', help='csv file with signature exposure per sample id')
parser.add_argument('metadata', help='Donor metadata info')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
signatures = args.signatures
metadata = args.metadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "signatures: ", signatures
print "metadata: ", metadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


# 1. Make dictionary relating sampleId with specimenId
#######################################################

info("1. Make dictionary relating sampleId with specimenId")
metadata = open(metadata, 'r')

medatadaDict = {}

for line in metadata:
    line = line.rstrip('\n')
    line = line.split("\t")
    specimenList = line[17].split(",")
    sampleList = line[18].split(",")

    for index, specimenId in enumerate(specimenList):
        sampleId = sampleList[index]  
        
        medatadaDict[sampleId] = specimenId         


# 2. Write reformatted line into the output file
##################################################
info("2. Write reformatted line into the output file")

outPath = outDir + '/pcawg7_sigs_samp_exposure.tsv'
outFile = open(outPath, 'w')

signatures = open(signatures, 'r')

# Per iteration, read a VCF
for line in signatures:
    
    ## a) Header:
    if line.startswith("#"):
        line = line.rstrip('\n')
        line = line.split(",")
        line.pop(0)

    ## b) Not header:
    else:

        ## I need to replace sample Id by specimen id
        line = line.rstrip('\n')
        line = line.split(",")
        sampleId = line[0]
        specimenId = medatadaDict[sampleId] 
        line[0] = specimenId
    
    ## Write tab separated line into the output file:
    row = '\t'.join(line)

    row = row + '\n'
    outFile.write(row)

## End ##
print
print "***** Finished! *****"
print
