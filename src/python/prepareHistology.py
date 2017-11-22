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
import pandas as pd

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: donorId vcfPath')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## 
###############################################

inputFile = open(inputPath, 'r')

histologyDict = {}

# Per iteration, read a VCF
for line in inputFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
        donorUniqueId = line[0]
        histology = line[1]       
        tier1 = line[2]
        tier2 = line[3]
    
        if (donorUniqueId not in histologyDict):
            histologyDict[donorUniqueId] = {}
            histologyDict[donorUniqueId]['histologyList'] = [ histology ]
            histologyDict[donorUniqueId]['tier1List'] = [ tier1 ]
            histologyDict[donorUniqueId]['tier2List'] = [ tier2 ]

        else:
            histologyDict[donorUniqueId]['histologyList'].append(histology)
            histologyDict[donorUniqueId]['tier1List'].append(tier1)
            histologyDict[donorUniqueId]['tier2List'].append(tier2)


histologyOkDict = {}

for donorUniqueId in histologyDict:

    histologyOkDict[donorUniqueId] = {}
    histologyOkDict[donorUniqueId]['histology_abbreviation'] = ",".join(set(histologyDict[donorUniqueId]['histologyList']))
    histologyOkDict[donorUniqueId]['histology_tier1'] = ",".join(set(histologyDict[donorUniqueId]['tier1List']))
    histologyOkDict[donorUniqueId]['histology_tier2'] = ",".join(set(histologyDict[donorUniqueId]['tier2List']))

histologyOkDf = pd.DataFrame(histologyOkDict) 
histologyOkDf = histologyOkDf.T

outFilePath = outDir + '/PCAWG_histology.ok.tsv'
histologyOkDf.to_csv(outFilePath, sep='\t') 




