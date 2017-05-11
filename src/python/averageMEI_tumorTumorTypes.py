#!/usr/bin/env python
#coding: utf-8

## Load modules/libraries
import sys
import argparse
import os
import numpy as np
import pandas as pd

## Get user's input
parser = argparse.ArgumentParser(description= "")
parser.add_argument('counts', help='')
parser.add_argument('metadata', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
counts = args.counts
metadata = args.metadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "counts: ", counts
print "metadata: ", metadata
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print


### 1) 
# Make dictionary with donorId as key and tumor types as values

metadata = open(metadata, 'r')
tumorTypeDict = {}

# Read file line by line
for line in metadata:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
       
        donorId = fieldsList[1]
        tumorType = fieldsList[9].split('-')[0]
        tumorTypeDict[donorId] = tumorType

#print tumorTypeDict

### 2) 

counts = open(counts, 'r')
countsDict = {}

# Read file line by line
for line in counts:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
        donorId = fieldsList[1]
        nbMEI = float(fieldsList[2])
        tumorType = tumorTypeDict[donorId]
        
        if tumorType not in countsDict:
            countsDict[tumorType] = []
            
        countsDict[tumorType].append(nbMEI)


averageDict = {}

for tumorType in countsDict:
    average = np.mean(countsDict[tumorType])
    averageDict[tumorType] = average
 
averageSeries = pd.Series(averageDict)
averageSeries.sort(ascending=False)

outFilePath = outDir + '/averageMEI_tumorTypes.tsv'
averageSeries.to_csv(outFilePath, sep='\t') 


print "***** Finished! *****"
print

