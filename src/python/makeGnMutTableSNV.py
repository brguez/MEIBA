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
import numpy as np

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='')
parser.add_argument('metadata', help='Tabular text file containing sample metadata info')
parser.add_argument('annot', help='Bed file containing gencode protein coding genes')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
metadata = args.metadata
annot = args.annot
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "metadata: ", metadata
print "annot: ", annot
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ##

### 1. Read metadata file and create list of included samples 
sampleMetadataFile = open(metadata, 'r')
allSampleIdList = []
sampleIdEqDict = {}

## For sample's metadata:
for line in sampleMetadataFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
       
        histology_exclusion_status = line[11]
        sampleIds = line[20]        
        aliquotIds = line[21]
        
        ### I will discard:
        #   - samples with excluded histology
        if (histology_exclusion_status == "included"):
        
            sampleIdList = sampleIds.split(",")
            aliquotIdList = aliquotIds.split(",")

            for index, sampleId in enumerate(sampleIdList):
                aliquotId = aliquotIdList[index]
                sampleIdEqDict[aliquotId] = sampleId

                allSampleIdList.append(sampleId)

### 2. Make list of protein coding genes
annotFile = open(annot, 'r')
protCodingGnList = []

## For each gene
for line in annotFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
        
        geneName = line[3]
        protCodingGnList.append(geneName)


### 3. Initialize dictionary containing SNV mutational status for gencode protein coding genes for each sample
mutStatusDict = {}

# For each sample
for sampleId in allSampleIdList:

    mutStatusDict[sampleId]  = {}

    # For each protein coding gene
    for geneName in protCodingGnList:

        mutStatusDict[sampleId][geneName] = 0

### 4. Read input file with exonic SNV calls and update mutational status dictionary
inputFile = open(inputPath, 'r')
samplesWithSNVList = []

## For sample's metadata:
for line in inputFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")

        tumor_wgs_aliquot_id = line[0]
        geneName = line[6]
        impact = line[14]

        ## If included sample     
        if (tumor_wgs_aliquot_id in sampleIdEqDict): 
            
            icgc_sample_id = sampleIdEqDict[tumor_wgs_aliquot_id] 
            samplesWithSNVList.append(icgc_sample_id)

            ## If non-synonymous SNV affecting gencode V19 gene
            if (impact != "Synonymous") and (geneName in mutStatusDict[icgc_sample_id]):        

                mutStatusDict[icgc_sample_id][geneName] = 1

### 5. Convert dictionary into dataframew
#            sampleId1     sampleId2     sampleId3....
# gene1      X1           Y1           Z1     
# gene2      X2           Y2           Z2
# ...
mutStatusDict = mutStatusDict
mutStatusDf =  pd.DataFrame(mutStatusDict)

## Convert to integers:
mutStatusDf = mutStatusDf.apply(pd.to_numeric)

## Remove samples with no mutational data (no mutation in protein coding genes)
samplesMutDataList = list(set(samplesWithSNVList))
mutStatusDf = mutStatusDf[samplesMutDataList]


# Save ordered dataframes into tsv
outFilePath = outDir + '/proteinCoding_nonSynonomous_SNVmut.tsv'
mutStatusDf.to_csv(outFilePath, sep='\t') 


## End ##
print
print "***** Finished! *****"
print
