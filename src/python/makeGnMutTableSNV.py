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
parser.add_argument('metadata', help='Tabular text file containing donor metadata info')
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

### 1. Read metadata file and create list of included donors 
donorMetadataFile = open(metadata, 'r')
donorIdList = []
sampleIdEqDict = {}

## For donor's metadata:
for line in donorMetadataFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
        
        icgc_donor_id = line[1]
        wgs_exclusion_trafic = line[3]
        histology_exclusion_status = line[11]
        tumor_wgs_specimen_count = line[18]
        tumor_wgs_aliquot_id = line[21]
        
        ### I will discard:
        #   - TraFiC blacklisted donors
        #   - Donors with excluded histology
        #   - Multi-sample donors      
        if (wgs_exclusion_trafic == "Whitelist") and (histology_exclusion_status == "included") and (tumor_wgs_specimen_count == "1") :
        
            donorIdList.append(icgc_donor_id)
            #donorIdList.append(tumor_wgs_aliquot_id)
            sampleIdEqDict[tumor_wgs_aliquot_id] = icgc_donor_id
           

# print len(donorIdList), len(sampleIdEqDict)
# 2681 2681, Ok!

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

#print "protCodingGnList: ", protCodingGnList

### 3. Initialize dictionary containing SNV mutational status for gencode protein coding genes for each donor
mutStatusDict = {}

# For each donor
for donorId in donorIdList:

    mutStatusDict[donorId]  = {}

    # For each protein coding gene
    for geneName in protCodingGnList:

        mutStatusDict[donorId][geneName] = 0


# print len(mutStatusDict)
# 2681

### 4. Read input file with exonic SNV calls and update mutational status dictionary
inputFile = open(inputPath, 'r')

## For donor's metadata:
for line in inputFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
        # print "TIOO: ", line

        tumor_wgs_aliquot_id = line[0]
        geneName = line[6]
        impact = line[14]

        #print tumor_wgs_aliquot_id, geneName, impact

        ## If included donor     
        if (tumor_wgs_aliquot_id in sampleIdEqDict): 
            
            icgc_donor_id = sampleIdEqDict[tumor_wgs_aliquot_id] 

            ## If non-synonymous SNV affecting gencode V19 gene
            if (impact != "Synonymous") and (geneName in mutStatusDict[icgc_donor_id]):        

                mutStatusDict[icgc_donor_id][geneName] = 1
                #mutStatusDict[tumor_wgs_aliquot_id][geneName] = 1


### 5. Convert dictionary into dataframew
#            donorId1     donorId2     donorId3....
# gene1      X1           Y1           Z1     
# gene2      X2           Y2           Z2
# ...

mutStatusDict = mutStatusDict
mutStatusDf =  pd.DataFrame(mutStatusDict)

## Convert to integers:
mutStatusDf = mutStatusDf.apply(pd.to_numeric)

## Remove donors with no non-synonymous mutation in protein coding gene:
mutStatusDf = mutStatusDf.loc[:, (mutStatusDf != 0).any(axis=0)]

# Tmp, select TP53
# mutStatusDf = mutStatusDf.loc[['TP53', 'ATM']]

# Save ordered dataframes into tsv
outFilePath = outDir + '/protein_coding_SNVmut.tsv'
mutStatusDf.to_csv(outFilePath, sep='\t') 


## End ##
print
print "***** Finished! *****"
print
