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
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('activity', help='')
parser.add_argument('metadata', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
activity = args.activity
metadata = args.metadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", activity
print "metadata: ", metadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read metadata file
##########################
# Initialize a dictionary with the following structure:
# - dict: key(donorId) -> projectCode
header("2. Read ancestry codes file")
metadataFile = open(metadata, 'r')

donorIdProjectCodeDict = {}

for line in metadataFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')
        
        exclusion = line[2]
        
        if (exclusion == "Whitelist"):
            donorId = line[0]
            tumorType = line[5].split('-')[0]
            donorIdProjectCodeDict[donorId] = tumorType

#### 2. Make a dictionary with the list of active source elements per donor
#############################################################################
## Generate a dictionary with the following format:
# - dict1: key(donorId) -> number_active_source_elements
# - dict2: key(projectCode) -> active_source_elements_list ** not redundant list

header("2. Compute parameters")

## Read input file
activityDf = pd.read_csv(activity, header=0, index_col=0, sep='\t')
activityDf = activityDf.T

outFilePath = outDir + '/source_elements_activity_transposed.tsv'
activityDf.to_csv(outFilePath, sep='\t') 

nbActiveSrcElementDonorDict = {}
activeSrcElementTumorDict = {}

## For each donor
for donorId, rowSerie in activityDf.iterrows():

    ## Only consider not excluded donors
    if donorId in donorIdProjectCodeDict:

        projectCode = donorIdProjectCodeDict[donorId]
        nbActiveSrcElementDonorDict[donorId] = 0
        
        # Initialize list if first donor from a given tumor type
        if projectCode not in activeSrcElementTumorDict:
            activeSrcElementTumorDict[projectCode] = []
    
        ## for each source element
        for sourceElement, activity in rowSerie.iteritems():

            ## Only consider active source elements
            if (activity>0):

                nbActiveSrcElementDonorDict[donorId] += 1
                
                ## Source not already reported as active in the current tumor type
                if sourceElement not in activeSrcElementTumorDict[projectCode]:
                    activeSrcElementTumorDict[projectCode].append(sourceElement)

#### 3. Generate output file with the number of active source elements per donor
#################################################################################
    
nbActiveSrcElementSeries =  pd.Series(nbActiveSrcElementDonorDict)

## Save output into tsv
outFilePath = outDir + '/donorId_nbActive_srcElements.tsv'
nbActiveSrcElementSeries.to_csv(outFilePath, sep='\t') 


#### 4. Generate output file with the number of distint active source elements per tumor type
###############################################################################################

nbActiveSrcElementTumorDict = {}

for projectCode in activeSrcElementTumorDict:

    nbActiveSrcElementTumorDict[projectCode] = len(activeSrcElementTumorDict[projectCode])    

nbActiveSrcElementTumorSeries =  pd.Series(nbActiveSrcElementTumorDict)

## Save output into tsv
outFilePath = outDir + '/tumorType_nbActive_srcElements.tsv'
nbActiveSrcElementTumorSeries.to_csv(outFilePath, sep='\t') 


####
header("Finished")


