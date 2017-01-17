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

        donorId = line[0]
        tumorType = line[1]
        donorIdProjectCodeDict[donorId] = tumorType

#### 2. Compute parameters 
##########################
# - Total number of transductions per source element and project code. 

## Generate dictionary with the following format:
# - dict: key(sourceElementId) -> dict2: key(projectCode) -> total_number_transductions

header("3. Compute parameters")

nbTransductionsPerTumorDict = {}

# Make donor Ids list
activityFile = open(activity, 'r')
donorIdList = activityFile.readline().rstrip().split("\t")
donorIdList = donorIdList[1:] # remove first list element

for line in activityFile:

    activityList = line.rstrip().split("\t")
    sourceElementId = activityList.pop(0)
    print "** source element ** ", sourceElementId

    ## Initialize dictionary for source element
    nbTransductionsPerTumorDict[sourceElementId] = {}

    # Initialize total number of source element transductions per project code to 0 values:
    for projectCode in donorIdProjectCodeDict.values():
        nbTransductionsPerTumorDict[sourceElementId][projectCode] = 0

    print "Initialized-dict: ", sourceElementId, len(nbTransductionsPerTumorDict), nbTransductionsPerTumorDict
    print "list-lengths (donorIds, activities): ", len(donorIdList), len(activityList)

    ## For each donor source element activity
    for donorIndex, donorId in enumerate(donorIdList):

        srcElementActivity = int(activityList[donorIndex])

        # a) Project code available
        if donorId in donorIdProjectCodeDict:
            
            projectCode = donorIdProjectCodeDict[donorId]
        
            nbTransductionsPerTumorDict[sourceElementId][projectCode] += srcElementActivity
        
            print "Source element activity (donorId, activity, tumorActivity): ", donorId, projectCode, srcElementActivity, nbTransductionsPerTumorDict[sourceElementId][projectCode]
        
        # b) Project code do not available
        else:
            print "[ERROR] No project code available for current donorId"      

#### 3. Make tables:
####################

### Source element total number of transuctions per tumor type 

# Create pandas dataframe from dictionary
df=pd.DataFrame(nbTransductionsPerTumorDict) 

# transpose dictionary to have source elements as rows and tumor types as columns
df=df.T 

# Save output into tsv
outFilePath = outDir + '/germline_source_element_activity_perTumorType.tsv'
df.to_csv(outFilePath, sep='\t') 

####
header("Finished")
