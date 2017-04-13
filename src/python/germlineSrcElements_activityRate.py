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


def binary(value):
    """
    Convert value into binary
    """

    # A) Value equal to 0
    if (value == 0): 
        boolean = 0

    # B) Value higher than 0
    else:
        boolean = 1

    return boolean


def activityStatus(activityRate):
    """
    Classify source element according its activity rate in one of these categories:
    Activity ranges:
        - none: 0
        - low: (0-3]
        - moderate: (3-9]
        - high: >9
    """

    # a) None: 0
    if (activityRate == 0):
        status = "none"        
    
    # b) Low: (0-3] 
    elif (activityRate > 0) and (activityRate <= 3):
        status = "low"

    # c) Moderate: (3-9]
    elif (activityRate > 3) and (activityRate <= 9):
        status = "moderate"
    
    # d) Hot: >9
    else:
        status = "high"  
    
    return status

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import seaborn as sns
import scipy

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('transductionsCountFile', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
transductionsCountFile = args.transductionsCountFile
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "transductionsCountFile: ", transductionsCountFile
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

### Read transductions count dataframe
transductionsCountDf = pd.read_csv(transductionsCountFile, header=0, index_col=0, sep='\t')

### Initialize activity rate dataframe
srcElementIds = transductionsCountDf.index
activityRateDf = pd.DataFrame(index=srcElementIds)

## Add the total number of transductions per source element
activityRateDf["totalNbTd"] = transductionsCountDf.sum(axis=1)

## Add the total number of active donors per source element
activeDonorBinaryDf = transductionsCountDf.applymap(binary)
activityRateDf["nbActiveDonors"] = activeDonorBinaryDf.sum(axis=1)

## Add the activity rate:
activityRateDf["activityRate"] = activityRateDf["totalNbTd"]/activityRateDf["nbActiveDonors"]

# for those elements active in 0 donors. Set the activity rate as 0
activityRateDf["activityRate"] = activityRateDf["activityRate"].fillna(0)

## Add the acticity status:
activityRateDf['activityStatus'] = activityRateDf["activityRate"].apply(activityStatus)

## Add the max activity values
activityRateDf["maxActivity"] = transductionsCountDf.max(axis=1)

## Save dataframe into output file
outFilePath = outDir + '/germline_srcElements_activityRate.tsv'
activityRateDf.to_csv(outFilePath, sep='\t') 


#### End
header("FINISH!!")


