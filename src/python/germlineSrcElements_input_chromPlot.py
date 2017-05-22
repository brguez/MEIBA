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
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import seaborn as sns
import scipy

## Get user's input ##
parser = argparse.ArgumentParser(description= "Generate input file for producing chromosome plot with source elements frequency and activity")
parser.add_argument('metadata', help='Source elements metadata info')
parser.add_argument('frequency', help='Source elements frequency info')
parser.add_argument('activity', help='Source elements activity info')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
metadata = args.metadata
frequency = args.frequency
activity = args.activity
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "metadata: ", metadata
print "frequency: ", frequency
print "activity: ", activity
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

### 1. Load input files into dataframes
#########################################
header("1. Load input files into dataframes")

## source elements metadata
metadataDf = pd.read_csv(metadata, header=0, index_col=0, sep='\t')

## source elements frequency
frequencyDf = pd.read_csv(frequency, header=0, index_col=0, sep='\t')

## source elements activity
activityDf = pd.read_csv(activity, header=0, index_col=0, sep='\t')


### 2. Collect relevant info in a dataframe 
#############################################
header("2. Collect relevant info in a dataframe")

colNames = ["chrom", "bkpA"]
finalDf = metadataDf.loc[:, colNames]
finalDf["alleleFreq"] = frequencyDf["PCAWG"]
finalDf["activityRate"] = activityDf["activityRate"]

## Reorder columns:
colNames = ["chrom", "bkpA", "alleleFreq", "activityRate"]
finalOrderedDf = finalDf.loc[:, colNames]

### 3. Print output file 
##########################
header("3. Print output file ")

## Save dataframe into output file
outFilePath = outDir + '/germline_source_elements_input_chromPlot.tsv'
finalOrderedDf.to_csv(outFilePath, sep='\t', index=False) 


#### End
header("FINISH!!")


