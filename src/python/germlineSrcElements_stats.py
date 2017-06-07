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

test = transductionsCountDf.loc[: ,"OCCAMS-AH-096"]
outFilePath = outDir + '/test.tsv'
test.to_csv(outFilePath, sep='\t') 

transductionsCountDf = transductionsCountDf.T 


### Compute the number of transductions per donor
nbTdPerDonorSeries = transductionsCountDf.sum(axis=1)

nbTdPerDonorSeries.sort_values(ascending=False, inplace=True)

# Save dataframe into output file
outFilePath = outDir + '/nb_transductions_germline_source_per_donor.tsv'
nbTdPerDonorSeries.to_csv(outFilePath, sep='\t') 


### Compute the number of active source elements per donor
nbActiveSrcElementsPerDonorSeries = transductionsCountDf.applymap(binary).sum(axis=1)

nbActiveSrcElementsPerDonorSeries.sort_values(ascending=False, inplace=True)

# Save dataframe into output file
outFilePath = outDir + '/nb_active_germline_source_elements_per_donor.tsv'
nbActiveSrcElementsPerDonorSeries.to_csv(outFilePath, sep='\t') 

### Average number of active source elements in those samples with retrotransposition activity
nbActiveSrcElementsPerDonorFilteredSeries = nbActiveSrcElementsPerDonorSeries[nbActiveSrcElementsPerDonorSeries > 0]

median = np.median(nbActiveSrcElementsPerDonorFilteredSeries.values)
average = np.average(nbActiveSrcElementsPerDonorFilteredSeries.values)

print "median_nb_active_germline_source_elements: ", median
print "median_nb_active_germline_source_elements: ", average


#### End
header("FINISH!!")


