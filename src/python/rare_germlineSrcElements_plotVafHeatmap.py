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
import scipy.stats as stats
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('frequencies', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
frequencies = args.frequencies
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "frequencies: ", frequencies
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Load input tables:
##########################
header("1. Load input tables")
frequenciesDf = pd.read_csv(frequencies, header=0, index_col=0, sep='\t')

# Convert allele frequencies from fractions to percentages
frequenciesDf = frequenciesDf * 100

print "frequenciesDf: ", frequenciesDf

#### 2. Select only rare source elements (VAF <1%)   
###################################################
header("2. Select only rare source elements (VAF <1%)")
frequenciesFilteredDf = frequenciesDf[frequenciesDf["PCAWG"] < 1]

print "frequenciesFilteredDf: ", frequenciesFilteredDf


#### 3. Make plots:
###################
header("3. Make plots")

## Set plotting style
sns.set_style("whitegrid")

# Transpose dataframe
frequenciesFilteredDf = frequenciesFilteredDf.T

print "frequenciesFilteredDf: ", frequenciesFilteredDf

### VAF per ancestry heatmap
fig = plt.figure(figsize=(10,2))
fig.suptitle('')

ax = sns.heatmap(frequenciesFilteredDf, vmin=0, vmax=1, annot=True, fmt=".2f", linewidths=.5, cmap=plt.cm.Blues, cbar=True, annot_kws={"size": 8}, square=True)

ax.set_xlabel('')
ax.set_ylabel('')
ax.xaxis.tick_top()

# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(90)

## Save figure 
fileName = outDir + "/Pictures/rare_germline_srcElements_VAF_heatmap.pdf"
plt.savefig(fileName)

####
header("Finished")
