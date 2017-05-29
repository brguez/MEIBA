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
parser.add_argument('contributions', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
frequencies = args.frequencies
contributions = args.contributions
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "frequencies: ", frequencies
print "contributions: ", contributions
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Load input tables:
##########################
header("1. Load input tables")
frequenciesDf = pd.read_csv(frequencies, header=0, index_col=0, sep='\t')

contributionsDf = pd.read_csv(contributions, header=0, index_col=0, sep='\t')

print "contributionsDf: ", contributionsDf

#### 1. Load input tables:
##########################
# Convert allele frequencies from fractions to percentages
frequenciesDf = frequenciesDf * 100

#### 3. Reorder ancestries in alphabetical order 
#################################################
## Make list of ancestries
ancestriesList = frequenciesDf.columns.values.tolist()

ancestriesList.pop(0) ## remove PCAWG group

## Order ancestries in alphabetically order
ancestriesListSorted = sorted(ancestriesList)

## Reorder dataframe
frequenciesSortedDf = frequenciesDf.loc[: , ancestriesListSorted]

#### 4. Select only those source elements contributing at least 1%  
####################################################################

srcElementsList = list(contributionsDf.index) # Remove other category
del srcElementsList[-1]


print srcElementsList
frequenciesSortedFilteredDf = frequenciesSortedDf.loc[srcElementsList , :]

print "frequenciesSortedFilteredDf: ", frequenciesSortedFilteredDf

#### 5. Make plots:
###################
header("5. Make plots")

## Set plotting style
sns.set_style("whitegrid")

# Transpose dataframe
frequenciesSortedFilteredDf = frequenciesSortedFilteredDf.T

### VAF per ancestry heatmap
fig = plt.figure(figsize=(10,2))
fig.suptitle('')

ax = sns.heatmap(frequenciesSortedFilteredDf, vmin=0, vmax=100, annot=True, fmt=".1f", linewidths=.5, cmap=plt.cm.Blues, cbar=True, annot_kws={"size": 8}, square=True)

ax.set_xlabel('')
ax.set_ylabel('')
ax.xaxis.tick_top()

# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(45)

## Save figure 
fileName = outDir + "/Pictures/germline_srcElements_VAF_heatmap.pdf"
plt.savefig(fileName)

####
header("Finished")
