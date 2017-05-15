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

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputFile', help='tsv with the number of different MEI per donor')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputFile = args.inputFile
fileName = args.fileName
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputFile: ", inputFile
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Load input table:
##########################
loadDf = pd.read_csv(inputFile, header=0, index_col=0, sep='\t')

#### 2. Make plots:
###################
header("2. Make plots")

#### 2.1 Number of elements per donor and ancestry
## Remove donors with unkown ancestry:
loadFilteredDf = loadDf[loadDf['ancestry'] != "UNK"]

## Order ancestries in alphabetically order
orderList = sorted(set(loadFilteredDf['ancestry'].tolist()))

### Plotting
fig = plt.figure(figsize=(5,6))

# Create the violin plot
ax = sns.violinplot(x='ancestry', y='nbMEI', data=loadFilteredDf, palette="muted", order=orderList)

# y limit
#sns.plt.ylim(0,21)

## Modify axis labels
ax.set(xlabel='', ylabel='# Events')
plt.xticks(rotation=60)

# Remove top and right axes
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax.set_axisbelow(True)

## Save figure
outFile = outDir + '/' + fileName + ".ancestries.pdf"
fig.savefig(outFile)



#### 2.2 Number of elements per donor and project code
## Order project codes in alphabetically order
orderList = sorted(set(loadDf['projectCode'].tolist()))

### Plotting

fig = plt.figure(figsize=(18,6))

# Create the violin plot
ax = sns.violinplot(x='projectCode', y='nbMEI', data=loadDf, palette="muted", order=orderList)

# y limit
#sns.plt.ylim(0,21)

## Modify axis labels
ax.set(xlabel='', ylabel='# Events')
plt.xticks(rotation=60)

# Remove top and right axes
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax.set_axisbelow(True)

## Save figure
outFile = outDir + '/' + fileName + ".projectCodes.pdf"
fig.savefig(outFile)



#### 2.3 Number of elements per donor and histology tumor type
## Order histology tumour types in alphabetically order
orderList = sorted(set(loadDf['tumorType'].tolist()))

### Plotting

fig = plt.figure(figsize=(16,5))

# Create the violin plot
ax = sns.violinplot(x='tumorType', y='nbMEI', data=loadDf, palette="muted", order=orderList)

# y limit
#sns.plt.ylim(0,21)

## Modify axis labels
ax.set(xlabel='', ylabel='# Events')
plt.xticks(rotation=60)

# Remove top and right axes
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax.set_axisbelow(True)

## Save figure
outFile = outDir + '/' + fileName + ".histologies.pdf"
fig.savefig(outFile)

####
header("Finished")
