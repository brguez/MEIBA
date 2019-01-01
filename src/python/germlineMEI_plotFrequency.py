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
parser.add_argument('alleleCounts', help='tsv with MEI allele counts')
parser.add_argument('alleleFreqs', help='tsv with MEI allele freqs')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
alleleCounts = args.alleleCounts
alleleFreqs = args.alleleFreqs
fileName = args.fileName
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "alleleCounts: ", alleleCounts
print "alleleFreqs: ", alleleFreqs
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Load input table:
##########################
alleleCountsDf = pd.read_csv(alleleCounts, header=0, index_col=0, sep='\t')
alleleFreqsDf = pd.read_csv(alleleFreqs, header=0, index_col=0, sep='\t')

#### 2. Make plots:
###################
header("2. Make plots")

#### 2.1 Make allele counts scatterplot
header("2.1 Make allele counts scatterplot")

### Organize the data for plotting
alleleCountList =  alleleCountsDf['COHORT'].tolist()

alleleCountTuple = [(i, alleleCountList.count(i)) for i in set(alleleCountList)]
tmpList = map(list, zip(*alleleCountTuple))
alleleCountList = tmpList[0]
nbMEIList = tmpList[1]

### Make plot
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1)
plt.scatter(alleleCountList, nbMEIList, color='#87CEFA', edgecolor='black', linewidth='0.5', s=20, alpha=0.7) 
ax.set_xscale('log', basex=10)
ax.set_yscale('log', basex=10)
plt.xlim(0.5, max(alleleCountList) + 100)
plt.ylim(0.5, max(nbMEIList) + 100)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize ticks
# X axis
xPosList = [ 1, 5, 10, 50, 100, 200, 300, 400, 500, 1000, 2500, max(alleleCountList) ]
ax.set_xticks(xPosList)
ax.set_xticklabels(xPosList)
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)

# y axis
yPosList = [ 1, 10, 100, 1000, max(nbMEIList) ]
ax.set_yticks(yPosList)
ax.set_yticklabels(yPosList)

## Save figure
outFile = outDir + '/' + fileName + ".alleleCount2.pdf"
fig.savefig(outFile)

outFile = outDir + '/' + fileName + ".alleleCount2.svg"
fig.savefig(outFile)

outFile = outDir + '/' + fileName + ".alleleCount2.png"
fig.savefig(outFile)

#### 2.2 Make allele frequencies histogram
header("2.2 Make allele frequencies histogram")

### Organize the data for plotting
alleleFreqsList =  alleleFreqsDf['COHORT'].tolist()

## Make plot
fig = plt.figure(figsize=(8,8))
fig.suptitle('Variant allele frequencies', fontsize=20)
ax = fig.add_subplot(1, 1, 1)
plt.hist(alleleFreqsList)
#plt.hist(alleleFreqsList, bins=40, color='#008000', alpha=0.75)
plt.xlabel("VAF", fontsize=14)
plt.ylabel("# MEI", fontsize=14)
plt.xlim(0, 1)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 1.01, 0.1))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

## Save figure
outFile = outDir + '/' + fileName + ".alleleFreq.pdf"
fig.savefig(outFile)


####
header("Finished")

