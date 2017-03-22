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
import matplotlib.patches as mpatches
import scipy.stats as stats

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('counts', help='')
parser.add_argument('activeSrc', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
counts = args.counts
activeSrc = args.activeSrc
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "counts: ", counts
print "activeSrc: ", activeSrc
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read file with the number of retrotransposition events per donor
#########################################################################

header("1. Read file with the number of retrotransposition events per donor")
countsFile = open(counts, 'r')

countsActiveSrcDict = {}

for line in countsFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')
       
        donorId = line[0]
        totalNbEvents = int(line[3])
        totalNbL1 = int(line[4])
        totalNbTD = int(line[5])

        countsActiveSrcDict[donorId] = {}
        countsActiveSrcDict[donorId]["totalNbEvents"] = totalNbEvents
        countsActiveSrcDict[donorId]["totalNbL1"] = totalNbL1
        countsActiveSrcDict[donorId]["totalNbTD"] = totalNbTD

#### 2. Read file with the number of active source elements per donor
#############################################################################

header("2. Read file with the number of active source elements per donor")

## Read input file
activeSrcFile = open(activeSrc, 'r')

for line in activeSrcFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')
        
        donorId = line[0]
        nbActiveSrcElements = int(line[1])
        countsActiveSrcDict[donorId]["nbActiveSrcElements"] = nbActiveSrcElements


#### 3. Generate output file with the number of active source elements per donor
#################################################################################
    
header("3. Compute correlations")

##### 3.1 Prepare data for plotting
## Convert into a dataframe
countsActiveSrcDataframe =  pd.DataFrame(countsActiveSrcDict).T

## Filter dataframe to select only those donors with at least one retrotransposition event. 
countsActiveSrcDataframeFiltered1 = countsActiveSrcDataframe[countsActiveSrcDataframe['totalNbEvents'] > 0]
countsActiveSrcDataframeFiltered2 = countsActiveSrcDataframe[countsActiveSrcDataframe['totalNbL1'] > 0]
countsActiveSrcDataframeFiltered3 = countsActiveSrcDataframe[countsActiveSrcDataframe['totalNbTD'] > 0]


nbActiveSrcElementList1 = countsActiveSrcDataframeFiltered1['nbActiveSrcElements'].tolist()
totalNbEventsList = countsActiveSrcDataframeFiltered1['totalNbEvents'].tolist()

nbActiveSrcElementList2 = countsActiveSrcDataframeFiltered2['nbActiveSrcElements'].tolist()
totalNbL1List = countsActiveSrcDataframeFiltered2['totalNbL1'].tolist()

nbActiveSrcElementList3 = countsActiveSrcDataframeFiltered3['nbActiveSrcElements'].tolist()
totalNbTDList = countsActiveSrcDataframeFiltered3['totalNbTD'].tolist()


### 3.2 Assess correlation between number of active source elements and the total number of events

# Compute correlation
corr = stats.pearsonr(nbActiveSrcElementList1, totalNbEventsList)
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient) + "; p_value: " + str(pvalue)  

# Make scatterplot
fig = plt.figure(figsize=(6,6))

ax1 = fig.add_subplot(1, 1, 1)
#ax1.set_title("SVA", fontsize=14)
plt.scatter(nbActiveSrcElementList1, totalNbEventsList, color='#008000', alpha=.4)
plt.xlim((-1, (max(nbActiveSrcElementList1) + 1)))
plt.ylim((0, (max(totalNbEventsList) + 10)))
plt.xlabel('Number of active germline source elements', fontsize=12)
plt.ylabel('Number of retrotransposition events', fontsize=12)
ax1.text(0, 1.02, text, transform = ax1.transAxes)

## Add trend line:
plt.plot(np.unique(nbActiveSrcElementList1), np.poly1d(np.polyfit(nbActiveSrcElementList1, totalNbEventsList, 1))(np.unique(nbActiveSrcElementList1)), color='#000000',)

## Save figure
fileName = outDir + "/nbActiveSrcElement_totalNbEvents_correlation.pdf"
plt.savefig(fileName)


### 3.3 Assess correlation between number of active source elements and number of L1 events

# Compute correlation
corr = stats.pearsonr(nbActiveSrcElementList2, totalNbL1List)
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient) + "; p_value: " + str(pvalue)  

# Make scatterplot
fig = plt.figure(figsize=(6,6))

ax1 = fig.add_subplot(1, 1, 1)
#ax1.set_title("SVA", fontsize=14)
plt.scatter(nbActiveSrcElementList2, totalNbL1List, color='#008000', alpha=.4)
plt.xlim((-1, (max(nbActiveSrcElementList2) + 1)))
plt.ylim((0, (max(totalNbL1List) + 10)))
plt.xlabel('Number of active germline source elements', fontsize=12)
plt.ylabel('Number of L1 events', fontsize=12)
ax1.text(0, 1.02, text, transform = ax1.transAxes)

## Add trend line:
plt.plot(np.unique(nbActiveSrcElementList2), np.poly1d(np.polyfit(nbActiveSrcElementList2, totalNbL1List, 1))(np.unique(nbActiveSrcElementList2)), color='#000000',)

## Save figure
fileName = outDir + "/nbActiveSrcElement_totalNbL1_correlation.pdf"
plt.savefig(fileName)


### 3.4 Assess correlation between number of active source elements and the total number of L1 transductions

# Compute correlation
corr = stats.pearsonr(nbActiveSrcElementList3, totalNbTDList)
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient) + "; p_value: " + str(pvalue)  

# Make scatterplot
fig = plt.figure(figsize=(6,6))

ax1 = fig.add_subplot(1, 1, 1)
#ax1.set_title("SVA", fontsize=14)
plt.scatter(nbActiveSrcElementList3, totalNbTDList, color='#008000', alpha=.4)
plt.xlim((-1, (max(nbActiveSrcElementList3) + 1)))
plt.ylim((0, (max(totalNbTDList) + 10)))
plt.xlabel('Number of active germline source elements', fontsize=12)
plt.ylabel('Number of L1-transductions', fontsize=12)
ax1.text(0, 1.02, text, transform = ax1.transAxes)

## Add trend line:
plt.plot(np.unique(nbActiveSrcElementList3), np.poly1d(np.polyfit(nbActiveSrcElementList3, totalNbTDList, 1))(np.unique(nbActiveSrcElementList3)), color='#000000',)

## Save figure
fileName = outDir + "/nbActiveSrcElement_totalNbTD_correlation.pdf"
plt.savefig(fileName)


####
header("Finished")


