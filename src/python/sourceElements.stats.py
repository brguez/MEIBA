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
parser.add_argument('inputVCF', help='multi-sample VCF file containing genotyped source elements')
parser.add_argument('ancestryFile', help='text file with the ancestry code per donor Id')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
ancestryFile = args.ancestryFile
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "ancestryFile: ", ancestryFile
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read input multi-sample VCF and generate a VCF object
#############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Read ancestry codes file
##################################
# Make two dictionary with the following structure:
# - dict1: key(ancestryCode) -> dict2: key(donorId) -> nbSourceElements
# - dict3: key(donorId) -> -> ancestryCode

header("2. Read ancestry codes file")
ancestryFile = open(ancestryFile, 'r')

nbSourceElementsDict = {}
donorIdAncestryDict = {}

for line in ancestryFile:

    line = line.rstrip('\r\n')

    #print "line", line

    if not line.startswith("#"):
        line = line.split('\t')

        donorId = line[0]
        ancestryCode = line[1]

        # print 'values: ', donorId, ancestryCode

        ## Dict1a and Dict1b
        if ancestryCode not in nbSourceElementsDict:

            # Create dictionary
            nbSourceElementsDict[ancestryCode] = {}


        # Initialize to 0 values:
        nbSourceElementsDict[ancestryCode][donorId] = 0

        ## dict3
        donorIdAncestryDict[donorId] = ancestryCode


#### 3. Compute parameters:
###########################
alleleFreqDict = {}

header("3. Compute parameters")

for MEIObj in VCFObj.lineList:

    sourceElementId = MEIObj.chrom + '_' + str(MEIObj.pos) + '_' + MEIObj.infoDict["CLASS"]

    ## Total number of chromosome copies in the population
    # Number of donors * 2 (diploid, two copies of a given chromosome)
    totalNbChrom = len(MEIObj.genotypesDict) * 2

    ## Compute MEI allele frequency per L1 source element
    alleleCount = 0

    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        ancestryCode = donorIdAncestryDict[donorId]
        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]
        nbReadsMEI = float(genotypeFieldList[1])
        totalNbReads = float(genotypeFieldList[2])

        ## Update counters and store VAF values
        # A) Heterozygous
        if (genotype == "0/1"):
            nbSourceElementsDict[ancestryCode][donorId] += 1
            alleleCount +=  1

        # B) Homozygous
        elif (genotype == "1/1"):
            nbSourceElementsDict[ancestryCode][donorId] += 1
            alleleCount += 2

    ## Compute MEI allele frequency:
    alleleFrequency = float(alleleCount) / totalNbChrom

    ## Save into list. One per MEI type
    alleleFreqDict[sourceElementId] = alleleFrequency

#### 4. Make plots:
#####################
header("4. Make plots")

# - Variant allele frequencies histogram across PCAWG donors (done)
# - Number of source elements per donor and ancestry. Boxplot (done)
# - Number of source elements per donor and ancestry. Violin Plot (done)

#### 4.1 Source element variant allele frequencies across PCAWG donors
header("4.1 Make variant allele frequencies plot")

alleleFreqList = alleleFreqDict.values()

fig = plt.figure(figsize=(5,6))
fig.suptitle('Variant allele frequencies (VAF)', fontsize=14)

## Make plot
ax1 = fig.add_subplot(1, 1, 1)
plt.hist(alleleFreqList, bins=20, color='#008000', alpha=0.75)
plt.xlabel("VAF", fontsize=14)
plt.ylabel("# L1 source elements", fontsize=12)
plt.ylim(0, 11)
plt.xlim(0, 1)

# Remove top and right axes
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 1.01, 0.1))
locs, labels = plt.xticks()
# plt.setp(labels, rotation=30)

## Save figure
fileName = outDir + "/PCAWG_sourceElements_VAF_hist.pdf"
plt.savefig(fileName)

#### 4.2 Number of source elements per donor and ancestry
header("4.2 Number of source elements per donor and ancestry")

### A) Boxplot
## Organize the data for plotting
tupleListNbSourceElements = []

for ancestryCode in sorted(nbSourceElementsDict):

    # Make tuple (ancestryCode, list with number of source elements per donor)
    nbDonors= len(nbSourceElementsDict[ancestryCode].values())
    xLabel = ancestryCode + '(' +  str(nbDonors) + ')'
    nbSourceElementsTuple = (xLabel, nbSourceElementsDict[ancestryCode].values())

    # Add tuple to the list
    tupleListNbSourceElements.append(nbSourceElementsTuple)


## Make nested list with the following format:
# [donor1_nbSourceElements, donor2_nbSourceElements, ..., donorN_nbSourceElements], [donor1_nbSourceElements, donor2_nbSourceElements, ..., donorN_nbSourceElements] , ... [donor1_nbSourceElements, donor2_nbSourceElements, ..., donorN_nbSourceElements]
#                               ancestry1_list                                                                  ancestry2_list                                                                          ancestryN_list

tmpList = map(list, zip(*tupleListNbSourceElements))
ancestryCodesList = tmpList[0]
nbSourceElementsPerDonor = tmpList[1]

### Plotting
fig = plt.figure(figsize=(5,6))
fig.suptitle('# Source elements per donor', fontsize=14)

ax1 = fig.add_subplot(1, 1, 1)

# Create the boxplot
bp = ax1.boxplot(nbSourceElementsPerDonor)
plt.ylabel("# Source L1", fontsize=12)

## Customize boxplot:
# change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#696969', linewidth=1)

# change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#696969', linewidth=1)

# change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#696969', linewidth=1)

# change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#8b0000', linewidth=2)

# Add the ancestry codes to the x-axis
ax1.set_xticklabels(ancestryCodesList, fontsize = 10)
locs, labels = plt.xticks()
# plt.setp(labels, rotation=25)

# Remove top and right axes
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_axisbelow(True)

## Save figure
fileName = outDir + "/PCAWG_nbSourceElementsPerDonor_boxplot.pdf"
fig.savefig(fileName)

### B) Violin plot
## Organize the data for plotting into a dictionary:
# - dict1:
#       nbSourceElements -> list[nbSourceElementsDonor1, nbSourceElementsDonor2, ..., nbSourceElementsDonorN]
#       ancestryPerDonor -> list[ancestryDonor1, ancestryDonor2, ..., nbSourceElementsDonorN]

dict4pandas = {}
dict4pandas['nbSourceElements'] = []
dict4pandas['ancestry'] = []

for ancestryCode in sorted(nbSourceElementsDict):

    nbSourceElementsPerDonorList = nbSourceElementsDict[ancestryCode].values()
    nbDonors = len(nbSourceElementsPerDonorList)

    xLabel = ancestryCode + '(' +  str(nbDonors) + ')'
    xLabelList = [ xLabel ] * nbDonors

    dict4pandas['nbSourceElements'] = dict4pandas['nbSourceElements'] + nbSourceElementsPerDonorList
    dict4pandas['ancestry'] = dict4pandas['ancestry'] + xLabelList

# Make pandas dataframe from dict:
dataframe = pd.DataFrame(dict4pandas)

### Plotting

fig = plt.figure(figsize=(5,6))
fig.suptitle('# Source elements per donor', fontsize=14)

# Create the violin plot
ax = sns.violinplot(x='ancestry', y='nbSourceElements', data=dataframe, palette="muted")

# y limit
sns.plt.ylim(0,21)

## Modify axis labels
ax.set(xlabel='', ylabel='# Source L1')

# Remove top and right axes
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax.set_axisbelow(True)

## Save figure
fileName = outDir + "/PCAWG_nbSourceElementsPerDonor_violinplot.pdf"
fig.savefig(fileName)

####
header("Finished")
