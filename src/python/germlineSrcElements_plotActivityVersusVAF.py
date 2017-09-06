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
from collections import Counter

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= "Plot the number of active source source elements per tumor genome across each tumor type")
parser.add_argument('activity', help='')
parser.add_argument('VAF', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
activity = args.activity
VAF = args.VAF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "activity: ", activity
print "VAF: ", VAF
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Generate dataframe containing for each source element its activity rate and VAF
########################################################################################
header("1. Generate dataframe containing for each source element its activity rate and VAF")

## Load input files
activityDf = pd.read_csv(activity, header=0, index_col=0, sep='\t')
VAFDf = pd.read_csv(VAF, header=0, index_col=0, sep='\t')

## Create dataframe
dataframe = pd.concat([activityDf["activityRate"], VAFDf["PCAWG"]], axis=1)

#### Bin the allele frequency estimates to variants <1%, 1-5%, and >5%? 
# For each of these 3 bins, what is the average L1 activity
range1 = dataframe[dataframe["PCAWG"] < 0.01]
range2 = dataframe[(dataframe["PCAWG"] >= 0.01) & (dataframe["PCAWG"] <= 0.05)]
range3 = dataframe[dataframe["PCAWG"] > 0.05]

print "*** VAF ***"
print "range1-average: ", range1["activityRate"].mean(axis=0)
print "range2-average: ", range2["activityRate"].mean(axis=0)
print "range3-average: ", range3["activityRate"].mean(axis=0)

print "range1-median: ", range1["activityRate"].median(axis=0)
print "range2-median: ", range2["activityRate"].median(axis=0)
print "range3-median: ", range3["activityRate"].median(axis=0)

#### Bin the L1 activity to variants 1-5, >5-10, and >10? 
# For each of these 3 bins, what is the average VAF
range1 = dataframe[(dataframe["activityRate"] >= 1) & (dataframe["activityRate"] <= 5)]
range2 = dataframe[(dataframe["activityRate"] > 5) & (dataframe["activityRate"] <= 10)]
range3 = dataframe[dataframe["activityRate"] > 10]

print "*** Activity ***"
print "range1-average: ", range1["PCAWG"].mean(axis=0)
print "range2-average: ", range2["PCAWG"].mean(axis=0)
print "range3-average: ", range3["PCAWG"].mean(axis=0)

print "range1-median: ", range1["PCAWG"].median(axis=0)
print "range2-median: ", range2["PCAWG"].median(axis=0)
print "range3-median: ", range3["PCAWG"].median(axis=0)


#### 2. Apply logarithm and filter the dataframe
###############################################
## Add pseudocount and apply logarithm to activity values
pseudocount = 1
dataframe["activityRate"] += pseudocount

dataframe["log10ActivityRate"] = np.log10(dataframe["activityRate"])

print "dataframe: ", dataframe

## Exclude 7 L1PA source elements from the analysis:
L1PAs = ["1p31.1-1", "11p15.2", "15q25.2-1", "18q21.32", "3p21.1", "5p13.1", "6p12.3"]

dataframe = dataframe.drop(L1PAs)

#### 3. Make the scatter plot
###############################

header("3. Make the scatter plot")

fig = plt.figure(figsize=(5,6))

### Compute correlation
corr = stats.spearmanr(dataframe["activityRate"], dataframe["PCAWG"])
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'spearman_corr: ' + str(coefficient) + "; p_value: " + str(pvalue)  

### Make the scatterplot
ax = sns.regplot(x="PCAWG", y="log10ActivityRate", data=dataframe, scatter=True, fit_reg=False)

### Add the p-value
ax.text(0, 1.02, text, transform = ax.transAxes)

### Axis labels
ax.set_xlabel('VAF')
ax.set_ylabel('Transductions/sample')


# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(90)

## Y ticks
ax.set(xlim=(0,1))

## Save figure 
#fileName = outDir + "/Pictures/activityRate_VAF_scatterPlot.ALL.ALL.pdf"
#fileName = outDir + "/Pictures/activityRate_VAF_scatterPlot.ALL.L1HS.pdf"
#fileName = outDir + "/Pictures/activityRate_VAF_scatterPlot.EUR.ALL.pdf"
fileName = outDir + "/Pictures/activityRate_VAF_scatterPlot.EUR.L1HS.pdf"
plt.savefig(fileName)


#### End
header("FINISH!!")


