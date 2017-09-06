#!/usr/bin/env python
#coding: utf-8

#### MAIN ####

## Import modules ##
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
parser.add_argument('activityRateFile', help='')
parser.add_argument('alleleFreqFile', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
activityRateFile = args.activityRateFile
alleleFreqFile = args.alleleFreqFile
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "activityRateFile: ", activityRateFile
print "alleleFreqFile: ", alleleFreqFile
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1) Make activity rate dataframe
#####################################
dfActivityRate = pd.read_csv(activityRateFile, header=0, index_col=0, sep='\t')

print "dfActivityRate: ", dfActivityRate

#### 2) Make allele frequency across ancestries dataframe
##########################################################
## Load allele count matrix into dataframe

dfAlleleFreqAncestry = pd.read_csv(alleleFreqFile, header=0, index_col=0, sep='\t')

print "dfAlleleFreqAncestry: ", dfAlleleFreqAncestry 

#### 3) none/low vs moderate/high activity VAF comparison
############################################################

### Prepare Data

dataframe = pd.concat([dfActivityRate["activityStatus"], dfAlleleFreqAncestry["PCAWG"]], axis=1)

print "dataframe: ", dataframe

dataframe1 = dataframe[(dataframe['activityStatus']=="none") | (dataframe['activityStatus']=="low")]
dataframe1['activityStatus'] = "none/low"

dataframe2 = dataframe[(dataframe['activityStatus']=="moderate") | (dataframe['activityStatus']=="high")]
dataframe2['activityStatus'] = "moderate/high"
dataframe3 = pd.concat([dataframe1, dataframe2])

x = dataframe3[(dataframe3['activityStatus']=="none/low")]["PCAWG"].values
y = dataframe3[(dataframe3['activityStatus']=="moderate/high")]["PCAWG"].values

mannWhitneyU = stats.mannwhitneyu(x, y, alternative='two-sided')

Ustat = mannWhitneyU[0]
pvalue = round(mannWhitneyU[1], 4) 

## Make boxplot:
#################
title = "Mann-Whitney U = " + str(Ustat) + ", P = " + str(pvalue)

fig = plt.figure(figsize=(5,6))
fig.suptitle(title, fontsize=14)

ax = sns.boxplot(x="activityStatus", y="PCAWG",data=dataframe3, width=0.5, showfliers=False)
ax = sns.stripplot(x="activityStatus", y="PCAWG", data=dataframe3, jitter=True, color=".3")
ax.set(ylim=(-0.1, 1.1))
ax.set_xlabel('')

# Add mann whitney U statistic and p-value to the plot:

## Save figure 
fileName = outDir + "/low_vs_high_activity_VAF_boxPlot.pdf"
plt.savefig(fileName)

sys.exit(1)



