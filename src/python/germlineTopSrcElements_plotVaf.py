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
parser.add_argument('alleleFreqs', help='')
parser.add_argument('alleleCounts', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
alleleFreqs = args.alleleFreqs
alleleCounts = args.alleleCounts
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "alleleFreqs: ", alleleFreqs
print "alleleCounts: ", alleleCounts
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Load input tables:
##########################
header("1. Load input tables")
frequenciesDf = pd.read_csv(alleleFreqs, header=0, index_col=0, sep='\t')
countsDf = pd.read_csv(alleleCounts, header=0, index_col=0, sep='\t')

# Convert allele frequencies from fractions to percentages
frequenciesDf = frequenciesDf * 100


#### 2. Select hot source L1:
#############################
#hotSourceList = ["22q12.1", "14q23.1", "6p22.1", "6p24.1", "Xp22.2-1", "9q32", "2q24.1", "3q21.1", "Xp22.2-2", "7p12.3", "3q26.1"]
hotSourceList = ["22q12.1", "14q23.1", "6p22.1", "6p24.1", "Xp22.2-1", "9q32", "2q24.1", "3q21.1", "Xp22.2-2", "7p12.3", "3q26.1", "1p12", "8q24.22", "1p31.1-2", "13q21.2-2", "1p22.3", "5q14.3", "7p14.3"]
frequenciesHotDf = frequenciesDf.loc[hotSourceList]
countsHotDf = countsDf.loc[hotSourceList]

print "frequenciesHotDf: ", frequenciesHotDf

print "countsHotDf: ", countsHotDf

#### 3. Make plot:
###################
header("3. Make plot")

## Set plotting style
sns.set_style("whitegrid")

# Initialize the matplotlib figure
fig = plt.figure(figsize=(6, 2))

# Plot source element allele frequency in the complete PCAWG cohort
ax = frequenciesHotDf["PCAWG"].plot(kind='bar')
ax.set(ylim=(0, 100), ylabel="VAF (%)")

# Add to the top of each bar the source elemenet allele count
rects = ax.patches
labels = countsHotDf["PCAWG"].values

for rect, label in zip(rects, labels):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2, height, label, ha='center', va='bottom', size='8')

## Save figure 
fileName = outDir + "/Pictures/topSourceElements_VAF_alleleCounts.pdf"
plt.savefig(fileName)

####
header("Finished")
