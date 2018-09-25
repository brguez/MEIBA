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

def scatterCorr(arrayA, arrayB, threshold, outPath):
    """

        Interpretation of strength of correlation

        very weak: < 0,15 
        weak: 0,15-0,25  
        moderate: 0,25-0,40 
        strong: 0,40-0,75
        very strong: >0,75

    """
    corr = stats.spearmanr(arrayA, arrayB)
    coefficient = float(format(corr[0], '.3f'))
    pvalue = float(corr[1])
    print "pvalue: ", pvalue

    ## Make scatterplot if rho >= threshold or <= -theshold
    if (coefficient >= threshold) or (coefficient <= -threshold):
    
        # Make scatterplot
        fig = plt.figure(figsize=(6,6))
        ax1 = fig.add_subplot(1, 1, 1)
        #plot = sns.jointplot(x=arrayA, y=arrayB, kind="hex", xlim=(0,40), gridsize=50, dropna=True, cmap="Blues", stat_func=spearmanr)
        plot = sns.jointplot(x=arrayA, y=arrayB, kind="kde", space=0, xlim=(0,30), gridsize=50, dropna=True, cmap="Blues", stat_func=spearmanr)
        plt.xlabel('# L1', fontsize=12)
        plt.ylabel('Replication time', fontsize=12)

#        sns.plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  # shrink fig so cbar is visible
#        cax = plot.fig.add_axes([.85, .25, .05, .4])  # x, y, width, height
#        sns.plt.colorbar(cax=cax)

        #sns.jointplot(x=arrayA, y=arrayB, kind="kde", space=0, color="b", xlim=(0,30))



        ## Save figure
        fileName = outPath + '_' + str(coefficient) + '_correlation.pdf' 

        plt.savefig(fileName)

    return coefficient, pvalue


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
from scipy.stats import spearmanr

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")
#sns.set(font="Verdana")

## Get user's input ##
parser = argparse.ArgumentParser(description="Compute correlation between L1 retrotransposition rate and diverse genomic features (replication time, gene expression and gene density)")
parser.add_argument('input', help='Genomic features table')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputFile = args.input
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputFile: ", inputFile
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Load input tables:
##########################
inputDf = pd.read_csv(inputFile, header=0, sep='\t')

print "inputDf: ", inputDf

#### 2. Number L1 insertions and L1 endonuclease motif correlation
###################################################################
## Remove bins with a number of L1 EN motifs == 0
# These will mostly  correspond to telomeric, centromeric regions
filteredDf = inputDf[inputDf["nbL1Motif"] > 0]

## Make plot
fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("nbL1", "nbL1Motif", data=filteredDf, xlim=(0,30), kind="kde", space=0, dropna=True, cmap="Blues", stat_func=spearmanr)

## Save figure
outPath = outDir + '/nbL1_nbL1Motif_corr.pdf'
plt.savefig(outPath)

outPath = outDir + '/nbL1_nbL1Motif_corr.svg'
plt.savefig(outPath)


#### 3. Number L1 insertions and median replication time correlation
#####################################################################
## Remove bins with NA values for replication timing (one bin with NA expression)
filteredDf = inputDf.dropna(subset=["medianRT"]) 

## Make plot
fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("nbL1", "medianRT", data=filteredDf, xlim=(0,30), kind="kde", space=0, dropna=True, cmap="Blues", stat_func=spearmanr)

## Save figure
outPath = outDir + '/nbL1_medianRT_corr.pdf'
plt.savefig(outPath)

outPath = outDir + '/nbL1_medianRT_corr.svg'
plt.savefig(outPath)



####
header("Finished")

