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

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('input', help='tsv...')
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


#### 2. Number L1 insertions and median replication time correlation
#####################################################################
fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("nbL1", "medianRT", data=inputDf, xlim=(0,30), kind="kde", space=0, dropna=True, cmap="Blues", stat_func=spearmanr)
#sns.regplot("nbL1", "medianRT", data=inputDf, ax=plot.ax_joint, scatter=False)

## Save figure
outPath = outDir + '/nbL1_medianRT_correlation.pdf'
plt.savefig(outPath)

#### 3. Number L1 insertions and median expression correlation
###############################################################

## Remove bins with NA values for gene expression (no bins with NA expression actually)
filteredDf = inputDf.dropna(subset=["medianExpr"]) 

## Add pseudocount to expression values (to avoid expr of 0) 
pseudocount = 10
filteredDf["medianExpr"] += pseudocount

## Apply logarithm to expression values
filteredDf["medianExpr"] += pseudocount

filteredDf["log10MedianExpr"] = np.log10(filteredDf["medianExpr"])

## Make plot
fig = plt.figure(figsize=(6,6))
ax2 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("nbL1", "log10MedianExpr", data=filteredDf, xlim=(0,30), ylim=(1,7), shade_lowest=True, kind="kde", space=0, dropna=True, cmap="Blues", stat_func=spearmanr)
#sns.regplot("nbL1", "log10MedianExpr", data=filteredDf, ax=plot.ax_joint, scatter=False)

## Save figure
outPath = outDir + '/nbL1_medianExpression_correlation.pdf'
plt.savefig(outPath)

#### 4. Number L1 insertions and gene density correlation
#####################################################################
fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("nbL1", "geneDensity", data=inputDf, xlim=(0,30), kind="kde", space=0, dropna=True, cmap="Blues", stat_func=spearmanr)
#sns.regplot("nbL1", "medianRT", data=inputDf, ax=plot.ax_joint, scatter=False)

## Save figure
outPath = outDir + '/nbL1_geneDensity_correlation.pdf'
plt.savefig(outPath)

####
header("Finished")

