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
from scipy.stats import spearmanr

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description="")
parser.add_argument('input', help='')
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
inputDf = pd.read_csv(inputFile, header=0, index_col=0, sep='\t')

inputDf = inputDf.convert_objects(convert_numeric=True)

print "inputDf: ", inputDf

#### 2. Number L1 insertions and number of active source elements
#####################################################################

## Add pseudocount to expression values (to avoid expr of 0) 
pseudocount = 1
inputDf["nbActive"] += pseudocount

## Apply logarithm to expression values
inputDf["nbActive"] += pseudocount

inputDf["log10nbActive"] = np.log10(inputDf["nbActive"])
inputDf["log10nbL1"] = np.log10(inputDf["nbL1"])

fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("log10nbL1", "log10nbActive", data=inputDf, kind="kde", space=0, cmap="Blues", stat_func=spearmanr)

#dropna=True, 
# xlim=(0,30), 

## Save figure
outPath = outDir + '/Pictures/nbL1_nbActiveSrcElements_correlation.pdf'
plt.savefig(outPath)


####
header("Finished")

