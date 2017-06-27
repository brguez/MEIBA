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
parser.add_argument('nbTD', help='')
parser.add_argument('VAF', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
nbTD = args.nbTD
VAF = args.VAF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "nbTD: ", nbTD
print "VAF: ", VAF
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 



#### 1. Generate dataframe containing for each source element its total number of transductions and the VAF
#############################################################################################################

header("1. Generate dataframe containing for each source element its total number of transductions and the VAF")

## Load input files
nbTdSeries = pd.read_csv(nbTD, header=0, index_col=0, sep='\t')
VAFDf = pd.read_csv(VAF, header=0, index_col=0, sep='\t')
VAFPCAWGSeries = VAFDf["PCAWG"]

VAFPCAWGSeries = VAFPCAWGSeries.rename("VAF")

## Create dataframe
nbTdVAFDf = pd.concat([nbTdSeries, VAFPCAWGSeries], axis=1)

print "nbTdVAFDf: ", nbTdVAFDf

#sys.exit(1)

#### 3. Make the scatter plot
###############################

header("3. Make the scatter plot")

fig = plt.figure(figsize=(5,6))

ax = sns.regplot(x="VAF", y="nbTd", data=nbTdVAFDf, scatter=True, fit_reg=False)

### Axis labels
#ax.set_xlabel('VAF')
#ax.set_ylabel('# Transductions PCAWG')

# turn the axis labels
#for item in ax.get_yticklabels():
#    item.set_rotation(0)

#for item in ax.get_xticklabels():
#    item.set_rotation(90)

## Y ticks
#ax.set(yticks=np.arange(0,26,2))

## Save figure 
fileName = outDir + "/Pictures/nbTd_VAF_scatterPlot.pdf"
plt.savefig(fileName)


#### End
header("FINISH!!")


