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
parser = argparse.ArgumentParser(description="Compute correlation between L1 retrotransposition rate and the number of active L1 source elements")
parser.add_argument('rtCounts', help='Retrotransposition couns per sample')
parser.add_argument('nbActive', help='Number of active source elements per sample')
parser.add_argument('palette', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
rtCounts = args.rtCounts
nbActive = args.nbActive
palette = args.palette
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "rtCounts: ", rtCounts
print "nbActive: ", nbActive
print "palette: ", palette
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### Read palette file
# Initialize a dictionary with the following structure:
# - dict: key(tumor_histology) -> RGB_colorI

paletteFile = open(palette, 'r')
colorTumorTypeDict = {}

for line in paletteFile:

    # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        tumorType = line[0]
        rgbColor = line[1]
 
        colorTumorTypeDict[tumorType] = rgbColor

#### Load input tables
rtCountsDf = pd.read_csv(rtCounts, header=0, sep='\t')
nbActiveDf = pd.read_csv(nbActive, header=0, sep='\t')
nbActiveDf.set_index('sampleId', inplace=True)

## Add the number of active source elements to the counts dataframe
rtCountsDf.insert(0, 'nbActiveSrc', nbActiveDf['nbActiveSrc'])


#### Make scatterplot
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(1, 1, 1)

for tumorType in sorted(colorTumorTypeDict.keys()):
    
    dataframe = rtCountsDf[rtCountsDf["tumorType"] == tumorType]

    color = colorTumorTypeDict[tumorType]    

    nbL1List = dataframe["nbL1"].values
    nbActiveList = dataframe["nbActiveSrc"].values
 
    plt.scatter(nbL1List, nbActiveList, c=color, edgecolor='black', linewidth='0.5', s=30, alpha=1, label=tumorType)

## Add axis labels
plt.xlabel('Number somatic L1 insertions', fontsize=12)
plt.ylabel('Number active L1 source elements', fontsize=12)

## Add legend
plt.legend(ncol=2, loc=9, bbox_to_anchor=(0.5, -0.1))


## Compute and add correlation to the plot
corr = stats.spearmanr(rtCountsDf["nbL1"], rtCountsDf["nbActiveSrc"])
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'spear_rho: ' + str(coefficient) + '; P = ' + str(pvalue)
ax.text(0.5, 0.1, text, transform = ax.transAxes)


## Save figure
outPath = outDir + '/nbL1_nbActiveSrc_corr.pdf'
plt.savefig(outPath)
outPath = outDir + '/nbL1_nbActiveSrc_corr.svg'
plt.savefig(outPath)



####
header("Finished")

