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
import seaborn as sns
import operator

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= "Plot the number of active source source elements per tumor genome across each tumor type")
parser.add_argument('activeSource', help='')
parser.add_argument('histologyOrder', help='File containing histology ordering. One row per histology')
parser.add_argument('palette', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
activeSource = args.activeSource
histologyOrder = args.histologyOrder
palette = args.palette
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "activeSource: ", activeSource
print "histologyOrder: ", histologyOrder
print "palette: ", palette
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

### 1. Read histology and create list with histology ordering
##############################################################
header("1. Read histology and create list with histology ordering")
histologyFile = open(histologyOrder, 'r')
histologyList = []

for line in histologyFile:
    line = line.rstrip('\n')
    line = line.split("\t")
    histology = line[0]
    histologyList.append(histology)

#### 2. Read palette file
##########################
# Initialize a dictionary with the following structure:
# - dict: key(tumor_histology) -> RGB_colorI
header("2. Read palette file")
paletteFile = open(palette, 'r')
colorHistologyDict = {}

for line in paletteFile:

    # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        tumorType = line[0]
        rgbColor = line[1]
 
        colorHistologyDict[tumorType] = rgbColor

#### 3. Load number of active source elements into a dataframe
################################################################
header("3. Load number of active source elements per donor into a dataframe")
activeSourceDf = pd.read_csv(activeSource, header=0, index_col=0, sep='\t')


#### 4. Make the strip plot
############################
header("4. Make the strip plot")
fig = plt.figure(figsize=(10,4))

ax = sns.stripplot(x='tumorHistology', y='nbActiveSrc', data=activeSourceDf, size=6, edgecolor="black", linewidth=0.5, alpha=1, jitter=0.25, palette=colorHistologyDict, order=histologyList)
### Axis labels
ax.set_xlabel('')
ax.set_ylabel('Active source elements / sample')

# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(90)

## Y ticks
ax.set(yticks=np.arange(0,26,2))

## Save figure 
fileName = outDir + "/nbActive_srcElements_perDonor_striplot.pdf"
plt.savefig(fileName)
fileName = outDir + "/nbActive_srcElements_perDonor_striplot.svg"
plt.savefig(fileName)

#### End
header("FINISH!!")


