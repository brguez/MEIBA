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
import operator

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= "Plot the number of active source source elements per tumor genome across each tumor type")
parser.add_argument('activeSource', help='')
parser.add_argument('donorMetadata', help='')
parser.add_argument('palette', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
activeSource = args.activeSource
donorMetadata = args.donorMetadata
palette = args.palette
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "activeSource: ", activeSource
print "donorMetadata: ", donorMetadata
print "palette: ", palette
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Read metadata file
##########################
# Initialize a dictionary with the following structure:
# - dict: key(donorId) -> tumorType

header("1. Read metadata file")

donorMetadataFile = open(donorMetadata, 'r')
donorIdTumorTypeDict = {}

for line in donorMetadataFile:

    # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        donorId = line[1]
        donorExclusion = line[3]

        histologyCount	= line[10]
        histologyExclusion = line[11]
        tumorHistology = line[12]

        ## Discard excluded donors for tumor types analysis (initial: 2813, after_excluding: 2743). 4 possible exclusion reasons:
        # - TraFiC excluded (22)
        # - Unknown histology (30)
        # - More than 1 possible histology (10)
        # - Excluded histology cohort (32)
        if (donorExclusion == 'Whitelist') and (tumorHistology != "UNK") and (histologyCount == "1") and (histologyExclusion == "included") :

            donorIdTumorTypeDict[donorId] = tumorHistology

# Convert into dataframe
donorIdTumorTypeSeries = pd.Series(donorIdTumorTypeDict, name='tumorType') 

#### 2. Read palette file
##########################
# Initialize a dictionary with the following structure:
# - dict: key(tumor_histology) -> RGB_colorI

header("2. Read palette file")

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

print colorTumorTypeDict

#### 3. Load number of active source elements into a dataframe
################################################################

header("3. Load number of active source elements per donor into a dataframe")

activeSourceSeries = pd.read_csv(activeSource, header=0, index_col=0, sep='\t')

## Filter out excluded donors:
donorIdList = donorIdTumorTypeSeries.keys()
activeSourceFilteredSeries = activeSourceSeries.loc[donorIdList]


#### 4. Merge series generated in 0) and 1) into a single dataframe
#####################################################################

header("4. Merge series generated in 0) and 1) into a single dataframe")
tumorTypeActiveSourceDf = pd.concat([donorIdTumorTypeSeries, activeSourceFilteredSeries], axis=1)

outFilePath = outDir + '/nbActive_germlineSource_tumourTypes.tsv'
tumorTypeActiveSourceDf.to_csv(outFilePath, sep='\t') 

## Select donors with at least one active source element in a single donor
tumorTypesDict = Counter(tumorTypeActiveSourceDf[tumorTypeActiveSourceDf["nbActiveSrc"] > 0]["tumorType"].tolist())

## Make list of tumor types with at least one active source element in a single donor 
selectedTumorTypesList =  [w for w in sorted(tumorTypesDict, key=tumorTypesDict.get, reverse=True)]

## Select donors of selected tumor types
tumorTypeActiveSourceFilteredDf = tumorTypeActiveSourceDf[tumorTypeActiveSourceDf["tumorType"].isin(selectedTumorTypesList)]


#### 5. Compute the average number of germline active source elements per tumor type
######################################################################################

header("5. Compute the average number of germline active source elements per tumor type")

tumorTypeNbActiveSrcDict = {k: g["nbActiveSrc"].tolist() for k,g in tumorTypeActiveSourceFilteredDf.groupby("tumorType")}

meanDict = {}
for tumorType in tumorTypeNbActiveSrcDict:

    meanDict[tumorType] = np.mean(tumorTypeNbActiveSrcDict[tumorType])


## Sort the tumour types according to their Z-score values
sortedTumorTypesList = [ "Eso-AdenoCA", "Lung-SCC", "Head-SCC", "ColoRect-AdenoCA", "Stomach-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Lung-AdenoCA", "Biliary-AdenoCA", "Panc-Endocrine", "Bone-Epith", "Lymph-BNHL", "Ovary-AdenoCA", "Bladder-TCC", "Bone-Osteosarc", "Kidney-ChRCC", "Prost-AdenoCA", "Panc-AdenoCA", "Breast-AdenoCA", "Skin-Melanoma", "Kidney-RCC", "Liver-HCC" ]

#### 6. Make the strip plot
############################

header("6. Make the strip plot")

fig = plt.figure(figsize=(12,4))

ax = sns.stripplot(x='tumorType', y='nbActiveSrc', data=tumorTypeActiveSourceFilteredDf, size=6, edgecolor="black", linewidth=1, alpha=.8, jitter=True, palette=colorTumorTypeDict, order=sortedTumorTypesList)

### Axis labels
ax.set_xlabel('')
ax.set_ylabel('Active source elements per sample')

# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(90)

## Y ticks
ax.set(yticks=np.arange(0,26,2))

## Save figure 
fileName = outDir + "/Pictures/nbActiveSourceElementsPerDonor_striplot.pdf"
plt.savefig(fileName)


#### End
header("FINISH!!")


