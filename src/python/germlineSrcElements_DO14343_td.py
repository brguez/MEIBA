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
parser.add_argument('tdGermline', help='')
parser.add_argument('tdSomatic', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
tdGermline = args.tdGermline
tdSomatic = args.tdSomatic
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "tdGermline: ", tdGermline
print "tdSomatic: ", tdSomatic
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. 
#############################################################################################################

header("1. ")

## Germline
tdGermlineDf = pd.read_csv(tdGermline, header=0, index_col=0, sep='\t')

# Select DO14343S donor
DO14343GermlineSerie = tdGermlineDf.loc[: , "d89b1fd6-bef4-4803-8ed3-3b442be600b6"]

# Select those source elements with more than 1 transductions
DO14343GermlineFilteredSerie = DO14343GermlineSerie[DO14343GermlineSerie > 0]

# Convert to dataframe
DO14343GermlineFilteredDf = pd.DataFrame([DO14343GermlineFilteredSerie])
DO14343GermlineFilteredDf = DO14343GermlineFilteredDf.T
DO14343GermlineFilteredDf["type"] = "germline"
DO14343GermlineFilteredDf = DO14343GermlineFilteredDf.rename(columns={'d89b1fd6-bef4-4803-8ed3-3b442be600b6': 'maxActicity'})

# print "DO14343GermlineFilteredDf: ", DO14343GermlineFilteredDf

## Somatic 
tdSomaticDf = pd.read_csv(tdSomatic, header=0, index_col=0, sep='\t')

# Select DO14343S donor
DO14343SomaticFilteredDf =  tdSomaticDf.loc[tdSomaticDf["sampleListStr"] == "d89b1fd6-bef4-4803-8ed3-3b442be600b6"]
DO14343SomaticFilteredSeries = DO14343SomaticFilteredDf["maxActicity"]

# Convert to dataframe
DO14343SomaticFilteredDf = pd.DataFrame([DO14343SomaticFilteredSeries])
DO14343SomaticFilteredDf = DO14343SomaticFilteredDf.T
DO14343SomaticFilteredDf["type"] = "somatic "

# print "DO14343SomaticFilteredDf: ", DO14343SomaticFilteredDf

## Concatenate germline and somatic dataframes into a single one

DO14343FilteredDf = pd.concat([DO14343SomaticFilteredDf, DO14343GermlineFilteredDf])

DO14343FilteredDf['maxActicity'] = DO14343FilteredDf[['maxActicity']].apply(pd.to_numeric)

print "DO14343FilteredDf: ", DO14343FilteredDf

#### 3. Make the strip plot
############################

header("3. Make the strip plot")

fig = plt.figure(figsize=(2,3))

ax = sns.stripplot(x='type', y='maxActicity', data=DO14343FilteredDf, size=2, edgecolor="gray", jitter=True)

### Axis labels
ax.set_xlabel('Source element type')
ax.set_ylabel('# transductions (DO14343)')

# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(90)

## Y ticks
ax.set(yticks=np.arange(0,48,4))

## Save figure 
fileName = outDir + "/DO14343_somatic_germline_sourceElements_striplot.pdf"
plt.savefig(fileName)


#### End
header("FINISH!!")


