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
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")


## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('transductionsCount', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
transductionsCount = args.transductionsCount
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "transductionsCount: ", transductionsCount
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### Read input 
transductionsCountDf = pd.read_csv(transductionsCount, header=0, index_col=0, sep='\t')

transductionsCountDf = transductionsCountDf.T

#### Select hot source L1
hotSourceList = ["22q12.1", "Xp22.2-1", "Xp22.2-2", "3q21.1", "6p22.1", "2q24.1", "7p12.3", "3q26.1"]
transductionsCountHotDf = transductionsCountDf[hotSourceList]

#### Generate dataframe with the number of transductions per donor and source element
tupleList = []

## For each L1 source element
for sourceL1 in transductionsCountHotDf:

    ## For each donor
    for donorId, nbTd in transductionsCountHotDf[sourceL1].iteritems():
        if (nbTd > 0):        
            nbTdTuple = (sourceL1, nbTd)
            tupleList.append(nbTdTuple)                       

## Convert into dataframe
hotL1Df =  pd.DataFrame(tupleList)
hotL1Df.columns = ['cytobandId', 'nbTransductions']

#### Plot data as swarmplot
fig = plt.figure(figsize=(3,5))

ax = sns.swarmplot(x='cytobandId', y='nbTransductions', data=hotL1Df, size=3, edgecolor="gray", order=hotSourceList)

### Axis labels
ax.set_xlabel('')
ax.set_ylabel('Number of transductions')

# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(90)

## Y ticks
ax.set(yticks=np.arange(0,91,10))

## Save figure 
fileName = outDir + "/strombolian_plinnian_examples.pdf"
plt.savefig(fileName)

####
header("Finished")
