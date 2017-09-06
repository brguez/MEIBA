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
parser.add_argument('activity', help='')
parser.add_argument('sortedSrc', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
activity = args.activity
sortedSrc = args.sortedSrc
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "activity: ", activity
print "sortedSrc: ", sortedSrc
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1.  
##########################
## Generate dictionary with the following format:
# - dict: key(sourceElementId) -> list[number_transductions_per_donor]
# Do not consider 0 cases. 

nbTransductionsPerDonorDict = {}

# Make donor Ids list
activityFile = open(activity, 'r')
donorIdList = activityFile.readline().rstrip().split("\t")
donorIdList = donorIdList[1:] # remove first list element

tupleList = []

for line in activityFile:

    nbTransductionsList = line.rstrip().split("\t")
    cytobandId = nbTransductionsList.pop(0)

    ## Initialize dictionary for source element

    ## For each donor source element activity
    for donorIndex, donorId in enumerate(donorIdList):

        nbTransductions = int(nbTransductionsList[donorIndex])

        # Not consider donors in which the source element is not active
        if (nbTransductions != 0):

            Tuple = (cytobandId, nbTransductions)

            tupleList.append(Tuple)                       


## Convert into dataframe
df =  pd.DataFrame(tupleList)
df.columns = ['cytobandId', 'nbTransductions']

#print dataframe

## Make source elements list
sortedSrcFile = open(sortedSrc, 'r')

sourceElementOrder = []

for line in sortedSrcFile:
    colList = line.rstrip().split("\t")
    cytobandId = colList[0]
    
    sourceElementOrder.append(cytobandId)


##### Make plot
#################



fig = plt.figure(figsize=(25,5))
#ax = sns.swarmplot(x='cytobandId', y='nbTransductions', data=hotL1Df, size=3, edgecolor="gray", order=sourceElementOrder)
ax = sns.swarmplot(x='cytobandId', y='nbTransductions', data=df, size=3, edgecolor="gray", order=sourceElementOrder)

### Axis labels
ax.set_xlabel('')
ax.set_ylabel('# transductions')

# turn the axis labels
for item in ax.get_yticklabels():
    item.set_rotation(0)

for item in ax.get_xticklabels():
    item.set_rotation(90)

## Y ticks
ax.set(yticks=np.arange(0,91,10))

## Save figure 
fileName = outDir + "/test.pdf"
plt.savefig(fileName)


####
header("Finished")
