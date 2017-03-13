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
parser.add_argument('metadata', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
activity = args.activity
metadata = args.metadata
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "activity: ", activity
print "metadata: ", metadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print



## Start ## 

#### 1. Read source elements equivalences 
##########################################
## Generate dictionary with the following format:
# - dict: key(sourceElementOldId) -> cytobandId
sourceIdEq = {}

metadataFile = open(metadata, 'r')

for line in metadataFile:

    # Skip header
    if not line.startswith('#'):
        fieldsList = line.rstrip().split("\t")
        cytobandId = fieldsList[0]
        sourceIdOld = fieldsList[2]
        sourceIdEq[sourceIdOld] = cytobandId

#### 2.  
##########################
## Generate dictionary with the following format:
# - dict: key(sourceElementId) -> list[number_transductions_per_donor]

# Do not consider 0 cases. 


nbTransductionsPerDonorDict = {}

# Make donor Ids list
activityFile = open(activity, 'r')
donorIdList = activityFile.readline().rstrip().split("\t")
donorIdList = donorIdList[1:] # remove first list element

for line in activityFile:

    nbTransductionsList = line.rstrip().split("\t")
    sourceIdOld = nbTransductionsList.pop(0)
    print "** source element ** ", sourceIdOld

    # Cytoband Id available
    if sourceIdOld in sourceIdEq:
        
        cytobandId = sourceIdEq[sourceIdOld]  

        ## Initialize dictionary for source element
        nbTransductionsPerDonorDict[cytobandId] = []

        ## For each donor source element activity
        for donorIndex, donorId in enumerate(donorIdList):

            nbTransductions = int(nbTransductionsList[donorIndex])

            # Not consider donors in which the source element is not active
            if (nbTransductions != 0):

                nbTransductionsPerDonorDict[cytobandId].append(nbTransductions)
        

#### Make strombolian dataframe
################################
strombolianTupleList = []

for sourceElement in nbTransductionsPerDonorDict:
    
    nbActiveSamples = len(nbTransductionsPerDonorDict[sourceElement])

    if (nbActiveSamples >= 50):

        print "sourceElement.ndDonors: ", sourceElement, len(nbTransductionsPerDonorDict[sourceElement])       
        for nbTransductions in nbTransductionsPerDonorDict[sourceElement]: 
            
           strombolianTuple = (sourceElement, nbTransductions)
           strombolianTupleList.append(strombolianTuple)                       

## Convert into dataframe
strombolianDf =  pd.DataFrame(strombolianTupleList)
strombolianDf.columns = ['cytobandId', 'nbTransductions']
print "strombolianDf: ", strombolianDf
#print dataframe


#### Make plinnean dataframe
################################
plinneanTupleList = []
for sourceElement in nbTransductionsPerDonorDict:
    
    nbActiveSamples = len(nbTransductionsPerDonorDict[sourceElement])

    if (nbActiveSamples > 0):
        maxNbTransductions = max(nbTransductionsPerDonorDict[sourceElement])

        if (maxNbTransductions >= 30) and (nbActiveSamples < 30):

            print "sourceElement.max: ", sourceElement, len(nbTransductionsPerDonorDict[sourceElement]), maxNbTransductions       

            for nbTransductions in nbTransductionsPerDonorDict[sourceElement]: 
                plinneanTuple = (sourceElement, nbTransductions)
                plinneanTupleList.append(plinneanTuple) 
  
## Convert into dataframe             
plinneanDf =  pd.DataFrame(plinneanTupleList)
plinneanDf.columns = ['cytobandId', 'nbTransductions']

print "plinneanDf: ", plinneanDf

#### Combine dataframes
################################

hotL1Df = pd.concat([strombolianDf, plinneanDf])


##### Make plot
#################
sourceElementOrder = ["22q12.1", "Xp22.2-1", "6p22.1", "2q24.1", "7p12.3"]

fig = plt.figure(figsize=(3,5))

#ax = sns.swarmplot(x='cytobandId', y='nbTransductions', data=hotL1Df, size=3, color=".3", edgecolor="gray", order=sourceElementOrder)
ax = sns.swarmplot(x='cytobandId', y='nbTransductions', data=hotL1Df, size=3, edgecolor="gray", order=sourceElementOrder)

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
fileName = outDir + "/strombolian_plinnian_examples.pdf"
plt.savefig(fileName)


####
header("Finished")
