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

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('transductionsCountFile', help='')
parser.add_argument('donorMetadata', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
transductionsCountFile = args.transductionsCountFile
donorMetadata = args.donorMetadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "transductionsCountFile: ", transductionsCountFile
print "donorMetadata: ", donorMetadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 0. Read metadata file
##########################
# Initialize a dictionary with the following structure:
# - dict: key(donorId) -> projectCode

header("1. Read metadata file")

donorMetadataFile = open(donorMetadata, 'r')
donorIdProjectCodeDict = {}

for line in donorMetadataFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')

        donorId = line[0]
        tumorType = line[5].split("-")[0]
        donorIdProjectCodeDict[donorId] = tumorType


#### 1. Compute the total number of transductions per source element in each tumor type
#########################################################################################

header("1. Compute the total number of transductions per source element in each tumor type")

nbTransductionsPerTumorTypeDict = {}

## Make donor Ids list
transductionsCountFile = open(transductionsCountFile, 'r')
donorIdList = transductionsCountFile.readline().rstrip().split("\t")
donorIdList = donorIdList[1:] # remove first list element

## Generate dictionary with the following format:
# - dict: key(sourceElementId) -> dict2: key(projectCode) -> total_number_transductions
for line in transductionsCountFile:

    activityList = line.rstrip().split("\t")
    sourceElementId = activityList.pop(0)
    
    ## Initialize dictionary for source element
    nbTransductionsPerTumorTypeDict[sourceElementId] = {}

    # Initialize total number of source element transductions per project code to 0 values:
    for projectCode in donorIdProjectCodeDict.values():
        nbTransductionsPerTumorTypeDict[sourceElementId][projectCode] = 0

    ## For each donor source element activity
    for donorIndex, donorId in enumerate(donorIdList):

        srcElementActivity = int(activityList[donorIndex])

        # a) Project code available
        if donorId in donorIdProjectCodeDict:
            
            projectCode = donorIdProjectCodeDict[donorId]
        
            nbTransductionsPerTumorTypeDict[sourceElementId][projectCode] += srcElementActivity
        
        # b) Project code do not available
        else:
            print "[ERROR] No project code available for current donorId"  


## Create pandas dataframe from dictionary
dfTransductionsCount = pd.DataFrame(nbTransductionsPerTumorTypeDict) 

## transpose dictionary to have source elements as rows and tumor types as columns
dfTransductionsCount=dfTransductionsCount.T 

# Save output into tsv
outFilePath = outDir + '/germline_srcElements_tdPerTumorType.tsv'
dfTransductionsCount.to_csv(outFilePath, sep='\t') 


#### 2. Compute the contribution (%) of each source element to the total number of
#################################################################################
# transductions in each tumor type
####################################
### Notes:
## contribution(%) = nb_transductions_src_element/total_nb_transductions_cancer_type

header("2. Compute the contribution (%) of each source element to the total number of transductions in each tumor type")

#### 2.1 Select those tumor types with at least 1 transduction
## Compute the total number of transductions per tumor type
projectCodeNbTransductionsSerie = dfTransductionsCount.sum(axis=0)

## Filter serie generated in the previous step to select only those tumor types with at least 1 transduction
projectCodeNbTransductionsFilteredSerie = projectCodeNbTransductionsSerie[projectCodeNbTransductionsSerie >= 1]

## Filter dataframe generated in 1) to select only those tumor types with at least 1 transduction
projectCodesList = list(projectCodeNbTransductionsFilteredSerie.index)
dfTransductionsCountFiltered = dfTransductionsCount[projectCodesList] 

#### 2.2 Compute the source element contribution (in %) for the selected tumor types to the total number of transductions
srcElementIds = dfTransductionsCountFiltered.index
projecCodes = dfTransductionsCountFiltered.columns

dfSrcElementContribution = pd.DataFrame(index=srcElementIds, columns=projecCodes)

# Iterate over row index labels (source element ids)
for srcElementId in srcElementIds:

    # Iterate over column index labels (project codes)
    for projectCode in projecCodes:
    
        transductionCountSrcElement = dfTransductionsCountFiltered.loc[srcElementId, projectCode]
        transductionCountTumorType = projectCodeNbTransductionsFilteredSerie.loc[projectCode]
        transductionPercTumorType = float(transductionCountSrcElement)/float(transductionCountTumorType) * 100

        ## Add source element contribution to dataframe
        dfSrcElementContribution.loc[srcElementId, projectCode] = transductionPercTumorType

#print "test: ", dfSrcElementContribution.sum(axis=0)
# Ok, all 100% as expected. 


#### 3. Compute the contribution (%) of each source element to 
##############################################################
# the total number of transductions in PCAWG
#############################################

header("3. Compute the contribution (%) of each source element to the total number of transductions in PCAWG")

#### 3.1 Incorporate the total number of transductions in PCAWG into the serie with the number of transductions per tumor type
## Compute the total number of transductions in the PCAWG cohort
totalNbTransductionsPCAWG = projectCodeNbTransductionsSerie.sum(axis=0)

print "total_nbTd_PCAWG: ", totalNbTransductionsPCAWG

## Add the total number of transductions in PCAWG to the tumor type serie
projectCodeNbTransductionsFilteredSerie["PCAWG"] = totalNbTransductionsPCAWG

## Sort serie in descending order:
projectCodeNbTransductionsSortedSerie = projectCodeNbTransductionsFilteredSerie.sort_values(ascending=False)

#### 3.2 Incorporate the source element contribution to the tumor type dataframe
## Compute the total number of transductions per source element in PCAWG cohort
srcElementNbTransductionsSerie = dfTransductionsCount.sum(axis=1)

## Compute source element contribution to the total number of transductions in PCAWG
contribution = lambda x: (float(x)/totalNbTransductionsPCAWG)*100

srcElementsTotalContributionSerie = srcElementNbTransductionsSerie.apply(contribution)

## Add total PCAWG contribution to the dataframe with source element contribution across tumor types
dfSrcElementContribution["PCAWG"] = srcElementsTotalContributionSerie

#### 3.3 Generate tsv with the total number of transductions per source element in PCAWG
## Sort source elements
srcElementNbTransductionsSerieSorted = srcElementNbTransductionsSerie.sort_values(ascending=False) # order the source elements in decreasing order. 

## Generate tsv
outFilePath = outDir + '/germline_srcElements_tdPCAWG.tsv'
srcElementNbTransductionsSerieSorted.to_csv(outFilePath, sep='\t')


#### 4. Sort and filter source element contribution dataframes
##############################################################

header("4. Sort and filter source element contribution dataframes")

#### 4.1 Sort dataframe
## Sort source elements (rows) in descending order of transductions contribution to PCAWG 
dfSrcElementContributionSorted = dfSrcElementContribution.sort_values('PCAWG', ascending=False)

## Sort project codes (columns) in descending order of total number of transductions
colOrder = projectCodeNbTransductionsSortedSerie.index.tolist()
dfSrcElementContributionSorted = dfSrcElementContributionSorted[colOrder]

#### 4.2 Filter dataframe
## Filter out those source elements contributing less than 1% to the total number of transductions in PCAWG
dfSrcElementContributionFiltered = dfSrcElementContributionSorted.loc[dfSrcElementContributionSorted['PCAWG'] >= 1]

## Group those elements contributing less than 1% in the category 'Other'
dfSrcElementContributionOther = dfSrcElementContributionSorted.loc[dfSrcElementContributionSorted['PCAWG'] < 1]

nbSrcElementsOther = len(dfSrcElementContributionOther.index)

otherContributionSerie = dfSrcElementContributionOther.sum(axis=0)

## Make final dataframe and add the other category as an additional row
rowName = "Other(" + str(nbSrcElementsOther)  + ")" 

dfSrcElementContributionFinal = dfSrcElementContributionFiltered
dfSrcElementContributionFinal.loc[rowName] = otherContributionSerie

#### 4.3 Save both unfiltered and filtered dataframes into a tsv

## Unfiltered dataframe
outFilePath = outDir + '/germline_srcElements_contribution.tsv'
dfSrcElementContributionSorted.to_csv(outFilePath, sep='\t') 

## Filtered dataframe
outFilePath = outDir + '/germline_srcElements_contribution.filtered.tsv'
dfSrcElementContributionFinal.to_csv(outFilePath, sep='\t') 


#### 5. Make plots
#######################

header("5. Make plots")

## Set plotting style
sns.set_style("whitegrid")

#### 5.1 Prepare Data

## transpose dataframe to have source elements as columns and tumor types as rows
dfSrcElementContributionFinal = dfSrcElementContributionFinal.T 

## Convert dataframe values into floats 
dfSrcElementContributionFinal = dfSrcElementContributionFinal.apply(lambda x: pd.to_numeric(x))

#### 5.2 Make heatmaps

### A) Source element contribution across tumor types
#fig = plt.figure(figsize=(17,7))
fig = plt.figure(figsize=(20,10))
fig.suptitle('')
ax1 = sns.heatmap(dfSrcElementContributionFinal, vmin=0, vmax=50, annot=True, fmt=".1f", linewidths=.5, cmap=plt.cm.Oranges, annot_kws={"size": 8}, square=True)

ax1.set_xlabel('')
ax1.set_ylabel('')
ax1.xaxis.tick_top()

# turn the axis labels
for item in ax1.get_yticklabels():
    item.set_rotation(0)

for item in ax1.get_xticklabels():
    item.set_rotation(90)

## Save figure 
fileName = outDir + "/germline_srcElements_contribution_tumorTypes_heatmap.pdf"
plt.savefig(fileName)


### B) Source element contribution in the PCAWG cohort
colList = [ "PCAWG" ] 
totalSrcElementContributionSeries = dfSrcElementContributionFinal.loc[colList, :]

fig = plt.figure(figsize=(8,4))
fig.suptitle('')

ax2 = sns.heatmap(totalSrcElementContributionSeries, vmin=0, vmax=10, annot=True, fmt=".1f", linewidths=.5, cmap=plt.cm.Oranges, annot_kws={"size": 8}, square=True)

ax2.set_xlabel('')
ax2.set_ylabel('')
ax2.xaxis.tick_top()

# turn the axis labels
for item in ax2.get_yticklabels():
    item.set_rotation(0)

for item in ax2.get_xticklabels():
    item.set_rotation(90)

## Save figure 
fileName = outDir + "/germline_srcElements_contribution_heatmap.pdf"
plt.savefig(fileName)

#### End
header("FINISH!!")


