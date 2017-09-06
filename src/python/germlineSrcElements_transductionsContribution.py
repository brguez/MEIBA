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
# - dict: key(donorId) -> tumorType

header("1. Read metadata file")

donorMetadataFile = open(donorMetadata, 'r')
donorIdTumorTypeDict = {}
donorIdAncestryDict = {}

for line in donorMetadataFile:

    # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        donorId = line[1]
        #donorId = line[0]
        donorExclusion = line[3]
        ancestry = line[4]
        histologyCount	= line[10]
        histologyExclusion = line[11]
        tumorHistology = line[12]

        ## Discard excluded donors for tumor types analysis (initial: 2813, after_excluding: 2704). 4 possible exclusion reasons:
        # - TraFiC excluded (60)
        # - Unknown histology (30)
        # - More than 1 possible histology (10)
        # - Excluded histology cohort (32)
        if (donorExclusion == 'Whitelist') and (tumorHistology != "UNK") and (histologyCount == "1") and (histologyExclusion == "included") :

            donorIdTumorTypeDict[donorId] = tumorHistology
            donorIdAncestryDict[donorId] = ancestry


#### 1. Compute the total number of transductions per source element in each tumor type
#########################################################################################
header("1. Compute the total number of transductions per source element in each tumor type")

transductionsCountDf = pd.read_csv(transductionsCountFile, header=0, index_col=0, sep='\t')
donorIdList = donorIdTumorTypeDict.keys()

transductionsCountFilteredDf = transductionsCountDf.loc[:, donorIdList]

## Generate dictionary with the following format:
# - dict: key(sourceElementId) -> dict2: key(tumorType) -> total_number_transductions
nbTdPerTumorTypeDict = {}

for sourceElementId, tdPerDonorSeries in transductionsCountFilteredDf.iterrows():

    nbTdPerTumorTypeDict[sourceElementId] = {}
        
    for donorId, nbTd in tdPerDonorSeries.iteritems():

        tumorType = donorIdTumorTypeDict[donorId]
        
        if tumorType not in nbTdPerTumorTypeDict[sourceElementId]:
            nbTdPerTumorTypeDict[sourceElementId][tumorType] = nbTd    

        else:
            nbTdPerTumorTypeDict[sourceElementId][tumorType] += nbTd


## Create pandas dataframe from dictionary
nbTdPerTumorTypeDf = pd.DataFrame(nbTdPerTumorTypeDict) 

## transpose dictionary to have source elements as rows and tumor types as columns
nbTdPerTumorTypeDf = nbTdPerTumorTypeDf.T 

# Save output into tsv
outFilePath = outDir + '/germline_srcElements_tdPerTumorType.tsv'
nbTdPerTumorTypeDf.to_csv(outFilePath, sep='\t') 


#### 2. Compute the total number of transductions per source element in each ancestry
#########################################################################################
header("2. Compute the total number of transductions per source element in each ancestry")


## Generate dictionary with the following format:
# - dict: key(sourceElementId) -> dict2: key(tumorType) -> total_number_transductions
nbTdPerAncestryDict = {}

for sourceElementId, tdPerDonorSeries in transductionsCountFilteredDf.iterrows():

    nbTdPerAncestryDict[sourceElementId] = {}
        
    for donorId, nbTd in tdPerDonorSeries.iteritems():

        ancestry = donorIdAncestryDict[donorId]
        #print "TIOO: ", donorId, ancestry, nbTd
        
        if ancestry not in nbTdPerAncestryDict[sourceElementId]:
            nbTdPerAncestryDict[sourceElementId][ancestry] = nbTd    

        else:
            nbTdPerAncestryDict[sourceElementId][ancestry] += nbTd

#print "nbTdPerAncestryDict: ", nbTdPerAncestryDict
## Create pandas dataframe from dictionary
nbTdPerAncestryDf = pd.DataFrame(nbTdPerAncestryDict) 

## transpose dictionary to have source elements as rows and ancestries as columns
nbTdPerAncestryDf = nbTdPerAncestryDf.T 

# Save output into tsv
outFilePath = outDir + '/germline_srcElements_tdPerAncestry.tsv'
nbTdPerAncestryDf.to_csv(outFilePath, sep='\t') 


#### 3. Compute the contribution (%) of each source element to the total number of
#################################################################################
# transductions in each tumor type
####################################
### Notes:
## contribution(%) = nb_transductions_src_element/total_nb_transductions_cancer_type

header("3. Compute the contribution (%) of each source element to the total number of transductions in each tumor type")

#### 3.1 Select those tumor types with at least 1 transduction
## Compute the total number of transductions per tumor type
totalNbTdPerTumorTypeSerie = nbTdPerTumorTypeDf.sum(axis=0)

## Filter serie generated in the previous step to select only those tumor types with at least 1 transduction
totalNbTdPerTumorTypeFilteredSerie = totalNbTdPerTumorTypeSerie[totalNbTdPerTumorTypeSerie >= 1]

## Filter dataframe generated in 1) to select only those tumor types with at least 1 transduction
tumorTypesList = list(totalNbTdPerTumorTypeFilteredSerie.index)
nbTdPerTumorTypeFilteredDf = nbTdPerTumorTypeDf[tumorTypesList] 

#### 3.2 Compute the source element contribution (in %) for the selected tumor types to the total number of transductions
srcElementIds = nbTdPerTumorTypeFilteredDf.index
tumorTypes = nbTdPerTumorTypeFilteredDf.columns

srcElementContributionDf = pd.DataFrame(index=srcElementIds, columns=tumorTypes)

# Iterate over row index labels (source element ids)
for srcElementId in srcElementIds:

    # Iterate over column index labels (tumor types)
    for tumorType in tumorTypes:
    
        nbTdSrcElement = nbTdPerTumorTypeFilteredDf.loc[srcElementId, tumorType]
        nbTdTumorType = totalNbTdPerTumorTypeFilteredSerie.loc[tumorType]
        tdPercTumorType = float(nbTdSrcElement)/float(nbTdTumorType) * 100

        ## Add source element contribution to dataframe
        srcElementContributionDf.loc[srcElementId, tumorType] = tdPercTumorType

#print "test: ", srcElementContributionDf.sum(axis=0)
# Ok, all 100% as expected. 


#### 3. Compute the contribution (%) of each source element to the total number of
#################################################################################
# transductions in each ancestry
####################################
### Notes:
## contribution(%) = nb_transductions_src_element/total_nb_transductions_ancestry

header("3. Compute the contribution (%) of each source element to the total number of transductions in each ancestry")

#### 3.1 Select those tumor types with at least 1 transduction
## Compute the total number of transductions per tumor type
totalNbTdPerAncestrySerie = nbTdPerAncestryDf.sum(axis=0)

print "totalNbTdPerAncestrySerie: ", totalNbTdPerAncestrySerie

#### 3.2 Compute the source element contribution (in %) for the selected tumor types to the total number of transductions
srcElementIds = nbTdPerAncestryDf.index
ancestries = nbTdPerAncestryDf.columns

srcElementContributionAncestryDf = pd.DataFrame(index=srcElementIds, columns=ancestries)

# Iterate over row index labels (source element ids)
for srcElementId in srcElementIds:

    # Iterate over column index labels (tumor types)
    for ancestry in ancestries:
    
        nbTdSrcElement = nbTdPerAncestryDf.loc[srcElementId, ancestry]
        nbTdAncestry = totalNbTdPerAncestrySerie.loc[ancestry]
        tdPercAncestry = float(nbTdSrcElement)/float(nbTdAncestry) * 100

        ## Add source element contribution to dataframe
        srcElementContributionAncestryDf.loc[srcElementId, ancestry] = tdPercAncestry

# print "test: ", srcElementContributionAncestryDf.sum(axis=0)
# Ok, all 100% as expected. 

#### 3. Compute the contribution (%) of each source element to 
##############################################################
# the total number of transductions in PCAWG
#############################################

header("3. Compute the contribution (%) of each source element to the total number of transductions in PCAWG")

#### 3.1 Incorporate the total number of transductions in PCAWG into the serie with the number of transductions per tumor type
## Compute the total number of transductions in the PCAWG cohort
totalNbTdPCAWG = transductionsCountDf.values.sum()

## Add the total number of transductions in PCAWG to the tumor type serie
totalNbTdPerTumorTypeFilteredSerie["PCAWG"] = totalNbTdPCAWG

## Sort serie in descending order:
totalNbTdPerTumorTypeSortedSerie = totalNbTdPerTumorTypeFilteredSerie.sort_values(ascending=False)

#### 3.2 Incorporate the PCAWG source element contribution to the tumor type contribution dataframe
## Compute the total number of transductions per source element in PCAWG cohort
nbTdPCAWGSerie = transductionsCountDf.sum(axis=1)

## Compute source element contribution to the total number of transductions in PCAWG
contribution = lambda x: (float(x)/totalNbTdPCAWG)*100

srcElementContributionPCAWGSerie = nbTdPCAWGSerie.apply(contribution)

## Add total PCAWG contribution to the dataframe with source element contribution across tumor types
srcElementContributionDf["PCAWG"] = srcElementContributionPCAWGSerie

#### 3.3 Generate tsv with the total number of transductions per source element in PCAWG
## Sort source elements
nbTdPCAWGSortedSerie = nbTdPCAWGSerie.sort_values(ascending=False) # order the source elements in decreasing order. 

## Generate tsv
outFilePath = outDir + '/germline_srcElements_tdPCAWG.tsv'
nbTdPCAWGSortedSerie.to_csv(outFilePath, sep='\t')


#### 4. Sort and filter source element contribution dataframes
##############################################################

header("4. Sort and filter source element contribution dataframes")

#### 4.2 Filter dataframes
### A) Group source elements contributing less than 1% into the category Other
## Filter out those source elements contributing less than 1% to the total number of transductions in PCAWG
srcElementContributionFilteredDf = srcElementContributionDf.loc[srcElementContributionDf['PCAWG'] >= 1]

## Group those elements contributing less than 1% into the category 'Other'
srcElementOtherContributionDf = srcElementContributionDf.loc[srcElementContributionDf['PCAWG'] < 1]

nbSrcElementsOther = len(srcElementOtherContributionDf.index)

otherContributionSerie = srcElementOtherContributionDf.sum(axis=0)

## Make dataframe with the other category as an additional row
rowName = "Other(" + str(nbSrcElementsOther)  + ")" 

srcElementContributionFilteredDf.loc[rowName] = otherContributionSerie

### B) Group tumor types with less than 100 transductions into the category Other
## Filter out those tumor types with less that 100 transductions 
nbTdPerTumorTypeFilteredSerie = totalNbTdPerTumorTypeSortedSerie[totalNbTdPerTumorTypeSortedSerie >= 100]
tumorTypeList = nbTdPerTumorTypeFilteredSerie.index.values.tolist()
srcElementContributionFilteredDf = srcElementContributionFilteredDf[tumorTypeList]

## Group those tumor types with less than 100 transductions into the category 'Other'
tumorTypeList.pop(0)
nbTdPerTumorTypeFilteredDf = nbTdPerTumorTypeDf.drop(tumorTypeList, axis=1)

## Total number of transductions in the 'Other' category
nbTdOther = nbTdPerTumorTypeFilteredDf.values.sum()

## Total number of transductions per source element in 'Other'
nbTdOtherSerie = nbTdPerTumorTypeFilteredDf.sum(axis=1)

## source element contribution to the total number of transductions in 'Other'
contribution = lambda x: (float(x)/nbTdOther)*100

otherContributionSerie = nbTdOtherSerie.apply(contribution)

## Add the 'Other' category 
nbTumorTypesOther = len(tumorTypeList)
colName = "Other(" + str(nbTumorTypesOther) + ")" 

srcElementContributionFilteredDf[colName] = otherContributionSerie
nbTdPerTumorTypeFilteredSerie[colName] = nbTdOther

#### 4.2 Sort dataframes
## Sort source elements (rows) in descending order of transductions contribution to PCAWG 
srcElementContributionSortedDf = srcElementContributionDf.sort_values('PCAWG', ascending=False)
srcElementContributionFilteredSortedDf = srcElementContributionFilteredDf.sort_values('PCAWG', ascending=False)

## Sort tumor types (columns) in descending order of total number of transductions
colOrder = totalNbTdPerTumorTypeSortedSerie.index.tolist()
srcElementContributionSortedDf = srcElementContributionSortedDf[colOrder]

colOrder = nbTdPerTumorTypeFilteredSerie.index.tolist()
srcElementContributionFilteredSortedDf = srcElementContributionFilteredSortedDf[colOrder]

## Sort source element in descending order of total number of transductions
nbTdPCAWGSortedSerie = nbTdPCAWGSerie.sort_values(ascending=False)

## Select those source elements contributing more than 1% from ancestries dataframe:
sourceElementList = srcElementContributionFilteredSortedDf.index.tolist()[1:]

srcElementContributionAncestryFilteredDf = srcElementContributionAncestryDf.loc[sourceElementList]
srcElementContributionAncestryFilteredDf["PCAWG"] = srcElementContributionFilteredSortedDf["PCAWG"]

nbTdPerAncestryFilteredDf = nbTdPerAncestryDf.loc[sourceElementList]
nbTdPerAncestryFilteredDf["PCAWG"] = nbTdPCAWGSortedSerie

#### 4.3 Save both unfiltered and filtered dataframes into a tsv

## Unfiltered dataframe
outFilePath = outDir + '/germline_srcElements_contribution.tsv'
srcElementContributionSortedDf.to_csv(outFilePath, sep='\t') 

## Filtered dataframe
outFilePath = outDir + '/germline_srcElements_contribution.filtered.tsv'
srcElementContributionFilteredSortedDf.to_csv(outFilePath, sep='\t') 


#### 5. Make plots
#######################

header("5. Make plots")

## Set plotting style
sns.set_style("whitegrid")

#### 5.1 Prepare Data

## transpose dataframes to have source elements as columns and tumor types as rows
srcElementContributionSortedDf = srcElementContributionSortedDf.T
srcElementContributionFilteredSortedDf = srcElementContributionFilteredSortedDf.T 
 
## Convert dataframe values into floats 
srcElementContributionSortedDf = srcElementContributionSortedDf.apply(lambda x: pd.to_numeric(x))
srcElementContributionFilteredSortedDf = srcElementContributionFilteredSortedDf.apply(lambda x: pd.to_numeric(x))

#### 5.2 Make heatmaps
### ** Not filtered dataframe **
## A) Source element contribution across tumor types
fig = plt.figure(figsize=(60,20))
fig.suptitle('')
ax1 = sns.heatmap(srcElementContributionSortedDf, vmin=0, vmax=50, annot=True, fmt=".1f", linewidths=.5, cmap=plt.cm.Oranges, cbar=True, annot_kws={"size": 9}, square=True)

ax1.set_xlabel('')
ax1.set_ylabel('')
ax1.xaxis.tick_top()

# turn the axis labels
for item in ax1.get_yticklabels():
    item.set_rotation(0)

for item in ax1.get_xticklabels():
    item.set_rotation(45)

## Save figure 
fileName = outDir + "/Pictures/germline_srcElements_contribution_tumorTypes_heatmap.pdf"
plt.savefig(fileName)


## B) Heatmap with total number of transductions per tumor type
fig = plt.figure(figsize=(5,60))
fig.suptitle('')

totalNbTdPerTumorTypeSortedDf = totalNbTdPerTumorTypeSortedSerie.to_frame(name="# transductions")
#totalNbTdPerTumorTypeSortedDf = totalNbTdPerTumorTypeSortedDf.T

ax3 = sns.heatmap(totalNbTdPerTumorTypeSortedDf, vmin=0, vmax=500, annot=True, fmt='.0f', linewidths=.5, cmap=plt.cm.Blues, cbar=True, annot_kws={"size": 9}, square=True)

ax3.set_xlabel('')
ax3.set_ylabel('')
ax3.xaxis.tick_top()

# turn the axis labels
for item in ax3.get_yticklabels():
    item.set_rotation(0)

for item in ax3.get_xticklabels():
    item.set_rotation(45)

## Save figure 
fileName = outDir + "/Pictures/totalNbTd_per_tumortype_heatmap.pdf"
plt.savefig(fileName)


## C) Heatmap with total number of transductions per source element
fig = plt.figure(figsize=(60,5))
fig.suptitle('')

nbTdPCAWGSortedDf = nbTdPCAWGSortedSerie.to_frame(name="# transductions")
nbTdPCAWGSortedDf = nbTdPCAWGSortedDf.T

ax4 = sns.heatmap(nbTdPCAWGSortedDf, vmin=0, vmax=500, annot=True, fmt='.0f', linewidths=.5, cmap=plt.cm.Greens, cbar=True, annot_kws={"size": 9}, square=True)

ax4.set_xlabel('')
ax4.set_ylabel('')
ax4.xaxis.tick_top()

# turn the axis labels
for item in ax4.get_yticklabels():
    item.set_rotation(0)

for item in ax4.get_xticklabels():
    item.set_rotation(45)

## Save figure 
fileName = outDir + "/Pictures/totalNbTd_per_sourceElement_heatmap.pdf"
plt.savefig(fileName)


### ** Filtered dataframe ** 
## A) Source element contribution across tumor types
fig = plt.figure(figsize=(11,5))
fig.suptitle('')
ax1 = sns.heatmap(srcElementContributionFilteredSortedDf, vmin=0, vmax=50, annot=True, fmt=".1f", linewidths=.5, cmap=plt.cm.Oranges, cbar=True, annot_kws={"size": 9}, square=True)

ax1.set_xlabel('')
ax1.set_ylabel('')
ax1.xaxis.tick_top()

# turn the axis labels
for item in ax1.get_yticklabels():
    item.set_rotation(0)

for item in ax1.get_xticklabels():
    item.set_rotation(45)

## Save figure 
fileName = outDir + "/Pictures/germline_srcElements_contribution_tumorTypes_heatmap.filtered.pdf"
plt.savefig(fileName)


## B) Heatmap with the total number of transductions per tumor type
fig = plt.figure(figsize=(5,11))
fig.suptitle('')

nbTdPerTumorTypeFilteredDf = nbTdPerTumorTypeFilteredSerie.to_frame(name="# transductions")

ax3 = sns.heatmap(nbTdPerTumorTypeFilteredDf, vmin=0, vmax=500, annot=True, fmt='.0f', linewidths=.5, cmap=plt.cm.Blues, cbar=True, annot_kws={"size": 9}, square=True)

ax3.set_xlabel('')
ax3.set_ylabel('')
ax3.xaxis.tick_top()

# turn the axis labels
for item in ax3.get_yticklabels():
    item.set_rotation(0)

for item in ax3.get_xticklabels():
    item.set_rotation(45)

## Save figure 
fileName = outDir + "/Pictures/totalNbTd_per_tumortype_heatmap.filtered.pdf"
plt.savefig(fileName)


#### 5.3 Make barplots

print "srcElementContributionAncestryFilteredDf: ", srcElementContributionAncestryFilteredDf

print "nbTdPerAncestryFilteredDf: ", nbTdPerAncestryFilteredDf

for ancestry, contributionSeries in srcElementContributionAncestryFilteredDf.iteritems():

    # Initialize the matplotlib figure
    fig = plt.figure(figsize=(6, 2))

    # Plot the contribution of each source element to the total number of retrotransposition insertions across PCAWG
    ax = contributionSeries.plot(kind='bar')
    ax.set(ylim=(0, 50), ylabel="Contribution (%)")

    # Add to the top of each bar the total number of retrotransposition events mediated by each src element
    rects = ax.patches
    labels = nbTdPerAncestryFilteredDf[ancestry].values

    for rect, label in zip(rects, labels):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2, height, label, ha='center', va='bottom', size='8')

    ## Save figure 
    fileName = outDir + "/Pictures/source_element_contribution_" + ancestry + ".pdf"
    plt.savefig(fileName)

### End
header("FINISH!!")


