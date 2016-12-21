#!/usr/bin/env python
#coding: utf-8


#### FUNCTIONS ####
def activity2binary(activity):
    """
    """

    # A) Active
    if (activity > 0):
        boolean = 1

    # B) Inactive
    else:
        boolean = 0

    return boolean


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
parser.add_argument('alleleFreqFile', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
transductionsCountFile = args.transductionsCountFile
alleleFreqFile = args.alleleFreqFile
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "transductionsCountFile: ", transductionsCountFile
print "alleleFreqFile: ", alleleFreqFile
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1) Make source element number of transductions per tumor type dataframe
#################################################################################
## Load number of transductions matrix into dataframe
# This dataframe contains the total number of transductions per source elements in each tumor type
dfTransductionsCount = pd.read_csv(transductionsCountFile, header=0, index_col=0, sep='\t')

# Make binary dataframe (source element active in a given cancer type: 1; inactive: 0)
dfTransductionsBinary = dfTransductionsCount.applymap(activity2binary)

#### 2) Make allele frequency across ancestries dataframe
##########################################################
## Load allele count matrix into dataframe
dfAlleleFreqAncestry = pd.read_csv(alleleFreqFile, header=0, index_col=0, sep='\t')

# Convert allele frequencies into percentages (easier to interpret...)
dfAlleleFreqAncestry = dfAlleleFreqAncestry * 100

#### 3) Compute the contribution (%) of each source element to the total number of
####################################################################################
# transductions in each tumor type
####################################
### Notes:
## 1. Tumor types with low retrotransposition activity discarded...
## 2. contribution(%) = nb_transductions_src_element/total_nb_transductions_cancer_type

### Compute the total number of transductions per tumor type
projectCodeNbTransductionsSerie = dfTransductionsCount.sum(axis=0)

### Filter serie generated in the previous step to select only those tumor types with at least 10 transductions
projectCodeNbTransductionsFilteredSerie = projectCodeNbTransductionsSerie[projectCodeNbTransductionsSerie >= 10]

## Filter dataframe generated in 1) to select only those tumor types with at least 10 transduction
projectCodesList = list(projectCodeNbTransductionsFilteredSerie.index)
dfTransductionsCountFiltered = dfTransductionsCount[projectCodesList] 

## Create dataframe with source element contribution (in %) to selected tumor type total number of transductions
srcElementIds = dfTransductionsCountFiltered.index
projecCodes = dfTransductionsCountFiltered.columns

dfSrcElementContribution = pd.DataFrame(index=srcElementIds, columns=projecCodes)

# Iterate over row index labels (source element ids)
for srcElementId in srcElementIds:

    # Iterate over column index labels (project codes)
    for projectCode in projecCodes:
    
        transductionCountSrcElement = dfTransductionsCountFiltered.loc[srcElementId, projectCode]
        transductionCountTumorType = projectCodeNbTransductionsFilteredSerie.loc[projectCode]

        print "transduction-src-cell", srcElementId, projectCode, transductionCountSrcElement        
        print "transduction-tumorType-cell", srcElementId, projectCode, transductionCountTumorType 
        
        transductionPercTumorType = float(transductionCountSrcElement)/float(transductionCountTumorType) * 100

        ## Add source element contribution to dataframe
        dfSrcElementContribution.loc[srcElementId, projectCode] = transductionPercTumorType

        print "transductionPercTumorType-cell", srcElementId, projectCode, transductionPercTumorType
        print "**************"

#print "test1: ", dfSrcElementContribution.sum(axis=0)
# Ok, all 100% as expected. 


#### 4) Compute the contribution (%) of each source element to 
##############################################################
# the total number of transductions in PCAWG
#############################################

### Compute the total number of transductions in PCAWG cohort
totalNbTransductionsPCAWG = projectCodeNbTransductionsSerie.sum(axis=0)

print "total-nb-td-PCAWG: ", totalNbTransductionsPCAWG
### Add the total number of transductions in PCAWG to the serie with total number of transductions per tumor type
projectCodeNbTransductionsFilteredSerie["PCAWG"] = totalNbTransductionsPCAWG

### Sort serie in descending order:
projectCodeNbTransductionsSortedSerie = projectCodeNbTransductionsFilteredSerie.sort_values(ascending=False)

### Compute the total number of transductions per source element in PCAWG cohort
srcElementNbTransductionsSerie = dfTransductionsCount.sum(axis=1)

### Compute source element contribution to the total number of transductions in PCAWG
contribution = lambda x: (float(x)/totalNbTransductionsPCAWG)*100

srcElementsTotalContributionSerie = srcElementNbTransductionsSerie.apply(contribution)

### Add total PCAWG contribution to the dataframe with source element contribution across tumor types:
dfSrcElementContribution["PCAWG"] = srcElementsTotalContributionSerie


#### 5) Sort and filter src element contribution and allele frequency dataframes
#################################################################################

##### 5.1 Source element contribution df 
### Sort source elements (rows) in descending order of transductions contribution to PCAWG 
dfSrcElementContributionSorted = dfSrcElementContribution.sort_values('PCAWG', ascending=False)

### Sort project codes (columns) in descending order of total number of transductions
colOrder = projectCodeNbTransductionsSortedSerie.index.tolist()
dfSrcElementContributionSorted = dfSrcElementContributionSorted[colOrder]

### Filter out those source elements contributing less than 1% to the total number of transductions in PCAWG
dfSrcElementContributionFiltered = dfSrcElementContributionSorted.loc[dfSrcElementContributionSorted['PCAWG'] >= 1]
#print  "src-elements-contributions-filtered: ", dfSrcElementContributionFiltered

### Group those elements contributing less than 1% in the category 'Other'
dfSrcElementContributionOther = dfSrcElementContributionSorted.loc[dfSrcElementContributionSorted['PCAWG'] < 1]

nbSrcElementsOther = len(dfSrcElementContributionOther.index)

otherContributionSerie = dfSrcElementContributionOther.sum(axis=0)

### Make final dataframe and add the other category as an additional row
rowName = "Other(" + str(nbSrcElementsOther)  + ")" 

dfSrcElementContributionFinal = dfSrcElementContributionFiltered
dfSrcElementContributionFinal.loc[rowName] = otherContributionSerie
print "final-contribution: ", dfSrcElementContributionFinal

##### 5.2 Source element allele frequency df 
#### Filter dataframe to only have those src elements selected in the contribution dataframe  
# Also, sort them in the same step to have the elements in the same order...
colOrder = dfSrcElementContributionFinal.index.tolist()
colOrder.pop() # Remove last element as it is 'Other' category...

dfAlleleFreqAncestryFiltered = dfAlleleFreqAncestry.loc[colOrder, :]

#### Discard UNK and use a custom ordering of ancestries 
rowOrder = [ "EUR", "ASN", "SAN", "AMR", "AFR" ]

dfAlleleFreqAncestryFiltered = dfAlleleFreqAncestryFiltered.loc[: , rowOrder]

print "final-ancestry: ", dfAlleleFreqAncestryFiltered

##### 5.3 Source element total number of transductions serie 
### Sort source elements (rows) in descending order of transductions contribution to PCAWG 
# Filter dataframe to only have those src elements selected in the contribution dataframe  
# Also, sort them in the same step to have the elements in the same order...
colOrder = dfSrcElementContributionFiltered.index.tolist()

colOrder.pop() # Remove last element as it is 'Other' category...

srcElementNbTransductionsSerieFiltered = srcElementNbTransductionsSerie.loc[colOrder]

### Group those elements contributing less than 1% in the category 'Other'
otherSrcElementNbTransductionsSerie = srcElementNbTransductionsSerie.drop(colOrder)

totalNbTransductionsOther = otherSrcElementNbTransductionsSerie.sum(axis=0)

### Add the 'Other' category to make a final series 
indexId = "Other(" + str(nbSrcElementsOther)  + ")" 

srcElementNbTransductionsSerieFinal = srcElementNbTransductionsSerieFiltered
srcElementNbTransductionsSerieFinal.loc[indexId] = totalNbTransductionsOther

print "totalNbTransductionsFinal: ", srcElementNbTransductionsSerieFinal

#df[~srcElementNbTransductionsSerie.isin(colOrder)]

# dfSrcElementContributionOther = dfSrcElementContributionSorted.loc[dfSrcElementContributionSorted['PCAWG'] < 1]

# nbSrcElementsOther = len(dfSrcElementContributionOther.index)

# otherContributionSerie = dfSrcElementContributionOther.sum(axis=0)

#### 6) Number of active source elements per tumor type
#######################################################

### Compute the number of active source elements per tumor type
projectCodeNbActiveSrcElementsSerie = dfTransductionsBinary.sum(axis=0)

#print "nb-active-src-elements: ", projectCodeNbActiveSrcElementsSerie

### Filter and order project codes
colOrder = projectCodeNbTransductionsSortedSerie.index.tolist()

colOrder.pop(0) # Remove PCAWG

projectCodeNbActiveSrcElementsSerieFiltered = projectCodeNbActiveSrcElementsSerie[colOrder]


#### 7. Make output matrices:
#############################
### Source element contribution to tumor type transductions total count

# Save output into tsv
outFilePath = outDir + '/germline_source_element_contribution_perTumorType2.tsv'
dfSrcElementContribution.to_csv(outFilePath, sep='\t') 

#### 8) Make plots
#######################

## Set plotting style
sns.set_style("whitegrid")

##### 8.1 Prepare Data

### transpose dataframes to have source elements as columns and tumor types as rows
dfSrcElementContributionFinal = dfSrcElementContributionFinal.T 
dfAlleleFreqAncestryFiltered = dfAlleleFreqAncestryFiltered.T


### Convert dataframe values into floats 
dfSrcElementContributionFinal = dfSrcElementContributionFinal.apply(lambda x: pd.to_numeric(x))
dfAlleleFreqAncestryFiltered =  dfAlleleFreqAncestryFiltered.apply(lambda x: pd.to_numeric(x))
dfAlleleFreqAncestryFiltered =  dfAlleleFreqAncestryFiltered.apply(lambda x: pd.to_numeric(x))

# Save ordered dataframes into tsv
outFilePath = outDir + '/germline_source_element_contribution_perTumorType_ordered2.tsv'
dfSrcElementContributionFinal.to_csv(outFilePath, sep='\t') 

outFilePath = outDir + '/germline_source_element_alleleFreq_perAncestry_ordered2.tsv'
dfAlleleFreqAncestryFiltered.to_csv(outFilePath, sep='\t') 

####### 8.2 Make heatmaps

##### A) Contribution heatmap
fig = plt.figure(figsize=(17,7))
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
fileName = outDir + "/source_element_contribution_tumorTypes_heatmap6.pdf"
plt.savefig(fileName)

##### B) Allele frequency per ancestry heatmap
fig = plt.figure(figsize=(7,4))
fig.suptitle('')

ax2 = sns.heatmap(dfAlleleFreqAncestryFiltered, vmin=0, vmax=100, annot=True, fmt=".1f", linewidths=.5, cmap=plt.cm.Blues, annot_kws={"size": 8}, square=True)

ax2.set_xlabel('')
ax2.set_ylabel('')
ax2.xaxis.tick_top()

# turn the axis labels
for item in ax2.get_yticklabels():
    item.set_rotation(0)

for item in ax2.get_xticklabels():
    item.set_rotation(90)

## Save figure 
fileName = outDir + "/source_element_alleleFreq_ancestry_heatmap6.pdf"
plt.savefig(fileName)

####### 8.2 Make barplots

##### A) Global contribution (%) to PCAWG

# Prepare series for plotting
totalSrcElementContributionSeries = dfSrcElementContributionFinal.loc['PCAWG', :]

print "test: ", totalSrcElementContributionSeries

# Initialize the matplotlib figure
fig = plt.figure(figsize=(12, 6))

# Plot the contribution of each source element to the total number of retrotransposition insertions across PCAWG
ax3 = totalSrcElementContributionSeries.plot(kind='bar')
ax3.set(ylim=(0, 26), ylabel="Somatic retrotransposition contribution (%)")

# turn the axis labels
for item in ax3.get_yticklabels():
    item.set_rotation(0)

for item in ax3.get_xticklabels():
    item.set_rotation(45)

# Add to the top of each bar the total number of retrotransposition events mediated by each src element
rects = ax3.patches
labels = srcElementNbTransductionsSerieFinal.values

for rect, label in zip(rects, labels):
    height = rect.get_height()
    ax3.text(rect.get_x() + rect.get_width()/2, height, label, ha='center', va='bottom', size='10')

## Save figure 
fileName = outDir + "/source_element_global_contribution_barplot.pdf"
plt.savefig(fileName)


##### B) Number of active source elements per tumor type barplot

# Initialize the matplotlib figure
fig = plt.figure(figsize=(6,8))

# Plot number of active source elements per tumor type
ax4 = projectCodeNbActiveSrcElementsSerieFiltered.plot(kind="barh")
ax4.invert_yaxis()

ax4.set(xlim=(0, 75), ylabel="",
       xlabel="# Active Source Elements")

## Save figure 
fileName = outDir + "/source_element_barplot.pdf"
plt.savefig(fileName)



