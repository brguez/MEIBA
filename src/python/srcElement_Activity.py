#!/usr/bin/env python
#coding: utf-8


#### FUNCTIONS ####

def activityStatus(row):

    ### Classify source element according its activity rate in one of these two categories:
    # a) None/Low activity
    if (row['activityRate'] <= 3):
        status = "none/low"        
    
    # b) Moderate/Hot activity
    else:
        status = "moderate/hot" 
    
    return status


def activityStatus2(row):
    """
    Classify source element according its activity rate in one of these categories:

    Activity ranges:
        - none: 0
        - low: (0-3]
        - moderate: (3-9]
        - hot: >9
    """

    # a) None: 0
    if (row['activityRate'] == 0):
        status = "none"        
    
    # b) Low: (0-3] 
    elif (row['activityRate'] > 0) and (row['activityRate'] <= 3):
        status = "low"

    # c) Moderate: (3-9]
    elif (row['activityRate'] > 3) and (row['activityRate'] <= 9):
        status = "moderate"
    
    # d) Hot: >9
    else:
        status = "hot"  
    

     # Moderate/Hot activity, "moderate/hot"
    return status

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
parser.add_argument('alleleCountFile', help='')
parser.add_argument('alleleFreqFile', help='')
parser.add_argument('transductionsCountFile', help='')
parser.add_argument('activityRateFile', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
alleleCountFile = args.alleleCountFile
alleleFreqFile = args.alleleFreqFile
transductionsCountFile = args.transductionsCountFile
activityRateFile = args.activityRateFile
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "alleleCountFile: ", alleleCountFile
print "alleleFreqFile: ", alleleFreqFile
print "transductionsCountFile: ", transductionsCountFile
print "activityRateFile: ", activityRateFile
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1) Make allele count dataframe
#####################################
## Load allele count matrix into dataframe

dfAlleleCount = pd.read_csv(alleleCountFile, header=0, index_col=0, sep='\t')
# print "allele-count-df: ", dfAlleleCount 

#### 2) Make allele frequency across ancestries dataframe
##########################################################
## Load allele count matrix into dataframe

dfAlleleFreqAncestry = pd.read_csv(alleleFreqFile, header=0, index_col=0, sep='\t')
# print "allele-freq-ancestry-df: ", dfAlleleFreqAncestry 

#### 3) Make number of transductions dataframe
###############################################
## Load number of transductions matrix into dataframe

dfTransductionsCount = pd.read_csv(transductionsCountFile, header=0, index_col=0, sep='\t')
# print "activity-df: ", dfTransductionsCount

#### 4) Make activity rate dataframe
#####################################
dfActivityRate = pd.read_csv(activityRateFile, header=0, index_col=0, sep='\t')

#### 5) Compute number of transductions per copy for each source element and tumor type as:
###########################################################################################
## transductions/copy = total_nb_transductions/allele_count
srcElementIds = dfAlleleCount.index
projecCodes = dfAlleleCount.columns

dfTransductionsPerCopy = pd.DataFrame(index=srcElementIds, columns=projecCodes)

# Iterate over row index labels (source element ids)
for srcElementId in srcElementIds:

    # Iterate over column index labels (project codes)
    for projectCode in projecCodes:
    
        transductionsCount = dfTransductionsCount.loc[srcElementId, projectCode]
        alleleCount = dfAlleleCount.loc[srcElementId, projectCode]

        #print "activity-cell", srcElementId, projectCode, transductionsCount        
        #print "alleleCount-cell", srcElementId, projectCode, alleleCount 
        
        ## Compute transductions per copy
        # a) denominator != 0
        if alleleCount != 0:
            transductionsPerCopy = float(transductionsCount)/float(alleleCount) 

        # b) numerator and denominator == 0
        elif (alleleCount == 0) and (transductionsCount == 0):
            transductionsPerCopy = float(0)
                
        # c) numerator != 0 and denominator == 0
        else:
            transductionsPerCopy = float(0)
            print "[ERROR] Inconsistent ALLELE and TRANSDUCTIONS PER COPY values"
            
        ## Add transductions per copy to dataframe
        dfTransductionsPerCopy.loc[srcElementId, projectCode] = transductionsPerCopy

        #print "transductionsPerCopy-cell", srcElementId, projectCode, transductionsPerCopy
        #print "**************"

# print "tranductions-per-copy: ", dfTransductionsPerCopy


#### 6) Plotting 
#################
## Set plotting style
sns.set_style("whitegrid")

##### 6.1) Activity rate histogram
####################################
activityRateList = dfActivityRate['activityRate'].values
fig = plt.figure(figsize=(5,6))
fig.suptitle('Source element activity-rate', fontsize=14)
# print "activityRateList: ", activityRateList 

## Make plot
ax1 = fig.add_subplot(1, 1, 1)
plt.hist(activityRateList, bins=14, color='#008000', alpha=0.75)
plt.xlabel("Transductions/sample", fontsize=14)
plt.ylabel("# L1 source elements", fontsize=12)

## Save figure
fileName = outDir + "/PCAWG_sourceElements_activityRate_hist.pdf"
plt.savefig(fileName)

##### 6.2) none/low vs moderate/hot activity VAF comparison
############################################################

# Prepare Data
###############
## Classify source elements as none/low or moderate/hot activity
dfActivityRate2 = dfActivityRate # temporal

dfActivityRate['activityStatus'] = dfActivityRate.apply(activityStatus, axis=1)


dfActivityRate2['activityStatus'] = dfActivityRate2.apply(activityStatus2, axis=1)
print "activityRate-and-status: ", dfActivityRate2

print "dfActivityRate: ", dfActivityRate

## Write classified elements into an output table:
dfActivityStatus = dfActivityRate["activityStatus"]

# Save output into tsv
outFilePath = outDir + '/germline_source_element_activityGroup.tsv'
dfActivityStatus.to_csv(outFilePath, sep='\t') 

outFilePath = outDir + '/sourceIdNew_alleleCount_VAF_activityRate_activityStatus.tsv'

dfActivityRate2.to_csv(outFilePath, sep='\t') 

## Make statistical test
x = dfActivityRate[(dfActivityRate['activityStatus']=="none/low")]["VAF"].values
y = dfActivityRate[(dfActivityRate['activityStatus']=="moderate/hot")]["VAF"].values

mannWhitneyU = stats.mannwhitneyu(x, y, alternative='two-sided')

Ustat = mannWhitneyU[0]
pvalue = round(mannWhitneyU[1], 4) 

## Make boxplot:
#################
title = "U: " + str(Ustat) + ", p-value: " + str(pvalue)

fig = plt.figure(figsize=(5,6))
fig.suptitle(title, fontsize=14)

ax = sns.boxplot(x="activityStatus", y="VAF",data=dfActivityRate, width=0.5, showfliers=False)
ax = sns.stripplot(x="activityStatus", y="VAF", data=dfActivityRate, jitter=True, color=".3")
ax.set(ylim=(-0.1, 1.1))
ax.set_xlabel('')

# Add mann whitney U statistic and p-value to the plot:

## Save figure 
fileName = outDir + "/low_vs_high_activity_VAF_boxPlot.pdf"
plt.savefig(fileName)

##### 6.3) Plotting heatmap
############################

# Prepare Data
##################
### Convert dataframes values into floats 
dfTransductionsPerCopy = dfTransductionsPerCopy.apply(lambda x: pd.to_numeric(x))
dfAlleleFreqAncestry = dfAlleleFreqAncestry.apply(lambda x: pd.to_numeric(x))

### Filter out those source elements with none/low activity 
# Make list with source elements with moderate/hot activity
moderateHotSrcList = list(dfActivityRate[(dfActivityRate['activityStatus']=="moderate/hot")]["VAF"].index)

# Select only these elements from the dataframes
dfTransductionsPerCopy = dfTransductionsPerCopy.loc[moderateHotSrcList] 
dfAlleleFreqAncestry = dfAlleleFreqAncestry.loc[moderateHotSrcList] 

print "moderateHotSrcList: ", moderateHotSrcList

### Filter out those project codes without any source element with moderate/hot activity
# Make list of project codes with total activity > 1
projectCodeTotalActivitySeries = dfTransductionsPerCopy.sum(axis=0)

# Select only these project codes from the dataframe
projectCodesSrcActivityList = list(projectCodeTotalActivitySeries[projectCodeTotalActivitySeries >= 1].index)

dfTransductionsPerCopyFiltered = dfTransductionsPerCopy[projectCodesSrcActivityList] 

### Reorder dataframes
# Make clustering at both row (source element) and column level
ax = sns.clustermap(dfTransductionsPerCopyFiltered)

rowIndexList = ax.dendrogram_row.reordered_ind
colIndexList = ax.dendrogram_col.reordered_ind

# Reorder dataframes according to clustering ordering
dfTransductionsPerCopyFilteredOrdered = dfTransductionsPerCopyFiltered.iloc[rowIndexList, colIndexList]

# For the allele frequencies heatmap use a custom ordering of ancestries and discard UNK
colIndexList = [ "EUR", "ASN", "SAN", "AMR", "AFR" ]

dfAlleleFreqAncestryFilteredOrdered = dfAlleleFreqAncestry.iloc[rowIndexList, :]
dfAlleleFreqAncestryFilteredOrdered = dfAlleleFreqAncestryFilteredOrdered.loc[: , colIndexList]

### [provisional] Select only those elements with a minimim number of 1 transduction per copy in at least one sample. 
dfTransductionsPerCopyFilteredOrdered = dfTransductionsPerCopyFilteredOrdered.head(n=9)
dfAlleleFreqAncestryFilteredOrdered = dfAlleleFreqAncestryFilteredOrdered.head(n=9)

print "dfTransductionsPerCopyFilteredOrdered: ", dfTransductionsPerCopyFilteredOrdered

print "dfAlleleFreqAncestryFilteredOrdered: ", dfAlleleFreqAncestryFilteredOrdered

### Make heatmap
####################
fig = plt.figure(figsize=(14,6))
fig.suptitle('')

### Activity per tumor type panel
fig.add_subplot(1, 2, 1)
ax1 = sns.heatmap(dfTransductionsPerCopyFilteredOrdered, vmin=0, vmax=5, annot=True, fmt=".2f", linewidths=.5, cmap=plt.cm.Oranges, cbar=False, annot_kws={"size": 8}, )

# cmap=plt.cm.Blues **use this color pallete for source element VAF accross ancestries...

ax1.set_xlabel('')
ax1.set_ylabel('')
ax1.xaxis.tick_top()

# turn the axis labels
for item in ax1.get_yticklabels():
    item.set_rotation(0)

for item in ax1.get_xticklabels():
    item.set_rotation(90)

### VAF per ancestry panel
fig.add_subplot(1, 2, 2)
ax2 = sns.heatmap(dfAlleleFreqAncestryFilteredOrdered, vmin=0, vmax=1, annot=True, fmt=".2f", linewidths=.5, cmap=plt.cm.Blues, cbar=False, annot_kws={"size": 8})

# cmap=plt.cm.Blues **use this color pallete for source element VAF accross ancestries...

ax2.set_xlabel('')
ax2.set_ylabel('')
ax2.xaxis.tick_top()

# turn the axis labels
for item in ax2.get_yticklabels():
    item.set_rotation(0)

for item in ax2.get_xticklabels():
    item.set_rotation(90)

## Save figure 
fileName = outDir + "/transductionsPerCopy_VAF_topSourceElements_heatmap.pdf"
plt.savefig(fileName)

##### 6.4 Plotting strip plot
###############################
fig = plt.figure(figsize=(10,6))
fig.suptitle('# Stripplot', fontsize=14)

ax = sns.stripplot(data=dfTransductionsPerCopy, jitter=True, linewidth=0.25)

for item in ax.get_xticklabels():
    item.set_rotation(90)

ax.set(ylim=(-0.5, 11))

#ax = sns.stripplot(x="day", y="total_bill", data=dfTransductionsPerCopy, jitter=True, figsize=(8,12), linewidth=1)

## Save figure 
fileName = outDir + "/stripPlot.pdf"
plt.savefig(fileName)


##### 7. Save Dataframe into tsv
#################################
outFilePath = outDir + '/germline_source_element_transductionsPerCopy.tsv'
dfTransductionsPerCopy.to_csv(outFilePath, sep='\t') 

####
