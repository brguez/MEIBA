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
import scipy.stats as stats
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputVCF', help='Multisample VCF with genotyped source elemetns')
parser.add_argument('metadata', help='Text file with the project and ancestry code per donor Id')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
metadata = args.metadata
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "metadata: ", metadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read input multi-sample VCF and generate a VCF object
#############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Read metadata file
##########################
# Initialize three dictionaries with the following structure:
# - dict1a: key(ancestryCode) -> dict1b: key(donorId) -> nbSourceElements
# - dict2a: key(donorId) -> ancestryCode
# - dict3a: key(donorId) -> projectCode
# - dict4a: key(ancestryCode) -> Total number of chromosome copies in the population (number of individuals with a given ancestry * 2 chromosomes)

header("2. Read ancestry codes file")
metadataFile = open(metadata, 'r')

nbSourceElementsDict = {}
donorIdAncestryDict = {}
donorIdProjectCodeDict = {}
ancestryTotalAlleleCountDict = {}

for line in metadataFile:

    line = line.rstrip('\r\n')

    # Skip header
    if not line.startswith("#"):
        line = line.split('\t')

        donorId = line[0]
        tumorType = line[1]
        ancestryCode = line[2]
           
        ### Dict1a and Dict1b
        if ancestryCode not in nbSourceElementsDict:

            # Create dictionary
            nbSourceElementsDict[ancestryCode] = {}

        # Initialize to 0 values:
        nbSourceElementsDict[ancestryCode][donorId] = 0

        ### dict2a
        donorIdAncestryDict[donorId] = ancestryCode

        ### dict3a
        donorIdProjectCodeDict[donorId] = tumorType

        ### dict4a
        # Initialize ancestry counter
        if ancestryCode not in ancestryTotalAlleleCountDict:
            ancestryTotalAlleleCountDict[ancestryCode] = 1
        
        # Update counter
        else:
            ancestryTotalAlleleCountDict[ancestryCode] += 1        

print "nbDonors-ancestry: ", ancestryTotalAlleleCountDict

## Multiply by 2 the number of donors per ancestry to obtain the total allele count:
for ancestry in ancestryTotalAlleleCountDict:

    ancestryTotalAlleleCountDict[ancestry] = ancestryTotalAlleleCountDict[ancestry] * 2

print "total-alleleCount-ancestry: ", ancestryTotalAlleleCountDict


#### 3. Compute parameters 
##########################
# - Number of source elements per donor. Donors classified according to their ancestry 
# - Global allele counts per source element across full PCAWG cohort 
# - Variant allele counts per source element and per tumor type == Source elements load per tumor type
# - Variant allele frequencies per source element and ancestry == Source elements load per ancestry (to do)

## Store info in three dictionaries:
# globalAlleleFreqDict - dict4a: key(sourceElementId) -> allele frequency across full PCAWG cohort
# alleleCountPerTumorDict - dict5a: key(sourceElementId) -> dict5b: key(projectCode) -> allele Count
# alleleCountPerAncestryDict - dict6a: key(sourceElementId) -> dict6b: key(ancestryCode) -> allele Count
# alleleFreqPerAncestryDict - dict7a: key(sourceElementId) -> dict7b: key(ancestryCode) -> allele frequency

globalAlleleFreqDict = {}
alleleCountPerTumorDict = {}
alleleCountPerAncestryDict = {}
alleleFreqPerAncestryDict = {}

header("3. Compute parameters")

# Open output file
outFilePath = outDir + '/germline_source_element_alleleCount_VAF.tsv'
outFile = open(outFilePath, 'a')

# Write header:
row = '#SourceElement' + "\t" + 'alleleCount' + "\t" + 'VAF' + "\n"
outFile.write(row)

## For each MEI:
for MEIObj in VCFObj.lineList:

    end = (MEIObj.infoDict["BKPB"] if "BKPB" in MEIObj.infoDict else "UNK")

    sourceElementId = MEIObj.chrom + ':' + str(MEIObj.pos) + '-' + str(end)
    print "** source element ** ", sourceElementId

    ## Initialize Dict5a and dict7a for source element
    alleleCountPerTumorDict[sourceElementId] = {}
    alleleFreqPerAncestryDict[sourceElementId] = {}

    # Initialize allele count per project code to 0 values:
    for projectCode in donorIdProjectCodeDict.values():
        alleleCountPerTumorDict[sourceElementId][projectCode] = 0        
    
    ## Initialize Dict6a for source element
    alleleCountPerAncestryDict[sourceElementId] = {}

    # Initialize allele count per ancestry to 0 values:
    for ancestry in ancestryTotalAlleleCountDict:
        alleleCountPerAncestryDict[sourceElementId][ancestry] = 0

    #print "Initialized-dict1: ", alleleCountPerTumorDict
    #print "Initialized-dict2: ", alleleCountPerAncestryDict
    
    ## Total number of chromosome copies in the population
    # Number of donors * 2 (diploid, two copies of a given chromosome)
    totalNbChrom = len(MEIObj.genotypesDict) * 2

    #print "total-nb-copies: ", totalNbChrom

    ## For each donor genotype:
    alleleCount = 0
    
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        ancestryCode = donorIdAncestryDict[donorId] 
        projectCode = donorIdProjectCodeDict[donorId]
        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        #### Update counters and store VAF values
        ## A) Insertion absent in reference genome
        if (MEIObj.alt == "<MEI>"):
            
            # print "Insertion absent in reference genome", donorId, genotype, projectCode          
    
            # a) Heterozygous
            if (genotype == "0/1"):
                # print "ancestryCode: ", ancestryCode

                nbSourceElementsDict[ancestryCode][donorId] += 1     
                
                alleleCount += 1
                alleleCountPerTumorDict[sourceElementId][projectCode] += 1
                alleleCountPerAncestryDict[sourceElementId][ancestryCode] += 1
                   
                # print "alleleCount: ", alleleCount
                # print "alleleCount-ProjectCode: ", alleleCountPerTumorDict[sourceElementId][projectCode]
                # print "alleleCount-ancestry: ", alleleCountPerAncestryDict[sourceElementId][ancestryCode]

            # b) Homozygous alternative
            elif (genotype == "1/1"):
                # print "ancestryCode: ", ancestryCode

                nbSourceElementsDict[ancestryCode][donorId] += 1                               
    
                alleleCount += 2
                alleleCountPerTumorDict[sourceElementId][projectCode] += 2
                alleleCountPerAncestryDict[sourceElementId][ancestryCode] += 2

                # print "alleleCount: ", alleleCount
                # print "alleleCount-ProjectCode: ", alleleCountPerTumorDict[sourceElementId][projectCode]
                # print "alleleCount-ancestry: ", alleleCountPerAncestryDict[sourceElementId][ancestryCode]

            # c) possibility would be missing allele (./.)

        ## B) Insertion in reference genome and absent in donor genome
        elif (MEIObj.ref == "<MEI>"):
        
            # print "Insertion in reference genome", donorId, genotype, projectCode        

            # a) Heterozygous
            if (genotype == "0/1"):
                # print "ancestryCode: ", ancestryCode

                nbSourceElementsDict[ancestryCode][donorId] += 1                
           
                alleleCount += 1
                alleleCountPerTumorDict[sourceElementId][projectCode] += 1
                alleleCountPerAncestryDict[sourceElementId][ancestryCode] += 1

                # print "alleleCount: ", alleleCount
                # print "alleleCount-ProjectCode: ", alleleCountPerTumorDict[sourceElementId][projectCode]
                # print "alleleCount-ancestry: ", alleleCountPerAncestryDict[sourceElementId][ancestryCode]

            # b) Homozygous reference
            elif (genotype == "0/0"):
                # print "ancestryCode: ", ancestryCode

                nbSourceElementsDict[ancestryCode][donorId] += 1
                                   
                alleleCount += 2
                alleleCountPerTumorDict[sourceElementId][projectCode] += 2
                alleleCountPerAncestryDict[sourceElementId][ancestryCode] += 2

                # print "alleleCount: ", alleleCount
                # print "alleleCount-ProjectCode: ", alleleCountPerTumorDict[sourceElementId][projectCode]
                # print "alleleCount-ancestry: ", alleleCountPerAncestryDict[sourceElementId][ancestryCode]

            # c) possibility would be missing allele (./.)

        ## C) Raise error...  
        else:
            msg="Incorrectly formated VCF line"
            info(msg)
    
    ## Compute global variant allele frequency (VAF)
    alleleFrequency = float(alleleCount) / totalNbChrom

    # Save VAF a dictionary. One per MEI type
    globalAlleleFreqDict[sourceElementId] = alleleFrequency

    # Write source element allele count and frequency in the output table
    chrom = MEIObj.chrom
    pos = MEIObj.pos
    element = MEIObj.infoDict["CLASS"]
    row = sourceElementId + "\t" + str(alleleCount)  + "\t" + str(alleleFrequency) + "\n"
    outFile.write(row)

    ## Compute source element variant allele frequency across ancestries
    for ancestryCode in alleleCountPerAncestryDict[sourceElementId]:
        totalAlleleCountAncestry = ancestryTotalAlleleCountDict[ancestryCode]
        alleleCountInAncestry = alleleCountPerAncestryDict[sourceElementId][ancestryCode]
        alleleFreqInAncestry = float(alleleCountInAncestry) / totalAlleleCountAncestry
        print "test: ", sourceElementId, ancestryCode, totalAlleleCountAncestry, alleleCountInAncestry, alleleFreqInAncestry
     
        alleleFreqPerAncestryDict[sourceElementId][ancestryCode] = alleleFreqInAncestry
    
    print "---------------------"

print "Finished-dict1: ", alleleCountPerTumorDict
print "Finished-dict2: ", alleleCountPerAncestryDict
print "Finished-dict3: ", alleleFreqPerAncestryDict

#### 4. Make tables:
####################

### 4.1 Source element allele count per tumor type 
# Create pandas dataframe from dictionary
df=pd.DataFrame(alleleCountPerTumorDict) 

# transpose dictionary to have source elements as rows and tumor types as columns
df=df.T 

# Save output into tsv
outFilePath = outDir + '/germline_source_element_alleleCount_perTumorType.tsv'
df.to_csv(outFilePath, sep='\t') 


### 4.2 Source element allele count per tumor type
# Create pandas dataframe from dictionary
df=pd.DataFrame(alleleCountPerTumorDict) 

# transpose dictionary to have source elements as rows and tumor types as columns
df=df.T 

# Save output into tsv
outFilePath = outDir + '/germline_source_element_alleleCount_perTumorType.tsv'
df.to_csv(outFilePath, sep='\t') 


### 4.3 Source element allele frequency per ancestry
# Create pandas dataframe from dictionary
df=pd.DataFrame(alleleFreqPerAncestryDict) 

# transpose dictionary to have source elements as rows and tumor types as columns
df=df.T 

# Save output into tsv
outFilePath = outDir + '/germline_source_element_alleleFreq_perAncestry.tsv'
df.to_csv(outFilePath, sep='\t') 

#### 5. Make plots:
###################
header("5. Make plots")

# - Variant allele frequencies histogram across PCAWG donors 
# - Number of source elements per donor and ancestry. Boxplot 
# - Number of source elements per donor and ancestry. Violin Plot 

#### 5.1 Source element variant allele frequencies across PCAWG donors
header("5.1 Make variant allele frequencies plot")

alleleFreqList = globalAlleleFreqDict.values()

fig = plt.figure(figsize=(5,6))
fig.suptitle('Variant allele frequencies (VAF)', fontsize=14)

## Make plot
ax1 = fig.add_subplot(1, 1, 1)
plt.hist(alleleFreqList, bins=10, color='#008000', alpha=0.75)
plt.xlabel("VAF", fontsize=14)
plt.ylabel("# L1 source elements", fontsize=12)
plt.xlim(0, 1)

# Remove top and right axes
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 1.01, 0.1))
locs, labels = plt.xticks()
# plt.setp(labels, rotation=30)

## Save figure 
fileName = outDir + "/germline_source_element_VAF_hist.pdf"
plt.savefig(fileName)

#### 5.2 Number of source elements per donor and ancestry
header("5.2 Number of source elements per donor and ancestry")

### A) Boxplot
## Organize the data for plotting
tupleListNbSourceElements = []

for ancestryCode in sorted(nbSourceElementsDict):

    # Make tuple (ancestryCode, list with number of source elements per donor)
    nbDonors= len(nbSourceElementsDict[ancestryCode].values())
    xLabel = ancestryCode + '(' +  str(nbDonors) + ')'
    nbSourceElementsTuple = (xLabel, nbSourceElementsDict[ancestryCode].values())

    # Add tuple to the list
    tupleListNbSourceElements.append(nbSourceElementsTuple)


## Make nested list with the following format:
# [donor1_nbSourceElements, donor2_nbSourceElements, ..., donorN_nbSourceElements], [donor1_nbSourceElements, donor2_nbSourceElements, ..., donorN_nbSourceElements] , ... [donor1_nbSourceElements, donor2_nbSourceElements, ..., donorN_nbSourceElements]
#                               ancestry1_list                                                                  ancestry2_list                                                                          ancestryN_list

tmpList = map(list, zip(*tupleListNbSourceElements))
ancestryCodesList = tmpList[0]
nbSourceElementsPerDonor = tmpList[1]

### Plotting
fig = plt.figure(figsize=(5,6))
fig.suptitle('# Source elements per donor', fontsize=14)

ax1 = fig.add_subplot(1, 1, 1)

# Create the boxplot
bp = ax1.boxplot(nbSourceElementsPerDonor)
plt.ylabel("# Source L1", fontsize=12)

## Customize boxplot:
# change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#696969', linewidth=1)

# change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#696969', linewidth=1)

# change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#696969', linewidth=1)

# change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#8b0000', linewidth=2)

# Add the ancestry codes to the x-axis
ax1.set_xticklabels(ancestryCodesList, fontsize = 10)
locs, labels = plt.xticks()
# plt.setp(labels, rotation=25)

# Remove top and right axes
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_axisbelow(True)

## Save figure
fileName = outDir + "/nb_germline_source_element_perDonor_boxplot.pdf"
fig.savefig(fileName)

### B) Violin plot
## Organize the data for plotting into a dictionary:
# - dict1:
#       nbSourceElements -> list[nbSourceElementsDonor1, nbSourceElementsDonor2, ..., nbSourceElementsDonorN]
#       ancestryPerDonor -> list[ancestryDonor1, ancestryDonor2, ..., nbSourceElementsDonorN]

dict4pandas = {}
dict4pandas['nbSourceElements'] = []
dict4pandas['ancestry'] = []

for ancestryCode in sorted(nbSourceElementsDict):

    nbSourceElementsPerDonorList = nbSourceElementsDict[ancestryCode].values()
    nbDonors = len(nbSourceElementsPerDonorList)

    xLabel = ancestryCode + '(' +  str(nbDonors) + ')'
    xLabelList = [ xLabel ] * nbDonors

    dict4pandas['nbSourceElements'] = dict4pandas['nbSourceElements'] + nbSourceElementsPerDonorList
    dict4pandas['ancestry'] = dict4pandas['ancestry'] + xLabelList

# Make pandas dataframe from dict:
dataframe = pd.DataFrame(dict4pandas)

### Plotting

fig = plt.figure(figsize=(5,6))
fig.suptitle('# Source elements per donor', fontsize=14)

# Create the violin plot
ax = sns.violinplot(x='ancestry', y='nbSourceElements', data=dataframe, palette="muted")

# y limit
#sns.plt.ylim(0,21)

## Modify axis labels
ax.set(xlabel='', ylabel='# Source L1')

# Remove top and right axes
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax.set_axisbelow(True)

## Save figure
fileName = outDir + "/nb_germline_source_element_perDonor_violinPlot.pdf"
fig.savefig(fileName)

####
header("Finished")
