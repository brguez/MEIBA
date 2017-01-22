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
import numpy as np
from matplotlib import pyplot as plt

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputVCF', help='multi-sample VCF file containing genotyped MEI')
parser.add_argument('metadata', help='Text file in tabular format containing donor_id/project_code equivalencies')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
metadata =  args.metadata
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

#### 2. Read metadata file
##########################
# Make four dictionaries with the following structure:

# - dict1a: key(projectId) -> dict2a: key(donorId) -> totalNbMEI
# - dict1b: key(projectId) -> dict2a: key(donorId) -> nbHeterozygousMEI
# - dict1c: key(projectId) -> dict2b: key(donorId) -> nbHomozygousMEI
# - dict1d: key(donorId) -> projectCode
# - dict1e: key(projectId) -> dict2a: key(donorId) -> totalNbMEIAlleles

# Note: nbHeterozygousMEI and nbHomozygousMEI are the number of heterozygous and homozygous MEI for a given donor
# Note: nbHeterozygousMEI and nbHomozygousMEI are iniciated with a value of 0
header("2. Read metadata file")
metadataFile = open(metadata, 'r')

totalNbMEIDict = {}
nbHeterozMEIDict = {}
nbHomozMEIDict = {}
donorIdProjectCodes = {}
totalNbMEIAllelesDict = {}

for line in metadataFile:

    line = line.rstrip('\n')
    line = line.split('\t')

    projectCode = line[0]
    donorId = line[1]

    ## Dict1a, Dict1b, Dict1c and Dict1e 
    if projectCode not in nbHeterozMEIDict:

        # Create dictionary
        totalNbMEIDict[projectCode] = {}
        nbHeterozMEIDict[projectCode] = {}
        nbHomozMEIDict[projectCode] = {}
        totalNbMEIAllelesDict[projectCode] = {}

    # Initialize to 0 values:
    totalNbMEIDict[projectCode][donorId] = 0
    nbHeterozMEIDict[projectCode][donorId] = 0
    nbHomozMEIDict[projectCode][donorId] = 0
    totalNbMEIAllelesDict[projectCode][donorId] = 0

    ## Dict1d
    donorIdProjectCodes[donorId] = projectCode

#print "totalNbMEIAllelesDict-empty: ", len(totalNbMEIAllelesDict), len(totalNbMEIDict), totalNbMEIAllelesDict


#### 3. Compute parameters:
###########################
# - Variant allele frequency per MEI across PCAWG cohort
# - Variant allele count per MEI across PCAWG cohort
# - Total number of MEI per donor (Save into dictionary generated in step 2)
# - Total number of MEI alleles per donor (Save into dictionary generated in step 2)
# - Number of heterozygous MEI per donor (Save into dictionary generated in step 2)
# - Number of homozygous MEI per donor (Save into dictionary generated in step 2)
# - Predicted zygosity based on number of supporting reads for heterozygous MEI
# - Predicted zygosity based on number of supporting reads for homozygous MEI

totalAlleleCountList = []
L1AlleleCountList = []
L1AlleleFreqList = []
AluAlleleCountList = []
AluAlleleFreqList = []
SvaAlleleCountList = []
SvaAlleleFreqList = []
ErvkAlleleCountList = []
ErvkAlleleFreqList = []
zygosityHeterozList = []
zygosityHomozList = []

### Make table containing for each germline insertions its VAF
## Two columns:
# insertionsId(chrom:pos)   VAF

# Open output file
outFilePath = outDir + '/MEI_CLASS_VAF.txt'
outFile = open(outFilePath, 'a')

# Write header:
row = '#MEI' + "\t" + 'CLASS' + "\t" +'VAF' + "\n"
outFile.write(row)


header("3. Compute parameters")

for MEIObj in VCFObj.lineList:

    ## Total number of chromosome copies in the population
    # Number of donors * 2 (diploid, two copies of a given chromosome)
    totalNbChrom = len(MEIObj.genotypesDict) * 2

    ## Compute MEI allele count, genotyping VAF and update counters to heterozygous and homozygous MEI per donor
    # MEI allele count: number of chromosomes in the population harvouring the MEI
    # Genotyping VAF: Ratio (number_reads_supporting_MEI)/(total_nb_reads_covering_MEI)
    alleleCount = 0

#    print "MEI: ", MEIObj.chrom, MEIObj.pos

    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        projectCode = donorIdProjectCodes[donorId]
        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]
        nbReadsMEI = float(genotypeFieldList[1])
        totalNbReads = float(genotypeFieldList[2])

        # Compute genotyping VAF:
        if totalNbReads == 0:
            zygosity = 0
        else:
            zygosity =  nbReadsMEI / totalNbReads

        ## Update counters and store VAF values
        # A) Heterozygous
        if (genotype == "0/1"):
            totalNbMEIDict[projectCode][donorId] += 1
            totalNbMEIAllelesDict[projectCode][donorId] += 1
            alleleCount += 1
            nbHeterozMEIDict[projectCode][donorId] += 1
            zygosityHeterozList.append(zygosity)

#            print "test: ", projectCode, donorId, genotype, totalNbMEIAllelesDict[projectCode][donorId]

        # B) Homozygous
        elif (genotype == "1/1"):
            totalNbMEIDict[projectCode][donorId] += 1
            totalNbMEIAllelesDict[projectCode][donorId] += 2
            alleleCount += 2
            nbHomozMEIDict[projectCode][donorId] += 1
            zygosityHomozList.append(zygosity)

#            print "test: ", projectCode, donorId, genotype, totalNbMEIAllelesDict[projectCode][donorId]

    ## Compute MEI allele frequency:
    alleleFrequency = float(alleleCount) / totalNbChrom

    ## Save into list. One per MEI type
    totalAlleleCountList.append(alleleCount)

    if (MEIObj.infoDict['CLASS'] == 'L1'):
        L1AlleleCountList.append(alleleCount)
        L1AlleleFreqList.append(alleleFrequency)

    elif (MEIObj.infoDict['CLASS'] == 'Alu'):
        AluAlleleCountList.append(alleleCount)
        AluAlleleFreqList.append(alleleFrequency)

    elif (MEIObj.infoDict['CLASS'] == 'SVA'):
        SvaAlleleCountList.append(alleleCount)
        SvaAlleleFreqList.append(alleleFrequency)

    else:
        ErvkAlleleCountList.append(alleleCount)
        ErvkAlleleFreqList.append(alleleFrequency)

    ## Add row to the table
    row = MEIObj.chrom + ":" + str(MEIObj.pos) + "\t" + MEIObj.infoDict['CLASS'] + "\t" + str(alleleFrequency) + "\n"
    outFile.write(row)

#    print "***********************"
 
print "totalNbMEIAllelesDict-filled: ", len(totalNbMEIAllelesDict), totalNbMEIAllelesDict

#### 4. Make plots:
#####################
header("4. Make plots")

# - Variant allele frequencies histogram across PCAWG donors (done)
# - Number of MEI per donor and tumor type boxplot (done)
# - Predicted zygosity histogram for homozygous and heterozygous variants (done)
# - Variant allele counts line graph (to do)

#### 4.1 Variant allele frequencies histogram across PCAWG donors
header("4.1 Make variant allele frequencies plot")

fig = plt.figure(figsize=(16,13))
fig.suptitle('Variant allele frequencies (VAF) across PCAWG donors', fontsize=20)

### a) L1
## Make plot
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("LINE-1", fontsize=16)
plt.hist(L1AlleleFreqList, bins=40, color='#008000', alpha=0.75)
plt.xlabel("VAF", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 1)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax1.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 1.01, 0.1))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

### b) ALU
## Make plot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=16)
plt.hist(AluAlleleFreqList, bins=40, color='#008000', alpha=0.75)
plt.xlabel("VAF", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 1)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax2.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 1.01, 0.1))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

### c) SVA
## Make plot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=16)
plt.hist(SvaAlleleFreqList, bins=40, color='#008000', alpha=0.75)
plt.xlabel("VAF", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 1)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax3.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 1.01, 0.1))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

### d) ERVK (no ERVK for now.. we need to improve script for identifying ERVK TSD)
## Make plot
#ax4 = fig.add_subplot(2, 2, 4)
#ax4.set_title("ERVK", fontsize=16)
#plt.hist(ErvkAlleleFreqList, bins=40, color='#008000', alpha=0.75)
#plt.xlabel("VAF", fontsize=14)
#plt.ylabel("# MEI")
#plt.xlim(0, 1)
#plt.ylim(0, max(ErvkAlleleFreqList) + 10)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
#ax4.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)
#ax4.set_axisbelow(True)

## Customize ticks
#plt.xticks(np.arange(0, 1.01, 0.1))
#locs, labels = plt.xticks()
#plt.setp(labels, rotation=30)
#plt.yticks(np.arange(0, max(ErvkAlleleFreqList)))

## Save figure
fileName = outDir + "/PCAWG_cohortVAFs_hist.pdf"
plt.savefig(fileName)

#### 4.2 Make total number and number of heterozygous and homozygous MEI per donor and tumor type boxplots
header("4.2 Make boxplots")

## Organize the data for plotting
tupleListTotalMEI = []
tupleListTotalMEIAlleles = []
tupleListHeterozMEI = []
tupleListHomozMEI = []

for projectCode in sorted(nbHeterozMEIDict):

    # Make tuple (projectCode, number of MEI per donor) for total number of MEI, number of heterozygous and homozygous MEI
    projectCodeNbTotalMEITuple = (projectCode, totalNbMEIDict[projectCode].values())
    projectCodeNbTotalMEIAllelesTuple = (projectCode, totalNbMEIAllelesDict[projectCode].values()) 
    projectCodeNbHeterozMEITuple = (projectCode, nbHeterozMEIDict[projectCode].values())
    projectCodeNbHomozMEITuple = (projectCode, nbHomozMEIDict[projectCode].values())

    # Add tuple to the list
    tupleListTotalMEI.append(projectCodeNbTotalMEITuple)
    tupleListTotalMEIAlleles.append(projectCodeNbTotalMEIAllelesTuple)
    tupleListHeterozMEI.append(projectCodeNbHeterozMEITuple)
    tupleListHomozMEI.append(projectCodeNbHomozMEITuple)

print "tuple-list: ", tupleListTotalMEIAlleles

## Make nested list with the following format:
# [donor1_nbMEI, donor2_nbMEI, ..., donorN_nbMEI], [donor1_nbMEI, donor2_nbMEI, ..., donorN_nbMEI] , ... [donor1_nbMEI, donor2_nbMEI, ..., donorN_nbMEI]
#                  project1_list                                  project2_list                                         projectN_list

# Total MEI
tmpList = map(list, zip(*tupleListTotalMEI))
projectCodesListTotal = tmpList[0]
nbMEIPerDonorTotal = tmpList[1]

# Total MEI alleles
tmpList = map(list, zip(*tupleListTotalMEIAlleles))
projectCodesListTotalAlleles = tmpList[0]
nbMEIPerDonorTotalAlleles = tmpList[1]

# Heterozygous MEI
tmpList = map(list, zip(*tupleListHeterozMEI))
projectCodesListHet = tmpList[0]
nbMEIPerDonorHet = tmpList[1]

# Homozygous MEI
tmpList = map(list, zip(*tupleListHomozMEI))
projectCodesListHom = tmpList[0]
nbMEIPerDonorHom = tmpList[1]

#### A) Make boxplot for total number of MEI
fig = plt.figure(figsize=(14,9))
fig.suptitle('Total # MEI per donor', fontsize=18)

ax1 = fig.add_subplot(1, 1, 1)

# Create the boxplot
bp = ax1.boxplot(nbMEIPerDonorTotal)
plt.ylabel("# MEI", fontsize=18)

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

# Add the project codes to the x-axis
ax1.set_xticklabels(projectCodesListTotal)
locs, labels = plt.xticks()
plt.setp(labels, rotation=75)

# Remove top and right axes
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax1.set_axisbelow(True)

## Save figure
fileName = outDir + "/PCAWG_totalNbMEIperDonor_boxplot.pdf"
fig.savefig(fileName)

#### B) Make boxplot for total number of MEI
fig = plt.figure(figsize=(14,9))
fig.suptitle('# MEI Alleles per donor', fontsize=18)

ax2 = fig.add_subplot(1, 1, 1)

# Create the boxplot
bp = ax2.boxplot(nbMEIPerDonorTotalAlleles)
plt.ylabel("# MEI Alleles", fontsize=18)

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

# Add the project codes to the x-axis
ax2.set_xticklabels(projectCodesListTotalAlleles)
locs, labels = plt.xticks()
plt.setp(labels, rotation=75)

# Remove top and right axes
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax2.set_axisbelow(True)

## Save figure
fileName = outDir + "/PCAWG_totalNbMEIAllelesPerDonor_boxplot.pdf"
fig.savefig(fileName)


#### C) Make boxplots for number of homozygous and heterozygous MEI
fig = plt.figure(figsize=(18,12))
fig.suptitle('# MEI per donor', fontsize=24)

### C.a) Heterozygous MEI
ax3 = fig.add_subplot(2, 1, 1)

# Create the boxplot
bp = ax3.boxplot(nbMEIPerDonorHet)
plt.ylabel("# Heterozygous MEI", fontsize=18)

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

# Add the project codes to the x-axis
ax3.set_xticklabels(projectCodesListHet)
locs, labels = plt.xticks()
plt.setp(labels, rotation=75)

# Remove top and right axes
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax3.set_axisbelow(True)

### C.b) Homozygous MEI
ax4 = fig.add_subplot(2, 1, 2)

# Create the boxplot
bp = ax4.boxplot(nbMEIPerDonorHom)
plt.ylabel("# Homozygous MEI", fontsize=18)

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

# Add the project codes to the x-axis
ax4.set_xticklabels(projectCodesListHom)
locs, labels = plt.xticks()
plt.setp(labels, rotation=75)

# Remove top and right axes
ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_left()

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax4.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax4.set_axisbelow(True)

## Save figure
fileName = outDir + "/PCAWG_nbMEIperDonor_boxplot.pdf"
fig.savefig(fileName)

#### 4.3 Predicted zygosity histogram for homozygous and heterozygous variants
header("4.3 Make zygosity histogram")
fig = plt.figure(figsize=(8,6))
fig.suptitle('Predicted zygosity', fontsize=14)

## Make plot
ax5 = fig.add_subplot(1, 1, 1)
plt.hist(zygosityHeterozList, bins=20, color='#008000', edgecolor='#000000', alpha=0.75, label='Heterozygous')
plt.hist(zygosityHomozList, bins=5, color='#A67D3D', edgecolor='#000000', alpha=0.75, label='Homozygous')
plt.xlabel("Zygosity", fontsize=12)
plt.ylabel("# MEI", fontsize=12)
plt.xlim(0, 1)

# Add an horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax5.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax5.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 1.01, 0.1))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

## Legend
plt.legend(fontsize=12, loc='upper left')

## Save figure
fileName = outDir + "/PCAWG_zygosity_histogram.pdf"
fig.savefig(fileName)

#### 4.4 Variant allele counts line chart
header("4.4 Make allele counts scatterplot")

### Organize the data for plotting

alleleCountNbMEITuple = [(i, totalAlleleCountList.count(i)) for i in set(totalAlleleCountList)]
tmpList = map(list, zip(*alleleCountNbMEITuple))
alleleCountList = tmpList[0]
nbMEIList = tmpList[1]

### Make plot
fig = plt.figure(figsize=(8,8))
fig.suptitle('Allele count spectrum', fontsize=16)
ax6 = fig.add_subplot(1, 1, 1)
plt.scatter(alleleCountList, nbMEIList, color='#008000', alpha=.4)
plt.xlabel("Variant allele count", fontsize=14)
plt.ylabel("# MEI")
ax6.set_xscale('log', basex=10)
ax6.set_yscale('log', basex=10)
plt.xlim(0.5, max(alleleCountList) + 100)
plt.ylim(0.5, max(nbMEIList) + 100)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax6.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax6.set_axisbelow(True)

## Customize ticks
# X axis
xPosList = [ 1, 5, 10, 50, 100, 200, 300, 400, 500, 1000, 2500, max(alleleCountList) ]
ax6.set_xticks(xPosList)
ax6.set_xticklabels(xPosList)
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)

# y axis
yPosList = [ 1, 10, 100, 1000, max(nbMEIList) ]
ax6.set_yticks(yPosList)
ax6.set_yticklabels(yPosList)

## Save figure
fileName = outDir + "/PCAWG_alleleCount_scatterplot.pdf"
fig.savefig(fileName)

####
header("Finished")
