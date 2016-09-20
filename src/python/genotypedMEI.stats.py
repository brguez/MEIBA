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
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Read metadata file
# Make three dictionaries with the following structure:
# - dict1a: key(projectId) -> dict2a: key(donorId) -> nbHeterozygousMEI
# - dict1b: key(projectId) -> dict2b: key(donorId) -> nbHomozygousMEI
# - dict1c: key(donorId) -> projectCode
# Note: nbHeterozygousMEI and nbHomozygousMEI are the number of heterozygous and homozygous MEI for a given donor
# Note: nbHeterozygousMEI and nbHomozygousMEI are iniciated with a value of 0
header("2. Read metadata file")
metadataFile = open(metadata, 'r')

nbHeterozMEIDict = {} 
nbHomozMEIDict = {}
donorIdProjectCodes = {}

for line in metadataFile:
	    
	line = line.rstrip('\n')
	line = line.split('\t')

	projectCode = line[0]
	donorId = line[1]

	## Dict1a and Dict1b
	if projectCode not in nbHeterozMEIDict:

		# Create dictionary
		nbHeterozMEIDict[projectCode] = {}
		nbHomozMEIDict[projectCode] = {}

	# Initialize to 0 values:
	nbHeterozMEIDict[projectCode][donorId] = 0
	nbHomozMEIDict[projectCode][donorId] = 0

	## Dict1c
	donorIdProjectCodes[donorId] = projectCode

#for projectCode in nbHeterozMEIDict:
#	print len(nbHeterozMEIDict[projectCode]), projectCode	



#### 3. Compute parameters:
# - Variant allele frequency per MEI
# - Variant allele count per MEI
# - Number of heterozygous MEI per donor (Save into dictionary generated in step 2)
# - Number of homozygous MEI per donor (Save into dictionary generated in step 2) 

L1AlleleCount = []
L1AlleleFreq = []
AluAlleleCount = []
AluAlleleFreq = []
SvaAlleleCount = []
SvaAlleleFreq = []
ErvkAlleleCount = []
ErvkAlleleFreq = []

header("3. Compute parameters")

for MEIObj in VCFObj.lineList:
	
	## Total number of chromosome copies in the population
	# Number of donors * 2 (diploid, two copies of a given chromosome)
	totalNbChrom = len(MEIObj.genotypesDict) * 2

	## Compute MEI allele count and update counters to heterozygous and homozygous MEI per donor
	# MEI allele count: number of chromosomes in the population harvouring the MEI
	alleleCount = 0

 	for donorId, genotypeField in MEIObj.genotypesDict.iteritems():		

		projectCode = donorIdProjectCodes[donorId]
		genotype = genotypeField.split(":")[0]

		#print "donor_genotype: ", projectCode, donorId, genotype

		# A) Heterozygous
		if (genotype == "0/1"):
			alleleCount +=  1
			nbHeterozMEIDict[projectCode][donorId] += 1 		
	
		# B) Homozygous		
		elif (genotype == "1/1"):
			alleleCount += 2
			nbHomozMEIDict[projectCode][donorId] += 1 

	## Compute MEI allele frequency:
	alleleFrequency = float(alleleCount) / totalNbChrom
	
	## Save into list. One per MEI type	
	if (MEIObj.infoDict['CLASS'] == 'L1'):		
		L1AlleleCount.append(alleleCount)
		L1AlleleFreq.append(alleleFrequency)

	elif (MEIObj.infoDict['CLASS'] == 'Alu'):
		AluAlleleCount.append(alleleCount)
		AluAlleleFreq.append(alleleFrequency)
	
	elif (MEIObj.infoDict['CLASS'] == 'SVA'):
		SvaAlleleCount.append(alleleCount)
		SvaAlleleFreq.append(alleleFrequency)
	
	else:
		ErvkAlleleCount.append(alleleCount)
		ErvkAlleleFreq.append(alleleFrequency)


#### 4. Make plots:
header("4. Make plots")

# - Variant allele frequencies histogram across PCAWG donors (done)
# - Variant allele counts line graph (to do)
# - Number of MEI per donor and tumor type boxplot (to do)
# - Variant allele frequencies histogram for homozygous and heterozygous variants  based on read support during variant calling (to do) 

### 4.1 Variant allele frequencies histogram across PCAWG donors 
header("4.1 Make variant allele frequencies plot")

fig = plt.figure(figsize=(16,13))                
fig.suptitle('Variant allele frequencies (VAF) across PCAWG donors', fontsize=20)

### A) L1
## Make plot
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("LINE-1", fontsize=16)
plt.hist(L1AlleleFreq, bins=40, color='#008000', alpha=0.75)
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

### B) ALU
## Make plot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=16)
plt.hist(AluAlleleFreq, bins=40, color='#008000', alpha=0.75)
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

### C) SVA
## Make plot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=16)
plt.hist(SvaAlleleFreq, bins=40, color='#008000', alpha=0.75)
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

### D) ERVK (no ERVK for now.. we need to improve script for identifying ERVK TSD)
## Make plot
#ax4 = fig.add_subplot(2, 2, 4)
#ax4.set_title("ERVK", fontsize=16)
#plt.hist(ErvkAlleleFreq, bins=40, color='#008000', alpha=0.75)
#plt.xlabel("VAF", fontsize=14)
#plt.ylabel("# MEI")
#plt.xlim(0, 1)
#plt.ylim(0, max(ErvkAlleleFreq) + 10)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
#ax4.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)
#ax4.set_axisbelow(True)

## Customize ticks
#plt.xticks(np.arange(0, 1.01, 0.1))
#locs, labels = plt.xticks()
#plt.setp(labels, rotation=30)
#plt.yticks(np.arange(0, max(ErvkAlleleFreq)))

##### Save figure
fileName = outDir + "/PCAWG_cohortVAFs_hist.pdf"
plt.savefig(fileName)

### 4.2 Make Number of heterozygous and homozigous MEI per donor and tumor type boxplots
header("4.2 Make boxplots")

#### Organize the data for plotting
tupleListHeterozMEI = []
tupleListHomozMEI = []

for projectCode in sorted(nbHeterozMEIDict):
	 
	# Make tuple (projectCode, number of MEI per donor) for heterozygous and homozygous variants	
	projectCodeNbHeterozMEITuple = (projectCode, nbHeterozMEIDict[projectCode].values())
	projectCodeNbHomozMEITuple = (projectCode, nbHomozMEIDict[projectCode].values())

	# Add tuple to the list
	tupleListHeterozMEI.append(projectCodeNbHeterozMEITuple)
	tupleListHomozMEI.append(projectCodeNbHomozMEITuple)

## Make nested list with the following format:
# [donor1_nbMEI, donor2_nbMEI, ..., donorN_nbMEI], [donor1_nbMEI, donor2_nbMEI, ..., donorN_nbMEI] , ... [donor1_nbMEI, donor2_nbMEI, ..., donorN_nbMEI] 
# 		   project1_list				  project2_list						projectN_list

# Heterozygous MEI
tmpList = map(list, zip(*tupleListHeterozMEI))
projectCodesListHet = tmpList[0]
nbMEIPerDonorHet = tmpList[1]

# Homozygous MEI
tmpList = map(list, zip(*tupleListHomozMEI))
projectCodesListHom = tmpList[0]
nbMEIPerDonorHom = tmpList[1]

#### Make boxplot

fig = plt.figure(figsize=(18,12))           
fig.suptitle('# MEI per donor', fontsize=24)

### a) Heterozygous MEI
ax1 = fig.add_subplot(2, 1, 1)

# Create the boxplot
bp = ax1.boxplot(nbMEIPerDonorHet)
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
ax1.set_xticklabels(projectCodesListHet)
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

### b) Homozygous MEI
ax2 = fig.add_subplot(2, 1, 2)

# Create the boxplot
bp = ax2.boxplot(nbMEIPerDonorHom)
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
ax2.set_xticklabels(projectCodesListHom)
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

#### Save figure

fileName = outDir + "/PCAWG_nbMEIperDonor_boxplot.pdf"
fig.savefig(fileName)

header("Finished")
