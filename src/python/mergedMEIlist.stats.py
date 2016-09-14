#!/usr/bin/env python
#coding: utf-8 

def header(string):
    """ 
        Display  header
    """ 
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, "****", string, "****"

def info(string):
    """ 
        Display basic information
    """ 
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string,

#### CLASSES ####


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
parser.add_argument('VCF', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 

###########
## Start ## 
###########

############################################
## 1. Create VCF object and read input VCF #
############################################
header("Read input VCF")
VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)

########################################
## 2. Plot information about MEI calls #
########################################

## 2.1 MEI total number of supporting discordant paired-ends histogram 
#######################################################################
header("Total number of discordant paired-ends histogram")

## Gather data
L1PElist = []
AluPElist = []
SVAPElist = []
ERVKPElist = []

for MEIObj in VCFObj.lineList:
	supportingPElist = MEIObj.genotype.split(':')
	totalSupportingPE = int(supportingPElist[0]) + int(supportingPElist[1])

	if (MEIObj.infoDict['CLASS'] == 'L1'):		
		L1PElist.append(totalSupportingPE)

	elif (MEIObj.infoDict['CLASS'] == 'Alu'):
		AluPElist.append(totalSupportingPE)
	
	elif (MEIObj.infoDict['CLASS'] == 'SVA'):
		SVAPElist.append(totalSupportingPE)
	else:
		ERVKPElist.append(totalSupportingPE)

#### Make plot
fig = plt.figure(figsize=(20,14))                
fig.suptitle('Total discordant paired-end support', fontsize=20)

### A) L1
## Make plot
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("LINE-1", fontsize=16)
plt.hist(L1PElist, bins=100, color='#008000', alpha=0.5)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 325)

## Customize ticks
plt.xticks(np.arange(0, 251, 10))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 326, 20))

### B) ALU
## Make plot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=16)
plt.hist(AluPElist, bins=125, color='#008000', alpha=0.5)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 2500)

## Customize ticks
plt.xticks(np.arange(0, 251, 10))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 2501, 200))

### C) SVA
## Make plot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=16)
plt.hist(SVAPElist, bins=100, color='#008000', alpha=0.5)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 75)

## Customize ticks
plt.xticks(np.arange(0, 251, 10))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 75, 5))

### D) ERVK
## Make plot
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("ERVK", fontsize=16)
plt.hist(ERVKPElist, bins=75, color='#008000', alpha=0.5)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 6)

## Customize ticks
plt.xticks(np.arange(0, 251, 14))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 6, 1))

##### Save figure
fileName = outDir + "/PCAWG_totalSupportingPE_hist.pdf"
plt.savefig(fileName)

## 2.2 Number of + cluster vs - cluster discordant paired-ends scatterplot (tomorrow)
###########################################################################
header("Number of + and - cluster supporting paired-ends scatterplots")


## Gather data
L1PEtuplelist = []
AluPEtuplelist = []
SVAPEtuplelist = []
ERVKPEtuplelist = []

for MEIObj in VCFObj.lineList:
	supportingPEtuple = tuple(map(int, MEIObj.genotype.split(':')))
	
	if (MEIObj.infoDict['CLASS'] == 'L1'):		
		L1PEtuplelist.append(supportingPEtuple)
	
	elif (MEIObj.infoDict['CLASS'] == 'Alu'):
		AluPEtuplelist.append(supportingPEtuple)
	
	elif (MEIObj.infoDict['CLASS'] == 'SVA'):
		SVAPEtuplelist.append(supportingPEtuple)
	
	else:
		ERVKPEtuplelist.append(supportingPEtuple)


### Make plot
fig = plt.figure(figsize=(10,14))                
fig.suptitle('Number of discordant paired-ends correlation', fontsize=18)

## L1
tmpList = map(list, zip(*L1PEtuplelist))
L1plusPE = tmpList[0]
L1minusPE = tmpList[1]

# Compute correlation
corr = stats.pearsonr(L1plusPE, L1minusPE)
coefficient = format(corr[0], '.3f') 
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

# Make scatterplot
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("LINE-1", fontsize=14)
plt.scatter(L1plusPE, L1minusPE, color='#008000', alpha=.4)
plt.xlim((0, (max(L1plusPE) + 10)))
plt.ylim((0, (max(L1minusPE) + 10)))
plt.xlabel('+ cluster', fontsize=12)
plt.ylabel('- cluster', fontsize=12)
ax1.text(0.5, 0.1, text, transform = ax1.transAxes)

## Alu
tmpList = map(list, zip(*AluPEtuplelist))
AluPlusPE = tmpList[0]
AluMinusPE = tmpList[1]

# Compute correlation
corr = stats.pearsonr(AluPlusPE, AluMinusPE)
coefficient = format(corr[0], '.3f') 
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

# Make scatterplot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=14)
plt.scatter(AluPlusPE, AluMinusPE, color='#008000', alpha=.4)
plt.xlim((0, (max(AluPlusPE) + 10)))
plt.ylim((0, (max(AluMinusPE) + 10)))
plt.xlabel('+ cluster', fontsize=12)
plt.ylabel('- cluster', fontsize=12)
ax2.text(0.5, 0.1, text, transform = ax2.transAxes)

## SVA
tmpList = map(list, zip(*SVAPEtuplelist))
SVAPlusPE = tmpList[0]
SVAMinusPE = tmpList[1]

# Compute correlation
corr = stats.pearsonr(SVAPlusPE, SVAMinusPE)
coefficient = format(corr[0], '.3f') 
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

# Make scatterplot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=14)
plt.scatter(SVAPlusPE, SVAMinusPE, color='#008000', alpha=.4)
plt.xlim((0, (max(SVAPlusPE) + 10)))
plt.ylim((0, (max(SVAMinusPE) + 10)))
plt.xlabel('+ cluster', fontsize=12)
plt.ylabel('- cluster', fontsize=12)
ax3.text(0.5, 0.1, text, transform = ax3.transAxes)

## ERVK
tmpList = map(list, zip(*ERVKPEtuplelist))
ERVKPlusPE = tmpList[0]
ERVKMinusPE = tmpList[1]

# Compute correlation
corr = stats.pearsonr(ERVKPlusPE, ERVKMinusPE)
coefficient = format(corr[0], '.3f') 
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

# Make scatterplot
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("ERVK", fontsize=14)
plt.scatter(ERVKPlusPE, ERVKMinusPE, color='#008000', alpha=.4)
plt.xlim((0, (max(ERVKPlusPE) + 10)))
plt.ylim((0, (max(ERVKMinusPE) + 10)))
plt.xlabel('+ cluster', fontsize=12)
plt.ylabel('- cluster', fontsize=12)
ax4.text(0.5, 0.1, text, transform = ax4.transAxes)

## Save figure
fileName = outDir + "/PCAWG_discordantPE_clusters_correlation.pdf"
plt.savefig(fileName)


## 2.3 MEI TSD length histogram 
################################
header("MEI TSD length histogram")

## Gather data
L1TSDlist = []
AluTSDlist = []
SVATSDlist = []

for MEIObj in VCFObj.lineList:
	
	if (MEIObj.infoDict['CLASS'] == 'L1'):		
		L1TSDlist.append(int(MEIObj.infoDict['TSLEN']))
	
	elif (MEIObj.infoDict['CLASS'] == 'Alu'):
		AluTSDlist.append(int(MEIObj.infoDict['TSLEN']))
	
	elif (MEIObj.infoDict['CLASS'] == 'SVA'):
		SVATSDlist.append(int(MEIObj.infoDict['TSLEN']))

#### Make plot
fig = plt.figure(figsize=(14,12))                
fig.suptitle('Target site duplication (TSD)', fontsize=20)

### A) L1
## Make plot
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("LINE-1", fontsize=16)
plt.hist(L1TSDlist, bins=75, color='#008000', alpha=0.5)
plt.xlabel("TSD length", fontsize=14)
plt.ylabel("# MEI")

## Customize ticks
plt.xticks(np.arange(0, 100, 5))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

### B) ALU
## Make plot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=16)
plt.hist(AluTSDlist, bins=75, color='#008000', alpha=0.5)
plt.xlabel("TSD length", fontsize=14)
plt.ylabel("# MEI")

## Customize ticks
plt.xticks(np.arange(0, 100, 5))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

### C) SVA
## Make plot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=16)
plt.hist(SVATSDlist, bins=75, color='#008000', alpha=0.5)
plt.xlabel("TSD length", fontsize=14)
plt.ylabel("# MEI")

## Customize ticks
plt.xticks(np.arange(0, 100, 5))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

##### Save figure
fileName = outDir + "/PCAWG_TSDlen_hist.pdf"
plt.savefig(fileName)

## 2.4 MEI strand bar plot (tomorrow)
###########################


## 2.5 MEI structure bar plot (tomorrow)
##############################


## 2.6 MEI insertion region pie chart
#######################################

header("MEI functional spectrum pie chart")

## Gather data

regionsList = []

for MEIObj in VCFObj.lineList:

	if (MEIObj.infoDict['REGION']=="splicing") or (MEIObj.infoDict['REGION']=="upstream,downstream") or (MEIObj.infoDict['REGION']=="upstream") or (MEIObj.infoDict['REGION']=="downstream"):
		region = "Other"
	
	elif (MEIObj.infoDict['REGION']=="UTR5") or (MEIObj.infoDict['REGION']=="UTR3"):
		region = "UTR"

	elif (MEIObj.infoDict['REGION']=="ncRNA_exonic") or (MEIObj.infoDict['REGION']=="ncRNA_intronic"):
		region = "ncRNA"

	else:
		region = MEIObj.infoDict['REGION']

	regionsList.append(region)

regionTuples = [(x, int(regionsList.count(x))) for x in set(regionsList)]
regionList =  [list(t) for t in zip(*regionTuples)]

labels = regionList[0]
sizes = regionList[1]

## Make plot
fig = plt.figure(figsize=(6,6))                
fig.suptitle('MEI functional spectrum', fontsize=16)
colors = ['#008000', '#A67D3D', '#87CEFA', '#ff0000', '#FFD700', '#FFA500']
patches, texts, perc = plt.pie(sizes, colors=colors, startangle=90, autopct='%1.1f%%', pctdistance=1.2, labeldistance=1)
plt.legend(patches, labels, loc="best", fontsize=11)

##### Save figure
fileName = outDir + "/PCAWG_funcSpectrum_piechart.pdf"
plt.savefig(fileName)


## End ##
print 
print "***** Finished! *****"
print 


