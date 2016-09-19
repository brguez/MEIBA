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
import matplotlib.patches as mpatches

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

## 2.1 Number of MEI per chromosome 
####################################
header("Number of MEI per chromosome")

# Iniciate lists
# 24 element lists. The 24 chromosomes (1, 2, 3, 4...Y)
L1List = [0] * 24
AluList = [0] * 24
SVAList = [0] * 24
ERVKList = [0] * 24

## Count number of MEI per chromosome
for MEIObj in VCFObj.lineList:
	supportingPElist = MEIObj.genotype.split(':')
	totalSupportingPE = int(supportingPElist[0]) + int(supportingPElist[1])

	if (MEIObj.chrom == "X"):
		index = 23 - 1

	elif (MEIObj.chrom == "Y"):
		index = 24 - 1

	else:
		index = int(MEIObj.chrom) - 1 

	
	if (MEIObj.infoDict['CLASS'] == 'L1'):		
		L1List[index] = L1List[index] + 1

	elif (MEIObj.infoDict['CLASS'] == 'Alu'):
		AluList[index] = AluList[index] + 1

	elif (MEIObj.infoDict['CLASS'] == 'SVA'):
		SVAList[index] = SVAList[index] + 1
	
	else:
		ERVKList[index] = ERVKList[index] + 1 

## Make figure
# Make bar plot 
xpos = np.arange(1, 25)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(8,6))
fig.suptitle('# MEI per chromosome')
plt.ylabel('# MEI', fontsize=12)

ax = fig.add_subplot(111)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

# Plot data
p1 = ax.bar(xpos, AluList, color='#A67D3D', alpha=0.75, edgecolor='#000000', width=width, align='center')
p2 = ax.bar(xpos, L1List, color='#008000', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=AluList)
p3 = ax.bar(xpos, SVAList, color='#87CEFA', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=[i+j for i,j in zip(AluList, L1List)])
p4 = ax.bar(xpos, ERVKList, color='#ff0000', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=[i+j+z for i,j,z in zip(AluList, L1List, SVAList)])

## Customize ticks
chrList = range(1, 23)
chrList = chrList + ['X', 'Y']

plt.yticks(np.arange(0, 2001, 100))
plt.xticks(xpos, chrList)

# Rotate X ticks:
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

## Make legend 
circle1 = mpatches.Circle((0, 0), 5, color='#A67D3D', alpha=0.75)
circle2 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.75)
circle3 = mpatches.Circle((0, 0), 5, color='#87CEFA', alpha=0.75)
circle4 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.75)
plt.figlegend((circle1, circle2, circle3, circle4), ('ALU', 'L1', 'SVA', 'ERVK'), loc = 'lower center', ncol=4, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/PCAWG_MEIperChr_barPlot.pdf"
plt.savefig(fileName)

## 2.2 MEI total number of supporting discordant paired-ends histogram 
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
plt.hist(L1PElist, bins=100, color='#008000', alpha=0.75)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 325)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax1.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 251, 10))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 326, 20))

### B) ALU
## Make plot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=16)
plt.hist(AluPElist, bins=125, color='#008000', alpha=0.75)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 2500)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax2.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 251, 10))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 2501, 200))

### C) SVA
## Make plot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=16)
plt.hist(SVAPElist, bins=100, color='#008000', alpha=0.75)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 75)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax3.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 251, 10))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 75, 5))

### D) ERVK
## Make plot
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("ERVK", fontsize=16)
plt.hist(ERVKPElist, bins=75, color='#008000', alpha=0.75)
plt.xlabel("# discordant paired-end", fontsize=14)
plt.ylabel("# MEI")
plt.xlim(0, 250)
plt.ylim(0, 6)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax4.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax4.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 251, 14))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)
plt.yticks(np.arange(0, 6, 1))

##### Save figure
fileName = outDir + "/PCAWG_totalSupportingPE_hist.pdf"
plt.savefig(fileName)

## 2.3 Number of + cluster vs - cluster discordant paired-ends scatterplot
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


## 2.4 MEI TSD length histogram 
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
plt.hist(L1TSDlist, bins=75, color='#008000', alpha=0.75)
plt.xlabel("TSD length", fontsize=14)
plt.ylabel("# MEI")

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax1.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 100, 5))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

### B) ALU
## Make plot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=16)
plt.hist(AluTSDlist, bins=75, color='#008000', alpha=0.75)
plt.xlabel("TSD length", fontsize=14)
plt.ylabel("# MEI")

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax2.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 100, 5))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

### C) SVA
## Make plot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=16)
plt.hist(SVATSDlist, bins=75, color='#008000', alpha=0.75)
plt.xlabel("TSD length", fontsize=14)
plt.ylabel("# MEI")

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax3.set_axisbelow(True)

## Customize ticks
plt.xticks(np.arange(0, 100, 5))
locs, labels = plt.xticks()
plt.setp(labels, rotation=30)

##### Save figure
fileName = outDir + "/PCAWG_TSDlen_hist.pdf"
plt.savefig(fileName)

## 2.5 MEI strand bar plot 
###########################
header("MEI strand bar plot")

## Gather data
L1StrandList = []
AluStrandList = []
SVAStrandList = []
ERVKStrandList = []

for MEIObj in VCFObj.lineList:
	
	if (MEIObj.infoDict['CLASS'] == 'L1'):		
		L1StrandList.append(MEIObj.infoDict['STRAND'])
	
	elif (MEIObj.infoDict['CLASS'] == 'Alu'):
		AluStrandList.append(MEIObj.infoDict['STRAND'])
	
	elif (MEIObj.infoDict['CLASS'] == 'SVA'):
		SVAStrandList.append(MEIObj.infoDict['STRAND'])

	else:
		ERVKStrandList.append(MEIObj.infoDict['STRAND'])

# Compute the number of times the MEI is inserted in + and - for each MEI type
nbPlusMinusL1 = [L1StrandList.count(i) for i in set(L1StrandList)]
nbPlusMinusAlu = [AluStrandList.count(i) for i in set(AluStrandList)]
nbPlusMinusSVA = [SVAStrandList.count(i) for i in set(SVAStrandList)]
nbPlusMinusERVK = [ERVKStrandList.count(i) for i in set(ERVKStrandList)]

# Convert into percentages. 
percPlusMinusL1 = [float(i)/sum(nbPlusMinusL1)*100 for i in nbPlusMinusL1]
percPlusMinusAlu = [float(i)/sum(nbPlusMinusAlu)*100 for i in nbPlusMinusAlu]
percPlusMinusSVA = [float(i)/sum(nbPlusMinusSVA)*100 for i in nbPlusMinusSVA]
percPlusMinusERVK = [float(i)/sum(nbPlusMinusERVK)*100 for i in nbPlusMinusERVK]

# Put percentages into two lists
tmpList = map(list, zip(percPlusMinusL1, percPlusMinusAlu, percPlusMinusSVA, percPlusMinusERVK)) 
plusList = tmpList[0]
minusList = tmpList[1]

## Make figure
# Make bar plot 
xpos = np.arange(4)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(5,6))
fig.suptitle('MEI DNA strand')
plt.ylabel('%', fontsize=12)
ax = fig.add_subplot(111)
p1 = ax.bar(xpos, plusList, color='#A67D3D', alpha=0.75, edgecolor='#000000', width=width, align='center')
p2 = ax.bar(xpos, minusList, color='#008000', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=plusList)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize ticks
plt.yticks(np.arange(0, 101, 10))
plt.xticks(xpos, ('L1', 'ALU', 'SVA', 'ERVK'))

## Make legend 
circle1 = mpatches.Circle((0, 0), 5, color='#A67D3D', alpha=0.75)
circle2 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.75)
plt.figlegend((circle1, circle2), ('+ strand', '- strand'), loc = 'lower center', ncol=2, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/PCAWG_strand_barPlot.pdf"
plt.savefig(fileName)

## 2.6 MEI structure bar plot 
##############################
header("MEI structure bar plot")

## Gather data
L1StructureList = []
AluStructureList = []
SVAStructureList = []

for MEIObj in VCFObj.lineList:
	
	if (MEIObj.infoDict['CLASS'] == 'L1'):		
		L1StructureList.append(MEIObj.infoDict['STRUCT'])
	
	elif (MEIObj.infoDict['CLASS'] == 'Alu'):
		AluStructureList.append(MEIObj.infoDict['STRUCT'])
	
	elif (MEIObj.infoDict['CLASS'] == 'SVA'):
		SVAStructureList.append(MEIObj.infoDict['STRUCT'])

# Compute the number of times the MEI are inserted in + and - for each MEI type
nbInvFullDelL1 = [L1StructureList.count(i) for i in set(L1StructureList)]
nbInvFullDelAlu = [AluStructureList.count(i) for i in set(AluStructureList)]
nbInvFullDelSVA = [SVAStructureList.count(i) for i in set(SVAStructureList)]

# Convert into percentages. 
percInvFullDelL1 = [float(i)/sum(nbInvFullDelL1)*100 for i in nbInvFullDelL1]
percInvFullDelAlu = [float(i)/sum(nbInvFullDelAlu)*100 for i in nbInvFullDelAlu]
percInvFullDelSVA = [float(i)/sum(nbInvFullDelSVA)*100 for i in nbInvFullDelSVA]

# Put percentages into two lists
tmpList = map(list, zip(percInvFullDelL1, percInvFullDelAlu, percInvFullDelSVA)) 

invList = tmpList[0]
delList = tmpList[2]
fullList = tmpList[1]

## Make figure
# Make bar plot 
xpos = np.arange(3)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(5,6))
fig.suptitle('MEI structure')
plt.ylabel('%', fontsize=12)
ax = fig.add_subplot(111)
p1 = ax.bar(xpos, invList, color='#A67D3D', alpha=0.75, edgecolor='#000000', width=width, align='center')
p2 = ax.bar(xpos, delList, color='#008000', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=invList)
p3 = ax.bar(xpos, fullList, color='#FFA500', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=[i+j for i,j in zip(invList, delList)])

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize ticks
plt.yticks(np.arange(0, 101, 10))
plt.xticks(xpos, ('L1', 'ALU', 'SVA'))

## Make legend 
circle1 = mpatches.Circle((0, 0), 5, color='#A67D3D', alpha=0.75)
circle2 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.75)
circle3 = mpatches.Circle((0, 0), 5, color='#FFA500', alpha=0.75)
plt.figlegend((circle1, circle2, circle3), ('5\' Inverted', '5\' Truncated', 'Full-length'), loc = 'lower center', ncol=3, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/PCAWG_structure_barPlot.pdf"
plt.savefig(fileName)

## 2.7 MEI in CPG 
###################
header("MEI CPG region bar plot")

## Gather data
TSGRegionList = []
oncogeneRegionList = []
bothRegionList = []

test = []
genes = []

for MEIObj in VCFObj.lineList:
# Nota: number of CPG genes inconsistent with input VCf. Investigate...	

	if 'CPG' in MEIObj.infoDict:

		test.append(MEIObj)

		# A) Tumor suppresor gene (TSG)		
		if (MEIObj.infoDict['ROLE'] == 'TSG'):
			
			TSGRegionList.append(MEIObj.infoDict['REGION'])

		# B) Oncogene
		elif (MEIObj.infoDict['ROLE'] == 'oncogene'):
			oncogeneRegionList.append(MEIObj.infoDict['REGION'])
		
		# C) Both
		else:
			bothRegionList.append(MEIObj.infoDict['REGION'])

# Compute the number of times the MEI are inserted in the different genic regions (Only introns for oncogene and both. For TSG two in UTR)
nbIntronUtrTSG = [TSGRegionList.count(i) for i in set(TSGRegionList)]
nbIntronUtrOncogene = [oncogeneRegionList.count(i) for i in set(oncogeneRegionList)]
nbIntronUtrBoth = [bothRegionList.count(i) for i in set(bothRegionList)]

## Add 0 UTR to the oncogene and both list
nbIntronUtrOncogene.append(0)
nbIntronUtrBoth.append(0)

# Put percentages into two lists
tmpList = map(list, zip(nbIntronUtrTSG, nbIntronUtrOncogene, nbIntronUtrBoth)) 

intronList = tmpList[0]
utrList = tmpList[1]

## Make figure
# Make bar plot 
xpos = np.arange(3)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(5,6))
fig.suptitle('MEI in CPG')
plt.ylabel('# MEI', fontsize=12)
ax = fig.add_subplot(111)
p1 = ax.bar(xpos, intronList, color='#A67D3D', alpha=0.75, edgecolor='#000000', width=width, align='center')
p2 = ax.bar(xpos, utrList, color='#008000', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=intronList)
plt.ylim(0, 50)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize ticks
plt.yticks(np.arange(0, 51, 5))
plt.xticks(xpos, ('TSG', 'ONCOGENE', 'BOTH'), fontsize=12)

## Make legend 
circle1 = mpatches.Circle((0, 0), 5, color='#A67D3D', alpha=0.75)
circle2 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.75)
plt.figlegend((circle1, circle2), ('INTRON', '3\' UTR'), loc = 'lower center', ncol=2, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/PCAWG_CPG_barPlot.pdf"
plt.savefig(fileName)

## 2.8 MEI insertion region pie chart
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


