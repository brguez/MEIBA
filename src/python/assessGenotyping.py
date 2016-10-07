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
    

def info(string):
    """ 
        Display basic information
    """ 
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


def genotypes2df(VCFObj):
	""" 
        
	""" 
	donorGtList = []

	## For each MEI in the VCF
	for MEIObj in VCFObj.lineList:
	
		# Create a series of genotype (donorId labeled) 
		MEIid = MEIObj.infoDict['CLASS'] + '_' + MEIObj.chrom + "_" + str(MEIObj.pos)
		donorGt =  pd.Series(MEIObj.genotypesDict, name=MEIid)

		# Add the series to the list of series
		donorGtList.append(donorGt)

	## Merge line series into dataframe (row <- donor_ids, columns <- MEI_ids):
	df1 = pd.concat(donorGtList, axis=1)

	## Transpose dataframe (row <- MEI_ids, columns <- donor_ids)
	df2 = df1.transpose()

	return df2

def gt2binary(gtString):
	"""

	"""
	genotype = gtString.split(':')[0]

	# A) Homozygous reference	
	if (genotype == '0') or (genotype == '0|0') or (genotype == '0/0'):
		boolean = 0
	
	# B) Heterozygous or homozygous MEI 
	else:
		boolean = 1

	return boolean


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

## Get user's input ## 
parser = argparse.ArgumentParser(description= "")
parser.add_argument('VCFPCAWG', help='Multi-sample VCF file containing genotyped MEI from PCAWG project')
parser.add_argument('VCF1KGP', help='Multi-sample VCF file containing genotyped MEI from 1000 genomes project')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCFPCAWG =  args.VCFPCAWG
VCF1KGP = args.VCF1KGP
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCFPCAWG: ", VCFPCAWG
print "VCF1KGP: ", VCF1KGP
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"


## Start ## 

#### 1. Read input VCF and generate VCF object
###############################################
header("1. Process input VCFs ")

## PCAWG VCF
VCFObjPCAWG = formats.VCF()
donorIdListPCAWG = VCFObjPCAWG.read_VCF_multiSample(VCFPCAWG)

## 1000 genomes multi-sample VCF
VCFObj1KGP = formats.VCF()
donorIdList1KGP = VCFObj1KGP.read_VCF_multiSample(VCF1KGP)

#### 2. Make genotyping binary dataframes
########################################
header("2. Make genotyping binary matrices")

#### PCAWG
gtDfPCAWG = genotypes2df(VCFObjPCAWG)
gtBinaryDfPCAWG = gtDfPCAWG.applymap(gt2binary)

#print "binary-PCAWG: ", gtBinaryDfPCAWG

#print '************************'

#### 1KGP
gtDf1KGP = genotypes2df(VCFObj1KGP)
gtBinaryDf1KGP = gtDf1KGP.applymap(gt2binary)

#print "binary-1KGP: ", gtBinaryDf1KGP

#### 3. Compare 1KGP and PCAWG genotyping binary matrices to
#############################################################
# a) compute the number of:
########################
# TP: true positives (PCAWG: 1, 1KGP: 1)
# TN: true negatives (PCAWG: 0, 1KGP: 0)
# FP: false positives (PCAWG: 1, 1KGP: 0)
# FN: false negatives (PCAWG: 0, 1KGP: 1)
# b) Make an output file per donor containing the false positive cases. One false positive case per row

header("3. Compare 1KGP and PCAWG genotyping binary matrices")

## Initialize counters
rowNamesList = list(gtBinaryDfPCAWG.index)
colNamesList = list(gtBinaryDfPCAWG.columns)

nbTPdict = {}
nbTNdict = {}
nbFPdict = {}
nbFNdict = {}

# For each ronorId 
for colName in colNamesList:
	nbTPdict[colName] = 0	
	nbTNdict[colName] = 0	 
	nbFPdict[colName] = 0	 
	nbFNdict[colName] = 0	 

print 'TPdict: ', nbTPdict 
print 'TNdict: ', nbTNdict 
print 'FPdict: ', nbFPdict 
print 'FNdict: ', nbFNdict 
print '************************'

## Create output files 
# Create dictionary with the following format:
# donorId_1 -> output filehandle
# ....
# donorId_n -> output filehandle

outFilesDict = {}
 
# For each donorId 
for colName in colNamesList:

	# Create output filehandler
	fileName = colName + '_FP.txt'
	outFilePath = outDir + '/' + fileName
	outFile = open(outFilePath, 'w')	
	
	# Write header
	row = 'chrom' + '\t' + 'pos' + '\t' + 'class' + '\n'
	outFile.write(row)	

	# Save filehandler into dictionary
	outFilesDict[colName] = outFile	

	
## Count number of TP, TN, FP, FN
# For each row name (rowName <- MEI_id)
for rowName in rowNamesList:

	#print "MEI: ", rowName

	# For each column name (columnName <- donor_id)
	for colName in colNamesList:

		#print "DONOR: ", colName		

		# a) True positive (TP)
		if (gtBinaryDfPCAWG.loc[rowName, colName] == 1) and (gtBinaryDf1KGP.loc[rowName, colName] == 1):
			# print 'TP', gtBinaryDfPCAWG.loc[rowName, colName], gtBinaryDf1KGP.loc[rowName, colName]			
			nbTPdict[colName] += 1
	
		# b) True negative (TN)
		elif (gtBinaryDfPCAWG.loc[rowName, colName] == 0) and (gtBinaryDf1KGP.loc[rowName, colName] == 0):
			#print 'TN', gtBinaryDfPCAWG.loc[rowName, colName], gtBinaryDf1KGP.loc[rowName, colName]			
			nbTNdict[colName] += 1

		# c) False positive (FP)
		elif (gtBinaryDfPCAWG.loc[rowName, colName] == 1) and (gtBinaryDf1KGP.loc[rowName, colName] == 0):
			#print 'FP', gtBinaryDfPCAWG.loc[rowName, colName], gtBinaryDf1KGP.loc[rowName, colName]			
			#print 'FP', rowName, colName 
			nbFPdict[colName] += 1
		
			# Save false positive into the output file:
			outFile = outFilesDict[colName] 
			rowNameList = rowName.split('_')
			row = rowNameList[1] + '\t' + rowNameList[2] + '\t' + rowNameList[0] + '\n'
			outFile.write(row)

		# d) False negative (FN)
		else:
			# print 'FN', gtBinaryDfPCAWG.loc[rowName, colName], gtBinaryDf1KGP.loc[rowName, colName]			
			nbFNdict[colName] += 1
			
print 'TPdict: ', nbTPdict 
print 'TNdict: ', nbTNdict 
print 'FPdict: ', nbFPdict 
print 'FNdict: ', nbFNdict 


#### 4. Compute recall and precision for TraFiC genotyping results 
##################################################################
# using 1000 genomes genotyping as reference 
############################################
### Recall (also known as sensitivity): is the fraction of relevant instances that are retrieved

# Recall = TP / (TP + FN)

### Precision: fraction of retrieved instances that are relevant

# Precision =  TP / (TP + FP)

# For each donor compute recall and precision
seriesList = []

for donorId in colNamesList:

	TP = nbTPdict[donorId]
	FN = nbFNdict[donorId]
	FP = nbFPdict[donorId]

	# Compute recall
	recall = float(TP) / (TP + FN)
	
	# Compute precision	
	precision = float(TP) / (TP + FP)

	print donorId, TP, FN, FP, recall, precision
	
	# Create series containing recall and precision for a given donor
	series =  pd.Series([recall, precision], index=['recall', 'precision'], name=donorId)

	# Add series to the list
	seriesList.append(series)


# print 'seriesList: ', seriesList

## Merge line series into dataframe (row <- recall and precision, columns <- donor_ids):
df1 = pd.concat(seriesList, axis=1)

## Transpose dataframe (row <- MEI_ids, columns <- donor_ids)
df2 = df1.transpose()

# print 'df1: ', df1
print 'df2: ', df2


#### 5. Make barplot with recall and precision
###############################################

fig = plt.figure()
ax = fig.add_subplot(111)

## Setting the positions and width for the bars
pos = np.array(range(len(df2['recall']))) 
width = 0.35

# Create a bar with recall data
recallPlt = ax.bar(pos,
        df2['recall'],
        # of width
        width,
        # with alpha 0.5
        alpha=0.75,
        # with color
        color='#008000',
)

# Create a bar with precision data
precisionPlt = ax.bar(pos+width,
        df2['precision'],
        # of width
        width,
        # with alpha 0.5
        alpha=0.75,
        # with color
        color='#A67D3D',
)

## axes and labels
ax.set_xlim(-width,len(pos)+width)
ax.set_ylim(0,1)
ax.set_ylabel('Measure', fontsize=12)
ax.set_title('Genotyping evaluation with 1KGP as reference', fontsize=14)
xTickMarks = list(df2.index)
ax.set_xticks(pos+width)
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, fontsize=12)

## add a legend
ax.legend( (recallPlt[0], precisionPlt[0]), ('recall', 'precision') )


## Save figure
fileName = outDir + "/PCAWGvs1KGP_genotyping_barplot.pdf"
plt.savefig(fileName)


#### END
header("Finished")

