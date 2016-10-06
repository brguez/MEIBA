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

#### 2. Make genotyping binary matrices
########################################
header("2. Make genotyping binary matrices")

#### PCAWG
gtDfPCAWG = genotypes2df(VCFObjPCAWG)
gtBinaryDfPCAWG = gtDfPCAWG.applymap(gt2binary)

print "binary-PCAWG: ", gtBinaryDfPCAWG

print '************************'

#### 1KGP
gtDf1KGP = genotypes2df(VCFObj1KGP)
gtBinaryDf1KGP = gtDf1KGP.applymap(gt2binary)

print "binary-1KGP: ", gtBinaryDf1KGP

#### 3. Compare 1KGP and PCAWG genotyping binary matrices
###########################################################
# Compute the number of:
# TP: true positives (PCAWG: 1, 1KGP: 1)
# TN: true negatives (PCAWG: 0, 1KGP: 0)
# FP: false positives (PCAWG: 1, 1KGP: 0)
# FN: false negatives (PCAWG: 0, 1KGP: 1)

header("3. Compare 1KGP and PCAWG genotyping binary matrices")

rowNamesList = list(gtBinaryDfPCAWG.index)
colNamesList = list(gtBinaryDfPCAWG.columns)


## Initialize counters
nbTPdict = {}
nbTNdict = {}
nbFPdict = {}
nbFNdict = {}

for colName in colNamesList:
	nbTPdict[colName] = 0	
	nbTNdict[colName] = 0	 
	nbFPdict[colName] = 0	 
	nbFNdict[colName] = 0	 

print 'TPdict: ', nbTPdict 
print 'TNdict: ', nbTNdict 
print 'FPdict: ', nbFPdict 
print 'FNdict: ', nbFNdict 

## Count number of TP, TN, FP, FN
# For each row name 
for rowName in rowNamesList:

	# For each column name
	for colName in colNamesList:

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
			# print 'FP', gtBinaryDfPCAWG.loc[rowName, colName], gtBinaryDf1KGP.loc[rowName, colName]			
			# print 'FP', rowName, colName 
			nbFPdict[colName] += 1

		# d) False negative (FN)
		else:
			# print 'FN', gtBinaryDfPCAWG.loc[rowName, colName], gtBinaryDf1KGP.loc[rowName, colName]			
			nbFNdict[colName] += 1
			

print 'TPdict: ', nbTPdict 
print 'TNdict: ', nbTNdict 
print 'FPdict: ', nbFPdict 
print 'FNdict: ', nbFNdict 

#### END
header("Finished")
