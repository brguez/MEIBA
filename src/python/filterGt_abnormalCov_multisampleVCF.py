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
import pandas as pd

## Get user's input ##
parser = argparse.ArgumentParser(description= "Set as unknown genotype those genotypes with an abnormally high number of reads supporting the alternative allele. This are most likely false positives due to mapping artifacts")
parser.add_argument('inputVCF', help='multi-sample VCF file containing genotyped MEI')
parser.add_argument('outFileName', help='Output VCF name. In PCAWG we use the submitted_donor_id.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
outFileName = args.outFileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "outFileName: ", outFileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Read input multi-sample VCF and generate a VCF object
###############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)


#### 2. Organize data into dictionaries
#########################################
header("2. Organize data into dictionaries")

### I need to make two dictionaries with the following format:

## Dictionary 1
# Nested dict containing for each donor and MEI a list with MEI genotype and the number of reads supporting the alternative allele. 
# Key 1 (donorId1):
#              - key 2 (MEI1) -> list[genotype, NV]
genotypesDict = {}

## Dictionary 2:
# Dictionary containing for each MEI its index in the VCF object MEI list. 
# Key (MEIid): list_index
indexDict = {}
index = 0

## For each MEI
for MEIObj in VCFObj.lineList:

    nbGt = len(MEIObj.genotypesDict)

    MEIid = MEIObj.infoDict["CLASS"] + '_' + MEIObj.chrom + '_' + str(MEIObj.pos)
    
    # Add index to dictionary
    indexDict[MEIid] = index

    # For each donor and genotype
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]
        NV = genotypeFieldList[1]

        # Add donor and MEI with its genotype to dictionary
        if donorId not in genotypesDict:
            genotypesDict[donorId] = {}
        
        genotypesDict[donorId][MEIid] = {}

        genotypesDict[donorId][MEIid] = [genotype, NV]

    index += 1


#### 3. Create pandas dataframes from genotype dictionary 
############################################################
header("3. Create pandas dataframes from genotype dictionary")

## Create pandas dataframes from dictionaries
genotypesDf = pd.DataFrame(genotypesDict) 


#### 4. Compute the median number of reads supporting the alternative allele
#############################################################################
header("4. Compute the median number of reads supporting the alternative allele")

## Compute the median number of reads supporting the alternative allele for homozygous, heterozygous and haploid variants in each PCAWG donor
# Save the median values in dictionaries with the following format:
# Key (donorId): Number of reads supporting the alternative allele (NV)
 
medianNVHomDict = {}
medianNVHetDict = {}
medianNVHaplDict = {}

## For each donor 
for donorId, gtSerie in genotypesDf.iteritems():

    genotypesList = gtSerie.values.tolist()

    # Compute the median number of reads supporting the alternative allele for homozygous, heterozygous and 
    # haploialt variants
    medianNVHom = np.median([float(genotype[1]) for genotype in genotypesList if genotype[0] == '1/1'])
    medianNVHet = np.median([float(genotype[1]) for genotype in genotypesList if genotype[0] == '0/1'])
    medianNVHapl = np.median([float(genotype[1]) for genotype in genotypesList if genotype[0] == '1'])

    # Add the median number of the dictionary
    medianNVHomDict[donorId] = medianNVHom
    medianNVHetDict[donorId] = medianNVHet
    medianNVHaplDict[donorId] = medianNVHapl

#### 5. For each genotype assess if it has an abnormally 
#########################################################
# high number of reads supporting the alternative allele 
##########################################################
header("5. For each genotype assess if it has an abnormally high number of reads supporting the alternative allele ")

nbWrongGt = 0
wrongGtDict = {}

# For each donor
for donorId, gtSerie in genotypesDf.iteritems():
    
    # For each MEI in a given donor 
    for MEIid, genotypeInfo in gtSerie.iteritems():
        genotype = genotypeInfo[0]
        NV = float(genotypeInfo[1])
    
        ## Donor carring at least one copy of the alternative allele
        # If the donor does not carry the alternative makes no sense to apply the filter
        if (genotype == '1/1') or (genotype == '0/1') or (genotype == '1'):

            ## Compute the ratio, NV/NVmedian, where: 
            # NV:  Number of reads supporting the alternative allele for the current MEI
            # NVmedian: Median number of reads supporting the alternative allele for current MEI genotype state
            # a) Homozygous alternative        
            if (genotype == '1/1'):
                ratioNV = NV/medianNVHomDict[donorId]

            # b) Heterozygous
            elif (genotype == '0/1'):
                ratioNV = NV/medianNVHetDict[donorId]
    
            # c) Haploid variant
            else:
                ratioNV = NV/medianNVHaplDict[donorId]
 
            ## If 5 times more reads than the median set genotype as unknown             
            if (ratioNV >= 5):

                genotype = './.'

                # Count the total number of unknown genotypes and number of unknown genotypes per each variant. 
                nbWrongGt += 1

                if MEIid not in wrongGtDict:
                    wrongGtDict[MEIid] = 1
                else:
                    wrongGtDict[MEIid] += 1

                
                # Replace genotype in the VCF object
                MEIObj = VCFObj.lineList[indexDict[MEIid]]

                genotypeField = MEIObj.genotypesDict[donorId]
                genotypeFieldList = genotypeField.split(":")
                genotypeFieldList[0] = genotype
                MEIObj.genotypesDict[donorId] = ':'.join(genotypeFieldList) 


#### 6. Report the number of genotypes set as unkown 
header("6. Report the number of genotypes set as unkown ")

print "total-nb-unknown: ", nbWrongGt

print "nb-unknown-per-variant: "

for key, value in wrongGtDict.iteritems():
    print key, value


#### 7. Make output multisample VCF file
header("7. Make output multisample VCF file")

fileName = outFileName + '.vcf'
outFilePath = outDir + '/' + fileName

# 7.1 Write header
VCFObj.write_header(outFilePath)

# 7.2 Write variants
VCFObj.write_variants_multiSample(donorIdList, outFilePath)

####
header("Finished")
