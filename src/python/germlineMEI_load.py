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


def gt2binary(genotype):
    """
    Convert the genotype of a germline variant absent in the reference genome into binary: 1 (carrier) and 0 (not carrier). 
    """
   
    # A) Homozygous alternative, heterozygous or haploid carrier
    if (genotype == '1/1') or (genotype == '0/1') or (genotype == '1'):
        boolean = 1

    # B) Homozygous reference, haploid not carrier or unknown genotype
    else:
        boolean = 0

    return boolean

def gt2binary_ref(genotype):
    """
    Convert the genotype of a germline variant in the reference genome into binary: 1 (carrier) and 0 (not carrier). 
    """
   
    # A) Homozygous reference, heterozygous or haploid carrier
    if (genotype == '0/0') or (genotype == '0/1') or (genotype == '0'):
        boolean = 1

    # B) Homozygous alternative, haploid not carrier or unknown genotype
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
import formats

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('vcf', help='Multisample VCF containing genotyped MEI')
parser.add_argument('metadata', help='PCAWG donor metadata')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.vcf
metadata = args.metadata
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "metadata: ", metadata
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 0. Make metadata dataframe
###############################
header("0. Make metadata dataframe")

metadataFile = open(metadata, 'r')

metadataDict = {}

for line in metadataFile:

     # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        donorId = line[1]
        ancestry = line[4]
        projectCode = line[9]

        histologyCount	= line[10]
        histologyExclusion = line[11]
        tumorHistology = line[12]

        ## Exclude donors for tumor type analysis (total: 30 + 10 + 32 = 72). 3 possible exclusion reasons:
        # - Unknown histology (30)
        # - More than 1 possible histology (10)
        # - Excluded histology cohort (32)
        if (tumorHistology != "UNK") and (histologyCount == "1") and (histologyExclusion == "included") :
         
            metadataDict[donorId] = {}
            metadataDict[donorId]['ancestry'] = ancestry
            metadataDict[donorId]['projectCode'] = projectCode
            metadataDict[donorId]['tumorHistology'] = tumorHistology       

metadataDf = pd.DataFrame(metadataDict) 
metadataDf = metadataDf.T


#### 1. Read input multi-sample VCF and generate a VCF object
###############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Build dictionaries with donor genotypes
################################################
# Split filtered variants into two different dictionaries depending on if they are absent or not in the reference genome. 
# Nested dictionary format:
# key1 (MEIid) -> value (dict2)
#                   key2 (donorId)     ->   value (genotype)      

header("2. Build dictionaries with donor genotypes")

genotypesDict = {} # genotypes for variants absent the reference genome
genotypesRefDict = {} # genotypes for variants in the reference genome

## For each MEI
for MEIObj in VCFObj.lineList:

    ## Select only those MEI that passes all the filters
    if (MEIObj.filter == "PASS"):

        ## MEI identifier
        # A) MEI corresponds to a germline source element -> use source element identifier
        if ('SRCID' in MEIObj.infoDict):

            MEIid = MEIObj.infoDict['SRCID']

        # B) MEI does not correspond a source element -> create coordinates based identifier
        else:

            MEIid = MEIObj.infoDict["CLASS"] + '_' + MEIObj.chrom + '_' + str(MEIObj.pos)

        ## Split variants in two different dictionaries:
        # A) MEI absent in reference genome
        if (MEIObj.alt == "<MEI>"):

            genotypesDict[MEIid] = {}
     
            # For each donor and genotype
            for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

                # Discard excluded donor donors:
                if (donorId in metadataDict):

                    genotypeFieldList = genotypeField.split(":")
                    genotype = genotypeFieldList[0]    
                    genotypesDict[MEIid][donorId] = genotype

        ## B) MEI in the reference genome 
        elif (MEIObj.ref == "<MEI>"):

            genotypesRefDict[MEIid] = {}
     
            # For each donor and genotype
            for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

                # Discard excluded donor donors:
                if (donorId in metadataDict):

                    genotypeFieldList = genotypeField.split(":")
                    genotype = genotypeFieldList[0]
                    genotypesRefDict[MEIid][donorId] = genotype

        ## C) Raise error...  
        else:
            msg="Incorrectly formated VCF line"
            info(msg)
 

#### 3. Convert dictionaries into dataframes specifying 
########################################################
# donor status (1:carrier, 0:not_carrier)
###########################################

header("3. Convert dictionaries into dataframes specifying donor status")

### 3.1 Variant absent in reference genome
# a) No absent variants 
if not genotypesDict:
    boolAbsent = False

# b) There are absent variants
else:
    boolAbsent = True
    genotypesDf = pd.DataFrame(genotypesDict) 
    genotypesDf = genotypesDf.T
    genotypesAbsBinaryDf = genotypesDf.applymap(gt2binary)

### 3.2 Variant in the reference genome 
# a) No variants in the reference genome
if not genotypesRefDict:
    boolRef = False

# b) There are variants in the reference genome
else:
    boolRef = True
    genotypesRefDf = pd.DataFrame(genotypesRefDict) 
    genotypesRefDf = genotypesRefDf.T
    genotypesRefBinaryDf = genotypesRefDf.applymap(gt2binary_ref)
    

#### 4. Compute the number of different insertions passing the filters each donor carries
##########################################################################################

header("4. Compute the number of different insertions passing the filters each donor carries")

# a) There are insetions both absent and in the reference genome
if (boolAbsent) and (boolRef):

    # Concatenate dataframes
    dataframeList = [genotypesAbsBinaryDf, genotypesRefBinaryDf]
    genotypesAllBinaryDf = pd.concat(dataframeList)

# b) There are only insertions absent in the reference genome
elif (boolAbsent):
    genotypesAllBinaryDf = genotypesAbsBinaryDf

# c) There are only insertions in the reference genome
elif (boolRef):
    genotypesAllBinaryDf = genotypesRefBinaryDf

# d) There are not insertions passing the filters in the VCF
else:
    info("There are not insertions passing the filters in the VCF. Exit")
    exit(0)

## Compute the total number of MEI per donor 
nbMEIperDonorSeries = genotypesAllBinaryDf.sum(axis=0)
outputDf = metadataDf.assign(nbMEI=nbMEIperDonorSeries.values)


#### 5. Save into output file
##############################

header("5. Save into output file")
outFilePath = outDir + '/' + fileName + '.tsv'
outputDf.to_csv(outFilePath, sep='\t') 

#### End
header("FINISH!!")


