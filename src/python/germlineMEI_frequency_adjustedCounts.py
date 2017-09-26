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


def variantFrequencies(MEIObj, donorIdAncestryDict):
    """
    For a variant absent in the reference genome compute its allele count and frequency for the complete PCAWG cohort and for each echnicity
    """
    ancestryCodesList = set(donorIdAncestryDict.values())

    # Initialize allele count dictionary   
    alleleCountDict = {}
    alleleCountDict['PCAWG'] = 0

    for ancestry in ancestryCodesList:
        alleleCountDict[ancestry] = 0

    ## Total number of chromosome copies in the population
    # Will be used for computing the allele frequency
    # If no missing genotypes would be == Number of donors * 2 (diploid, two copies of a given chromosome)
    nbChromDict = {}
    nbChromDict['PCAWG'] = 0 

    for ancestry in ancestryCodesList:
        nbChromDict[ancestry] = 0

    # For each donor and genotype
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        # Select only whitelisted donors
        # if True:
        if (donorId not in blackDonorsList):

            genotypeFieldList = genotypeField.split(":")
            genotype = genotypeFieldList[0]

            ancestry = donorIdAncestryDict[donorId]
    
            # print "TIO: ", donorId, ancestry, genotype
           
            # a) Heterozygous
            if (genotype == "0/1"):

                #print donorId, ancestry, genotype, genotypeField
                nbChromDict['PCAWG'] += 2
                nbChromDict[ancestry] += 2

                alleleCountDict['PCAWG'] += 1     
                alleleCountDict[ancestry] += 1     

            # b) Homozygous alternative
            elif (genotype == "1/1"):
    
                #print donorId, ancestry, genotype, genotypeField
                nbChromDict['PCAWG'] += 2
                nbChromDict[ancestry] += 2

                alleleCountDict['PCAWG'] += 2
                alleleCountDict[ancestry] += 2

            # c) Homozygous reference
            elif (genotype == "0/0"):
    
                nbChromDict['PCAWG'] += 2
                nbChromDict[ancestry] += 2

            # d) Haploid carrier (males X and Y outside PAR region)
            elif (genotype == "1"):

                nbChromDict['PCAWG'] += 1
                nbChromDict[ancestry] += 1

                alleleCountDict['PCAWG'] += 1     
                alleleCountDict[ancestry] += 1     
         
            # e) Haploid not carrier (males X and Y outside PAR region)
            elif (genotype == "0"):

                nbChromDict['PCAWG'] += 1
                nbChromDict[ancestry] += 1 
                        

    ## Compute overall and per echnicity variant allele frequencies
    alleleFreqDict = {}
     
    # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
    if (nbChromDict['PCAWG'] == 0):
        alleleFreqDict['PCAWG'] = "UNK" 
   
    # b) Allele freq. estimation available 
    else:
        alleleFreqDict['PCAWG'] = float(alleleCountDict['PCAWG'])/float(nbChromDict['PCAWG'])
    
    for ancestry in ancestryCodesList:

        # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
        if (nbChromDict[ancestry] == 0):
            alleleFreqDict[ancestry] = "UNK" 
   
        # b) Allele freq. estimation available 
        else:    
            alleleFreqDict[ancestry] = float(alleleCountDict[ancestry])/float(nbChromDict[ancestry])

    ## A) Novel insertion
    if ('GERMDB' not in MEIObj.infoDict):

        alleleFreqDict["novelty"] = "novel"
  
    ## B) Not novel insertion
    else:

        alleleFreqDict["novelty"] = "known"


    return alleleCountDict, alleleFreqDict



def variantFrequencies_ref(MEIObj, donorIdAncestryDict):
    """
    For a variant in the reference genome compute its allele count and frequency for the complete PCAWG cohort and for each echnicity
    """

    ancestryCodesList = set(donorIdAncestryDict.values())

    # Initialize allele count dictionary   
    alleleCountDict = {}
    alleleCountDict['PCAWG'] = 0

    for ancestry in ancestryCodesList:
        alleleCountDict[ancestry] = 0

    ## Total number of chromosome copies in the population
    # Will be used for computing the allele frequency
    # If no missing genotypes would be == Number of donors * 2 (diploid, two copies of a given chromosome)
    nbChromDict = {}
    nbChromDict['PCAWG'] = 0 

    for ancestry in ancestryCodesList:
        nbChromDict[ancestry] = 0

    # For each donor and genotype
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():
    
        # Select only whitelisted donors
        # if True:
        if (donorId not in blackDonorsList):

            genotypeFieldList = genotypeField.split(":")
            genotype = genotypeFieldList[0]
            ancestry = donorIdAncestryDict[donorId]

            # a) Heterozygous
            if (genotype == "0/1"):
            
                #print donorId, ancestry, genotype, genotypeField

                nbChromDict['PCAWG'] += 2
                nbChromDict[ancestry] += 2

                alleleCountDict['PCAWG'] += 1     
                alleleCountDict[ancestry] += 1     

            # b) Homozygous alternative (MEI absent)
            elif (genotype == "1/1"):
    
                nbChromDict['PCAWG'] += 2
                nbChromDict[ancestry] += 2

            # c) Homozygous reference (MEI present)
            elif (genotype == "0/0"):

                #print donorId, ancestry, genotype, genotypeField
    
                nbChromDict['PCAWG'] += 2
                nbChromDict[ancestry] += 2
    
                alleleCountDict['PCAWG'] += 2
                alleleCountDict[ancestry] += 2

            # d) Haploid not carrier (males X and Y outside PAR region)
            elif (genotype == "1"):

                nbChromDict['PCAWG'] += 1
                nbChromDict[ancestry] += 1
            
            # e) Haploid carrier (males X and Y outside PAR region)
            elif (genotype == "0"):
    
                #print donorId, ancestry, genotype, genotypeField
    
                nbChromDict['PCAWG'] += 1
                nbChromDict[ancestry] += 1 
                
                alleleCountDict['PCAWG'] += 1     
                alleleCountDict[ancestry] += 1     

    ## Compute overall and per echnicity variant allele frequencies
    alleleFreqDict = {}
     
    # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
    if (nbChromDict['PCAWG'] == 0):
        alleleFreqDict['PCAWG'] = "UNK" 
   
    # b) Allele freq. estimation available 
    else:
        alleleFreqDict['PCAWG'] = float(alleleCountDict['PCAWG'])/float(nbChromDict['PCAWG'])
    
    for ancestry in ancestryCodesList:

        # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
        if (nbChromDict[ancestry] == 0):
            alleleFreqDict[ancestry] = "UNK" 
   
        # b) Allele freq. estimation available 
        else:    
            alleleFreqDict[ancestry] = float(alleleCountDict[ancestry])/float(nbChromDict[ancestry])


    return alleleCountDict, alleleFreqDict


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

#### 0. Create dictionary with donor id ancestry equivalences
###############################################################
header("0. Create dictionary with donor id ancestry equivalences")
metadataFile = open(metadata, 'r')

donorIdAncestryDict = {}

for line in metadataFile:

     # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        donorId = line[0]
        ancestry = line[4]

        donorIdAncestryDict[donorId] = ancestry

#print "donorIdAncestryDict: ", donorIdAncestryDict


### Generate list with PCAWG whitelisted donors
metadataFile = open(metadata, 'r')
blackDonorsList = []

for line in metadataFile:

     # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        donorId = line[0]
        status = line[1]

        if (status == "Excluded"):
            blackDonorsList.append(donorId)

#### 1. Read input multi-sample VCF and generate a VCF object
###############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Compute for each germline MEI that passes the filters its allele count and frequency
#############################################################################################
# Allele count and frequency computed overall and across the different echnicities. 

header("2. Compute for each germline MEI that passes the filters its allele count and frequency")

alleleCountDict = {}
alleleFreqDict = {}

## For each MEI
for MEIObj in VCFObj.lineList:

    ## Select only those MEI that passes all the filters
    if (MEIObj.filter == "PASS"):

        ## MEI identifier
        ## A) MEI correspond to a germline source element -> use source element identifier
        if ('SRCID' in MEIObj.infoDict):

            MEIid = MEIObj.infoDict['SRCID']
    
        # B) MEI does not correspond a source element -> create coordinates based identifier
        else:
    
            MEIid = MEIObj.infoDict["CLASS"] + '_' + MEIObj.chrom + '_' + str(MEIObj.pos)
            
        #print "MEIid: ", MEIid

        ## Compute MEI allele count and frequencies
        alleleCountDict[MEIid] = {}
        alleleFreqDict[MEIid] = {}
    
        # A) MEI absent in reference genome
        if (MEIObj.alt == "<MEI>"):
                     
            alleleCountDictTmp, alleleFreqDictTmp = variantFrequencies(MEIObj, donorIdAncestryDict)
    
            # Add MEI allele counts and frequencies to the dict
            alleleCountDict[MEIid] = alleleCountDictTmp
            alleleFreqDict[MEIid] = alleleFreqDictTmp
    
        # B) MEI in the reference genome 
        elif (MEIObj.ref == "<MEI>"):
            
            alleleCountDictTmp, alleleFreqDictTmp = variantFrequencies_ref(MEIObj, donorIdAncestryDict)
     
            # Add MEI allele counts and frequencies to the dict
            alleleCountDict[MEIid] = alleleCountDictTmp
            alleleFreqDict[MEIid] = alleleFreqDictTmp
    
        # C) Raise error...  
        else:
            msg="Incorrectly formated VCF line"
            info(msg)
     
    
#### 3. Convert dictionaries into dataframes and generate output table
########################################################################
# For allele count and frequency generate a table with the following format:
#                     PCAWG EUR ASN ...
# source_element1     X1    Y1  Z1
# source_element2     X2    Y2  Z2
# ...

header("3. Convert dictionaries into dataframes and generate output table")

### 3.1 MEI allele count 
# Create pandas dataframe from dictionary
alleleCountDf = pd.DataFrame(alleleCountDict) 

# transpose to have MEI as rows 
alleleCountDf = alleleCountDf.T 

# Reorder columns and remove UNK columns:
colOrder = ['PCAWG', 'EUR', 'ASN', 'AFR', 'SAN', 'AMR']
#colOrder = ['PCAWG', 'EUR', 'EAS', 'AFR', 'SAS', 'AMR'] ** used for 1KGP MEI genotypes
alleleCountDf = alleleCountDf[colOrder]

# Save output into tsv
outFilePath = outDir + '/' + fileName + '.alleleCount.tsv'
alleleCountDf.to_csv(outFilePath, sep='\t') 

### 3.1 MEI allele frequency 
# Create pandas dataframe from dictionary
alleleFreqDf = pd.DataFrame(alleleFreqDict) 

# transpose to have MEI as rows 
alleleFreqDf = alleleFreqDf.T 

# Reorder columns and remove UNK columns:
colOrder = ['PCAWG', 'EUR', 'ASN', 'AFR', 'SAN', 'AMR', 'novelty']
#colOrder = ['PCAWG', 'EUR', 'EAS', 'AFR', 'SAS', 'AMR' , 'novelty'] ** used for 1KGP MEI genotypes
alleleFreqDf = alleleFreqDf[colOrder]

# Save output into tsv
outFilePath = outDir + '/' + fileName + '.alleleFreq.tsv'
alleleFreqDf.to_csv(outFilePath, sep='\t') 


#### End
header("FINISH!!")


