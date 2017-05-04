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

    #print "** Insertion absent in reference genome"          

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

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        ancestry = donorIdAncestryDict[donorId]
    
        # print "TIO: ", donorId, ancestry, genotype
           
        # a) Heterozygous
        if (genotype == "0/1"):

            nbChromDict['PCAWG'] += 2
            nbChromDict[ancestry] += 2

            alleleCountDict['PCAWG'] += 1     
            alleleCountDict[ancestry] += 1     

            # print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]                
            # print "alleleCounts: ", alleleCountDict['PCAWG'], alleleCountDict[ancestry]

        # b) Homozygous alternative
        elif (genotype == "1/1"):
    
            nbChromDict['PCAWG'] += 2
            nbChromDict[ancestry] += 2

            alleleCountDict['PCAWG'] += 2
            alleleCountDict[ancestry] += 2

            # print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]               
            # print "alleleCounts: ", alleleCountDict['PCAWG'], alleleCountDict[ancestry]

        # c) Homozygous reference
        elif (genotype == "0/0"):
    
            nbChromDict['PCAWG'] += 2
            nbChromDict[ancestry] += 2
    
            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]    

        # d) Haploid carrier (males X and Y outside PAR region)
        elif (genotype == "1"):

            nbChromDict['PCAWG'] += 1
            nbChromDict[ancestry] += 1

            alleleCountDict['PCAWG'] += 1     
            alleleCountDict[ancestry] += 1     
                
            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]               
            #print "alleleCounts: ", alleleCountDict['PCAWG'], alleleCountDict[ancestry]
            
        # e) Haploid not carrier (males X and Y outside PAR region)
        elif (genotype == "0"):

            nbChromDict['PCAWG'] += 1
            nbChromDict[ancestry] += 1 
                
            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]               

    ## Compute overall and per echnicity variant allele frequencies
    alleleFreqDict = {}
     
    alleleFreqDict['PCAWG'] = float(alleleCountDict['PCAWG'])/float(nbChromDict['PCAWG'])
    
    for ancestry in ancestryCodesList:
        alleleFreqDict[ancestry] = float(alleleCountDict[ancestry])/float(nbChromDict[ancestry])

    return alleleCountDict, alleleFreqDict



def variantFrequencies_ref(MEIObj, donorIdAncestryDict):
    """
    For a variant in the reference genome compute its allele count and frequency for the complete PCAWG cohort and for each echnicity
    """

    #print "** Insertion in the reference genome"       

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

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        ancestry = donorIdAncestryDict[donorId]
    
        # print "TIO: ", donorId, ancestry, genotype

        # a) Heterozygous
        if (genotype == "0/1"):

            nbChromDict['PCAWG'] += 2
            nbChromDict[ancestry] += 2

            alleleCountDict['PCAWG'] += 1     
            alleleCountDict[ancestry] += 1     

            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]                
            #print "alleleCounts: ", alleleCountDict['PCAWG'], alleleCountDict[ancestry]

        # b) Homozygous alternative (MEI absent)
        elif (genotype == "1/1"):
    
            nbChromDict['PCAWG'] += 2
            nbChromDict[ancestry] += 2

            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]               

        # c) Homozygous reference (MEI present)
        elif (genotype == "0/0"):
    
            nbChromDict['PCAWG'] += 2
            nbChromDict[ancestry] += 2
    
            alleleCountDict['PCAWG'] += 2
            alleleCountDict[ancestry] += 2

            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]    
            #print "alleleCounts: ", alleleCountDict['PCAWG'], alleleCountDict[ancestry]

        # d) Haploid not carrier (males X and Y outside PAR region)
        elif (genotype == "1"):

            nbChromDict['PCAWG'] += 1
            nbChromDict[ancestry] += 1
 
            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]               
            
        # e) Haploid carrier (males X and Y outside PAR region)
        elif (genotype == "0"):

            nbChromDict['PCAWG'] += 1
            nbChromDict[ancestry] += 1 
                
            alleleCountDict['PCAWG'] += 1     
            alleleCountDict[ancestry] += 1     

            #print "chromCounts: ", nbChromDict['PCAWG'], nbChromDict[ancestry]               
            #print "alleleCounts: ", alleleCountDict['PCAWG'], alleleCountDict[ancestry]

    ## Compute overall and per echnicity variant allele frequencies
    alleleFreqDict = {}
     
    alleleFreqDict['PCAWG'] = float(alleleCountDict['PCAWG'])/float(nbChromDict['PCAWG'])
    
    for ancestry in ancestryCodesList:
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
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.vcf
metadata = args.metadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "metadata: ", metadata
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


#### 1. Read input multi-sample VCF and generate a VCF object
###############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Compute for each source element its allele count and frequency
########################################################################
# Allele count and frequency computed overall and across the different echnicities. 

header("2. Compute for each source element its allele count and frequency")

alleleCountDict = {}
alleleFreqDict = {}

## For each MEI
for MEIObj in VCFObj.lineList:

    ## MEI identifier
    ## A) MEI correspond to a germline source element -> use source element identifier
    if ('SRCID' in MEIObj.infoDict):

        MEIid = MEIObj.infoDict['SRCID']

    # B) MEI does not correspond a source element -> create coordinates based identifier
    else:

        MEIid = MEIObj.infoDict["CLASS"] + '_' + MEIObj.chrom + '_' + str(MEIObj.pos)
    
    # print "MEI-ID: ",  MEIid, MEIObj.ref, MEIObj.alt 

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
alleleCountDf = alleleCountDf[colOrder]

# Save output into tsv
outFilePath = outDir + '/germline_source_elements.alleleCount.tsv'
alleleCountDf.to_csv(outFilePath, sep='\t') 

### 3.1 MEI allele frequency 
# Create pandas dataframe from dictionary
alleleFreqDf = pd.DataFrame(alleleFreqDict) 

# transpose to have MEI as rows 
alleleFreqDf = alleleFreqDf.T 

# Reorder columns and remove UNK columns:
colOrder = ['PCAWG', 'EUR', 'ASN', 'AFR', 'SAN', 'AMR']
alleleFreqDf = alleleFreqDf[colOrder]

# Save output into tsv
outFilePath = outDir + '/germline_source_elements.alleleFreq.tsv'
alleleFreqDf.to_csv(outFilePath, sep='\t') 


#### End
header("FINISH!!")


