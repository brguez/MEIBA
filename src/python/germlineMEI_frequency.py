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
    For a variant absent in the reference genome compute its allele count and frequency for the complete cohort and for each echnicity
    """

    ancestryCodesList = set(donorIdAncestryDict.values())

    ## Initialize alternative allele count dictionary   
    altAlleleCountDict = {}
    altAlleleCountDict['COHORT'] = 0

    for ancestry in ancestryCodesList:
        altAlleleCountDict[ancestry] = 0

    ## Initialize total allele count dictionary
    # If no missing genotypes, for variants in automosomal chromosomes would be == Number of donors * 2 (diploid, two alleles in a given chromosome)
    totalAlleleCountDict = {}
    totalAlleleCountDict['COHORT'] = 0 

    for ancestry in ancestryCodesList:
        totalAlleleCountDict[ancestry] = 0

    ## For each donor and genotype
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        ancestry = donorIdAncestryDict[donorId]
           

        # a) Homozygous alternative
        if (genotype == "1/1"):
    
            totalAlleleCountDict['COHORT'] += 2
            totalAlleleCountDict[ancestry] += 2

            altAlleleCountDict['COHORT'] += 2
            altAlleleCountDict[ancestry] += 2

        # b) Heterozygous
        elif (genotype == "0/1"):

            totalAlleleCountDict['COHORT'] += 2
            totalAlleleCountDict[ancestry] += 2

            altAlleleCountDict['COHORT'] += 1     
            altAlleleCountDict[ancestry] += 1     

       
        # c) Haploid carrier (males X and Y outside PAR region)
        elif (genotype == "1"):

            totalAlleleCountDict['COHORT'] += 1
            totalAlleleCountDict[ancestry] += 1

            altAlleleCountDict['COHORT'] += 1     
            altAlleleCountDict[ancestry] += 1     
         
        # d) Homozygous reference
        elif (genotype == "0/0"):
    
            totalAlleleCountDict['COHORT'] += 2
            totalAlleleCountDict[ancestry] += 2

        # e) Haploid not carrier (males X and Y outside PAR region)
        elif (genotype == "0"):

            totalAlleleCountDict['COHORT'] += 1
            totalAlleleCountDict[ancestry] += 1 
                

    ## Compute overall and per echnicity variant allele frequencies
    alleleFreqDict = {}
     
    # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
    if (totalAlleleCountDict['COHORT'] == 0):
        alleleFreqDict['COHORT'] = "UNK" 
   
    # b) Allele freq. estimation available 
    else:
        alleleFreqDict['COHORT'] = float(altAlleleCountDict['COHORT'])/float(totalAlleleCountDict['COHORT'])
    
    for ancestry in ancestryCodesList:

        # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
        if (totalAlleleCountDict[ancestry] == 0):
            alleleFreqDict[ancestry] = "UNK" 
   
        # b) Allele freq. estimation available 
        else:    
            alleleFreqDict[ancestry] = float(altAlleleCountDict[ancestry])/float(totalAlleleCountDict[ancestry])

    return altAlleleCountDict, totalAlleleCountDict, alleleFreqDict


def variantFrequencies_ref(MEIObj, donorIdAncestryDict):
    """
    For a variant in the reference genome compute its allele count and frequency for the complete cohort and for each echnicity
    """

    ancestryCodesList = set(donorIdAncestryDict.values())

    ## Initialize alternative allele count dictionary   
    altAlleleCountDict = {}
    altAlleleCountDict['COHORT'] = 0

    for ancestry in ancestryCodesList:
        altAlleleCountDict[ancestry] = 0

    ## Initialize total allele count dictionary
    # If no missing genotypes, for variants in automosomal chromosomes would be == Number of donors * 2 (diploid, two alleles in a given chromosome)
    totalAlleleCountDict = {}
    totalAlleleCountDict['COHORT'] = 0 

    for ancestry in ancestryCodesList:
        totalAlleleCountDict[ancestry] = 0

    # For each donor and genotype
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        ancestry = donorIdAncestryDict[donorId]

        # a) Homozygous reference (MEI present)
        if (genotype == "0/0"):
    
            totalAlleleCountDict['COHORT'] += 2
            totalAlleleCountDict[ancestry] += 2
    
            altAlleleCountDict['COHORT'] += 2
            altAlleleCountDict[ancestry] += 2

        # b) Heterozygous
        elif (genotype == "0/1"):

            totalAlleleCountDict['COHORT'] += 2
            totalAlleleCountDict[ancestry] += 2

            altAlleleCountDict['COHORT'] += 1     
            altAlleleCountDict[ancestry] += 1     

        # c) Haploid carrier (males X and Y outside PAR region)
        elif (genotype == "0"):

            totalAlleleCountDict['COHORT'] += 1
            totalAlleleCountDict[ancestry] += 1 
                
            altAlleleCountDict['COHORT'] += 1     
            altAlleleCountDict[ancestry] += 1  

        # d) Homozygous alternative (MEI absent)
        elif (genotype == "1/1"):
    
            totalAlleleCountDict['COHORT'] += 2
            totalAlleleCountDict[ancestry] += 2

        # e) Haploid not carrier (males X and Y outside PAR region)
        elif (genotype == "1"):

            totalAlleleCountDict['COHORT'] += 1
            totalAlleleCountDict[ancestry] += 1
            

    ## Compute overall and per echnicity variant allele frequencies
    alleleFreqDict = {}
     
    # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
    if (totalAlleleCountDict['COHORT'] == 0):
        alleleFreqDict['COHORT'] = "UNK" 
   
    # b) Allele freq. estimation available 
    else:
        alleleFreqDict['COHORT'] = float(altAlleleCountDict['COHORT'])/float(totalAlleleCountDict['COHORT'])
    
    for ancestry in ancestryCodesList:

        # a) Allele freq. estimation not available for those insertions with unknown genotype in all the donors
        if (totalAlleleCountDict[ancestry] == 0):
            alleleFreqDict[ancestry] = "UNK" 
   
        # b) Allele freq. estimation available 
        else:    
            alleleFreqDict[ancestry] = float(altAlleleCountDict[ancestry])/float(totalAlleleCountDict[ancestry])


    return altAlleleCountDict, totalAlleleCountDict, alleleFreqDict


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
parser.add_argument('metadata', help='COHORT donor metadata')
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

        donorId = line[1]
        ancestry = line[4]

        donorIdAncestryDict[donorId] = ancestry


#### 1. Read input multi-sample VCF and generate a VCF object
###############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)


#### 2. Compute for each germline MEI that passes the filters its allele count and frequency
#############################################################################################
# Allele count and frequency computed overall and for the different echnicities. 

header("2. Compute for each germline MEI that passes the filters its allele count and frequency")

altAlleleCountDict = {}
alleleFreqDict = {}
totalAlleleCountDict = {}

## For each MEI
for MEIObj in VCFObj.lineList:

    ## Select only those MEI that passes all the filters or without filtering info
    if (MEIObj.filter == "PASS") or (MEIObj.filter == "."):

        ## MEI identifier
        MEIid = ''

        ## A) MEI correspond to a germline source element -> use source element identifier
        if ('SRCID' in MEIObj.infoDict):

            MEIid = MEIObj.infoDict['SRCID']
    
        # B) MEI does not correspond a source element -> create coordinates based identifier
        else:
    
            MEIid = MEIObj.infoDict["CLASS"] + '_' + MEIObj.chrom + '_' + str(MEIObj.pos)

        ## Compute MEI allele counts and frequencies
        altAlleleCountDict[MEIid] = {}
        totalAlleleCountDict[MEIid] = {}
        alleleFreqDict[MEIid] = {}

        # A) MEI absent in reference genome
        if (MEIObj.alt == "<MEI>"):
                     
            altAlleleCountDict[MEIid], totalAlleleCountDict[MEIid], alleleFreqDict[MEIid] = variantFrequencies(MEIObj, donorIdAncestryDict)
    
        # B) MEI in the reference genome 
        elif (MEIObj.ref == "<MEI>"):
            
            altAlleleCountDict[MEIid], totalAlleleCountDict[MEIid], alleleFreqDict[MEIid] = variantFrequencies_ref(MEIObj, donorIdAncestryDict)
     
        # C) Raise error 
        else:
            msg="Incorrectly formated VCF line"
            info(msg)
            continue

        ### Add allele counts and frequencies to the info dictionary in the VCF entry
        MEIObj.infoDict["AC"] = altAlleleCountDict[MEIid]['COHORT']
        MEIObj.infoDict["AN"] = totalAlleleCountDict[MEIid]['COHORT']   
        MEIObj.infoDict["AF"] = alleleFreqDict[MEIid]['COHORT']                
        MEIObj.infoDict["AFR_AF"] = alleleFreqDict[MEIid]['AFR']
        MEIObj.infoDict["AMR_AF"] = alleleFreqDict[MEIid]['AMR']
        MEIObj.infoDict["EAS_AF"] = alleleFreqDict[MEIid]['ASN']
        MEIObj.infoDict["EUR_AF"] = alleleFreqDict[MEIid]['EUR']
        MEIObj.infoDict["SAS_AF"] = alleleFreqDict[MEIid]['SAN']
        
        ## Update info string
        MEIObj.info = MEIObj.make_info()


#### 3. Write multisample VCF with new fields into a file
############################################################

header("3. Write multisample VCF with new fields into a file")

outFilePath = outDir + '/' + fileName + '.freq.vcf'

## Write header
VCFObj.write_header(outFilePath)

## Write variants
VCFObj.write_variants_multiSample(donorIdList, outFilePath)
    

#### 4. Write summary tables with allele counts and frequencies
#################################################################
# For allele count and frequency generate a table with the following format:
#                     COHORT EUR ASN ...
# source_element1     X1    Y1  Z1
# source_element2     X2    Y2  Z2
# ...

header("4. Write summary tables with allele counts and frequencies")

### 4.1 MEI allele count 
# Create pandas dataframe from dictionary
alleleCountDf = pd.DataFrame(altAlleleCountDict) 

# transpose to have MEI as rows 
alleleCountDf = alleleCountDf.T 

print 'alleleCountDf: ', alleleCountDf

# Reorder columns and remove UNK columns:
colOrder = ['COHORT', 'EUR', 'ASN', 'AFR', 'SAN', 'AMR'] # PCAWG

alleleCountDf = alleleCountDf[colOrder]

# Save output into tsv
outFilePath = outDir + '/' + fileName + '.alleleCount.tsv'
alleleCountDf.to_csv(outFilePath, sep='\t') 

### 4.2 MEI allele frequency 
# Create pandas dataframe from dictionary
alleleFreqDf = pd.DataFrame(alleleFreqDict) 

# transpose to have MEI as rows 
alleleFreqDf = alleleFreqDf.T 

# Reorder columns and remove UNK columns:
colOrder = ['COHORT', 'EUR', 'ASN', 'AFR', 'SAN', 'AMR']
alleleFreqDf = alleleFreqDf[colOrder]

# Save output into tsv
outFilePath = outDir + '/' + fileName + '.alleleFreq.tsv'
alleleFreqDf.to_csv(outFilePath, sep='\t') 

#### End
header("FINISH!!")


