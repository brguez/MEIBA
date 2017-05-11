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


def genotypes2df(VCFObj):
    """
    """
    donorGtList = []

    ## For each MEI in the VCF
    for MEIObj in VCFObj.lineList:

        # Create a series of genotype (donorId labeled)
        end = (MEIObj.infoDict["BKPB"] if "BKPB" in MEIObj.infoDict else "UNK")
        sourceElementId = MEIObj.chrom + ':' + str(MEIObj.pos) + '-' + str(end)

        donorGt =  pd.Series(MEIObj.genotypesDict, name=sourceElementId)

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

    # A) Homozygous reference (for unknown genotypes con)
    if (genotype == '0') or (genotype == '0|0') or (genotype == '0/0') or (genotype == './.'):
        boolean = 0

    # B) Heterozygous or homozygous MEI (carrier/no_carrier)
    else:
        boolean = 1

    return boolean

def series2binary(integer):
    """
    """
    if (integer > 0):
        boolean = 1

    else:
        boolean = 0

    return boolean

def selectDonorSet(nbAbsentSrc, binaryGenotypes):
    """
    """
    nbDonors = binaryGenotypes.shape[1]
    nbSrcElements = binaryGenotypes.shape[0]

    percCovered = 0
    accumulatedSeries = pd.Series([0] * nbSrcElements, index=binaryGenotypes.index)
      
    for iteration in range(1, (nbDonors + 1)):
    
        print "** Iteration nb. ",  iteration, " **"   
        bestDonor, newPercCovered, accumulatedSeries = selectBestDonor(nbAbsentSrc, percCovered, accumulatedSeries, binaryGenotypes)

        # b)
        if (newPercCovered > percCovered):
            percCovered = newPercCovered
            selectedDonorList.append(bestDonor)

            ## Discard selected donor from the dataframe;
            binaryGenotypes = binaryGenotypes.drop(bestDonor, axis=1)

            print "tioo: ", percCovered, len(selectedDonorList), selectedDonorList

        # b)
        else:
            print "Stop! No percentage increase"            
    
            break

def selectBestDonor(nbAbsentSrc, percCovered, accumulatedSeries, binaryGenotypes):
    """
    """
    # A key per donorId. The value will be the percentage of source elements covered after adding the candidate donors to the list of selected donors. 
    percCoveredDict = {}

    # A key per donorId. The value will be a series containing for each source element its binary status (1: covered by the selected set of donors and 0: not covered )
    unionSeriesDict = {}
    
    for donorId in binaryGenotypes:
    
        candidateDonorSeries = binaryGenotypes[donorId]
        tmpSeries = accumulatedSeries.add(candidateDonorSeries)

        unionSeries = tmpSeries.apply(series2binary) 
        unionSeriesDict[donorId] = unionSeries

        nbCoveredCandidate = candidateDonorSeries.sum() # not needed
        nbCoveredAccumulated = unionSeries.sum()
        percCoveredDict[donorId] = float(nbCoveredAccumulated)/float(nbAbsentSrc)*100
            
    ## Select the donor contributing to the highest percentage of covered source elements
    # Note: if there are several possible donors, select one randomly
    bestDonor =  max(percCoveredDict, key=percCoveredDict.get)

    print "bestDonor: ", bestDonor, percCoveredDict[bestDonor]

    return(bestDonor, percCoveredDict[bestDonor], unionSeriesDict[bestDonor])
    


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('sourceElementGt', help='')
parser.add_argument('donorMetadata', help='')
parser.add_argument('sourceElementMetadata', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
sourceElementGt = args.sourceElementGt
donorMetadata = args.donorMetadata
sourceElementMetadata = args.sourceElementMetadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", sourceElementGt
print "donorMetadata: ", donorMetadata
print "sourceElementMetadata: ", sourceElementMetadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Read donor metadata file
#################################
# Initialize a dictionary with the following structure:
# - dict: key(donorId) -> projectCode

header("1. Read donor metadata file")
metadataFile = open(donorMetadata, 'r')

donorIdProjectCodeDict = {}

for line in metadataFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')

        donorId = line[0]
        tumorType = line[9]
        donorIdProjectCodeDict[donorId] = tumorType

        #print "test: ", donorId, tumorType

#print "donorIdProjectCodeDict: ", donorIdProjectCodeDict

#### 2. Compute the allele count of each source element in EOPC-DE
###################################################################
## EOPC-DE is the tumor type with available samples for the validation of L1 source elements. 
# Initialize a dictionary with the following structure:
# - dict1: key(sourceElementId) -> dict2: key1("alleleCount") -> value1(alleleCount value)
#                                         key2("donorIdList") -> list of donor ids containing the insertion
# sourceElementId: chr:beg-end

header("2. Compute the allele count of each source element in EOPC-DE")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(sourceElementGt)

alleleCountsDict = {}


## For each MEI:
for MEIObj in VCFObj.lineList:

    end = (MEIObj.infoDict["BKPB"] if "BKPB" in MEIObj.infoDict else "UNK")

    sourceElementId = MEIObj.chrom + ':' + str(MEIObj.pos) + '-' + str(end)
    print "** source element ** ", sourceElementId

    ## Initialize source element dictionary
    alleleCountsDict[sourceElementId] = {}
    alleleCountsDict[sourceElementId]["alleleCount"] = 0
    alleleCountsDict[sourceElementId]["donorIdList"] = []


    ## For each donor:
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        #print donorId

        ## Project code available for current donors
        if (donorId in donorIdProjectCodeDict):
            projectCode = donorIdProjectCodeDict[donorId]

            # Select EOPC-DE tumor types
            if (projectCode == "EOPC-DE"):

                #print "insertion-Gt: ", genotype
                ## A) Insertion absent in reference genome
                if (MEIObj.alt == "<MEI>"):
            
                    # a) Heterozygous
                    if (genotype == "0/1"):
                        alleleCountsDict[sourceElementId]["alleleCount"] += 1   
                        alleleCountsDict[sourceElementId]["donorIdList"].append(donorId)

                    # b) Homozygous alternative
                    elif (genotype == "1/1"):
                        alleleCountsDict[sourceElementId]["alleleCount"] += 2   
                        alleleCountsDict[sourceElementId]["donorIdList"].append(donorId)

                    # Note c) possibility would be missing allele (./.)

                    #print "Insertion absent in reference genome", sourceElementId, donorId, projectCode, alleleCountsDict[sourceElementId]      

                ## B) Insertion in reference genome and absent in donor genome
                elif (MEIObj.ref == "<MEI>"):
        
                    #print "Insertion in reference genome", donorId, genotype, projectCode        

                    # a) Heterozygous
                    if (genotype == "0/1"):
                
                        alleleCountsDict[sourceElementId]["alleleCount"] += 1
                        alleleCountsDict[sourceElementId]["donorIdList"].append(donorId)

                    # b) Homozygous reference
                    elif (genotype == "0/0"):
                                      
                        alleleCountsDict[sourceElementId]["alleleCount"] += 2
                        alleleCountsDict[sourceElementId]["donorIdList"].append(donorId)

           
        # b) Project code not available for current donors (affects to few donors) 
        # I don't know why this happens...  check later     
        else:
         
            print "[ERROR] Unknown donor tumor type: ", donorId       


#### 3. Make output file containing source element metadata + allele count
##############################################################################
header("3. Make output file containing source element metadata + allele count")
metadataFile = open(sourceElementMetadata, 'r')

# Open output file
outFilePath = outDir + '/sourceElements_alleleCount_EOPCDE.tsv'
outFile = open(outFilePath, 'w')

# Write header:
row = '#cytobandId' + "\t" + 'sourceIdNew' + "\t" + 'sourceIdOld' + "\t" + 'refStatus' + "\t" +  'strand' + "\t" + 'TSDlen' + "\t" + 'novelty' + "\t" + 'activityStatus' + "\t" + 'alleleCount' + "\t" + 'donorIdList' + "\n"
outFile.write(row)

for line in metadataFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')
    
        print "line: ", line
        cytobandId, sourceIdNew, sourceIdOld, novelty, activityStatus = line

        ## If inconsistencies between source element identifier raise an error an skip
        # Problem only affects one element
        if sourceIdNew not in alleleCountsDict:
            continue
            print "[ERROR] source element coordenate not found "

        alleleCount = alleleCountsDict[sourceIdNew]["alleleCount"] 
        
        donorIdList = (','.join(alleleCountsDict[sourceIdNew]["donorIdList"]) if len(alleleCountsDict[sourceIdNew]["donorIdList"]) > 0 else "-")

        row = cytobandId + "\t" + sourceIdNew + "\t" + sourceIdOld + "\t" + novelty + "\t" + activityStatus + "\t" + str(alleleCount) + "\t" + donorIdList + "\n"
        outFile.write(row)


#### 4. Make genotyping binary matrix for EOPC-DE donors
##########################################################
# The binary matrix will only contain those source elements
# that are absent in the reference genome
## carrier: 1, no_carrier: 0
header("4. Make genotyping binary matrix for EOPC-DE donors")

#### Select MEI absent in the reference genome
VCFAbsentObj = formats.VCF()

## For each MEI:
for MEIObj in VCFObj.lineList:
    
    if (MEIObj.alt == "<MEI>"):
        VCFAbsentObj.addLine(MEIObj)

#### Make binary matrix for all the donors
gtAbsentDfPCAWG = genotypes2df(VCFAbsentObj)

gtAbsentBinaryDfPCAWG = gtAbsentDfPCAWG.applymap(gt2binary)

#### Filter binary matrix selecting EOPC-DE donors
# print "gtBinaryDfPCAWG: ", gtBinaryDfPCAWG

## Make list with EOPC donor ids
EOPCdonorIdList = []

for donorId in donorIdProjectCodeDict:
    
    projectCode = donorIdProjectCodeDict[donorId]

    # Select EOPC-DE tumor types
    if (projectCode == "EOPC-DE"):
        EOPCdonorIdList.append(donorId)


binaryGtAbsentSrcEOPCdf = gtAbsentBinaryDfPCAWG[EOPCdonorIdList]


#### 5. Find the collection of donors maximizing the number of 
################################################################
# source elements absent in the reference genome
##################################################
header("5. Find the collection of donors maximizing the number of source elements absent in the reference genome")

nbAbsentSrc = len(VCFAbsentObj.lineList)
selectedDonorList = []

print "nbAbsentSrc: ", nbAbsentSrc

selectDonorSet(nbAbsentSrc, binaryGtAbsentSrcEOPCdf)


#### 6. Report the number of source L1 in each EOPC 
####################################################

# Open output file
outFilePath = outDir + '/donorId_nbSourceL1_EOPCDE.tsv'
outFile = open(outFilePath, 'w')

for donorId in binaryGtAbsentSrcEOPCdf:

    sourceL1list = [ sourceId for sourceId, active in binaryGtAbsentSrcEOPCdf[donorId].iteritems() if active == 1]
    nbSourceL1 = len(sourceL1list)
    sourceL1Str = ",".join(sourceL1list)
    
    row = donorId + "\t" + str(nbSourceL1) + "\t" + sourceL1Str + "\n"
    outFile.write(row)

#    medianNVHom = np.median([float(genotype[1]) for genotype in genotypesList if genotype[0] == '1/1'])



#binaryGtAbsentSrcEOPCdf = binaryGtAbsentSrcEOPCdf

#print "binaryGtAbsentSrcEOPCdf: ", binaryGtAbsentSrcEOPCdf

####
header("Finished")
