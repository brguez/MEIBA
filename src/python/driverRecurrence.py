#!/usr/bin/env python
#coding: utf-8

def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"

def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
from operator import itemgetter, attrgetter, methodcaller
import pandas as pd

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: donorId vcfPath')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Count the number of MEI affecting each gene
###############################################

inputFile = open(inputPath, 'r')

totalNbL1 = 0 # Total number L1 insertions
totalNbExonic = 0 # Total number L1 insertions in exons
totalNbGene = 0 # Total number L1 insertions affecting genes (exonic, splicing, UTR5, UTR3 or intronic regions)

recurrencyExonicDict = {}
recurrencyGeneDict = {}

outFilePath = outDir + '/L1_events_PCAWG.bed'
outFile = open(outFilePath, 'w')

# Per iteration, read a VCF
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")

    donorId = line[0]
    vcfPath = line[1]

    #print "Processing: ", donorId, vcfPath

    # Create VCF object
    VCFObj = formats.VCF()

    ## Raise error if the input VCF is not available
    if not os.path.isfile(vcfPath):

        print "[ERROR] Input file is not available"

    ## Input VCF available
    else:
     
        # Read VCF and add information to VCF object
        VCFObj.read_VCF(vcfPath)
        
        ## For each MEI
        for MEIObj in VCFObj.lineList:  

            MEIClass = MEIObj.infoDict["CLASS"] if "CLASS" in MEIObj.infoDict else "NA"

            ## Select only those L1 insertions that that pass all the filters:
            if (MEIObj.filter == "PASS") and (MEIClass == "L1"):
                totalNbL1 += 1

                row = MEIObj.chrom + "\t" + str(MEIObj.pos) + "\t" + str(MEIObj.pos + 1) + "\n"
                outFile.write(row)


                # Select insertions affecting genes
                if ("GENE" in MEIObj.infoDict): 

                    geneName = MEIObj.infoDict["GENE"]
                    MEIRole = MEIObj.infoDict["ROLE"] if "ROLE" in MEIObj.infoDict else "-"    

                    # a) Exonic MEI
                    if (MEIObj.infoDict["REGION"] == "exonic"):

                        totalNbExonic += 1     

                        if (geneName not in recurrencyExonicDict):
                            recurrencyExonicDict[geneName] = {}
                            recurrencyExonicDict[geneName]['gnRole'] = MEIRole
                            recurrencyExonicDict[geneName]['nbMEI'] = 1
                        else:
                            recurrencyExonicDict[geneName]['nbMEI'] += 1

                    # b) MEI within gene boundaries
                    elif (MEIObj.infoDict["REGION"] in ['exonic', 'splicing', 'UTR5', 'UTR3', 'intronic']):

                        totalNbGene += 1

                        if (geneName not in recurrencyGeneDict):
                            recurrencyGeneDict[geneName] = {}
                            recurrencyGeneDict[geneName]['gnRole'] = MEIRole
                            recurrencyGeneDict[geneName]['nbMEI'] = 1
                        else:
                            recurrencyGeneDict[geneName]['nbMEI'] += 1

## Convert into dataframe
###########################
## Exonic MEI
recurrencyExonicDf = pd.DataFrame(recurrencyExonicDict) 
recurrencyExonicDf = recurrencyExonicDf.T

recurrencyExonicSortedDf = recurrencyExonicDf.sort_values('nbMEI', ascending=False) # sort in decreasing order

## MEI within gene boundaries
recurrencyGeneDf = pd.DataFrame(recurrencyGeneDict) 
recurrencyGeneDf = recurrencyGeneDf.T
recurrencyGeneSortedDf = recurrencyGeneDf.sort_values('nbMEI', ascending=False) # sort in decreasing order


## Print dataframes as tsv
###########################
## Exonic MEI
outFilePath = outDir + '/recurrency_exonic.tsv'
recurrencyExonicSortedDf.to_csv(outFilePath, sep='\t') 

## MEI within gene boundaries
outFilePath = outDir + '/recurrency_gene.tsv'
recurrencyGeneSortedDf.to_csv(outFilePath, sep='\t') 


## Report basic stats
######################
print "*** Summary ***"
print "totalNbL1: ", totalNbL1
print "totalNbExonic: ", totalNbExonic
print "totalNbGeneBoundaries: ", totalNbGene

