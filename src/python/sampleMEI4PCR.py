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

def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string


#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import os.path
import os.path
import random


## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('input', help='')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputFile = args.input
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputFile: ", inputFile
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 


## Organize insertions passing the filters into a dictionary
#############################################################
inputFile = open(inputFile, 'r')

### Initialize dictionaries
insertionsDict = {}
insertionsDict["L1"] = {}
insertionsDict["L1"]["TD0"] = []
insertionsDict["L1"]["TD1"] = []
insertionsDict["L1"]["TD2"] = []
insertionsDict["L1"]["EI"] = []

insertionsDict["Alu"] = {}
insertionsDict["Alu"]["TD0"] = []

# Read file line by line
for line in inputFile:
    line = line.rstrip('\r\n')
    
    ## Discard header
    if not line.startswith("#"):        
        fieldsList = line.split("\t")

        rtClass = fieldsList[5]
        insertionType = fieldsList[6]
        mechanism = fieldsList[7]
        filterField = fieldsList[8]

        ## Select insertions passing the filters:
        if (filterField == "PASS"):

            ## Add insertion to the list
            insertionsDict[rtClass][insertionType].append(line)

            ## Add solo EI insertions to an independent list
            if (mechanism == "EI") and (insertionType == "TD0"):
                insertionsDict[rtClass]["EI"].append(line)

#print "insertionsDict: ", len(insertionsDict["L1"]["TD0"]), len(insertionsDict["L1"]["TD1"]), len(insertionsDict["L1"]["TD2"]), len(insertionsDict["L1"]["EI"]), len(insertionsDict["Alu"]["TD0"])


## Pick random samples to do PCR
#################################
## Sample 30 of each category. For solo EI L1 and Alu select all of them as there are few events (<10)
sampleSize = 30

## Solo L1
nbSolo = len(insertionsDict["L1"]["TD0"])
finalSize = nbSolo if (sampleSize > nbSolo) else sampleSize ## if less than sample size use the maximum possible
indexList = random.sample(xrange(0, nbSolo), finalSize)
soloL1List = list(insertionsDict["L1"]["TD0"][index] for index in indexList)

## Partnered L1
nbPartnered = len(insertionsDict["L1"]["TD1"])
finalSize = nbPartnered if (sampleSize > nbPartnered) else sampleSize ## if less than sample size use the maximum possible
indexList = random.sample(xrange(0, nbPartnered), finalSize)
partneredL1List = list(insertionsDict["L1"]["TD1"][index] for index in indexList)

## Orphan L1
nbOrphan = len(insertionsDict["L1"]["TD2"])
finalSize = nbOrphan if (sampleSize > nbOrphan) else sampleSize ## if less than sample size use the maximum possible
indexList = random.sample(xrange(0, nbOrphan), finalSize)
orphanL1List = list(insertionsDict["L1"]["TD2"][index] for index in indexList)

## Solo EI L1 
EIL1List = insertionsDict["L1"]["EI"]

## Solo Alu
aluList = insertionsDict["Alu"]["TD0"]


## Write output files
#######################

## Solo L1
filePath = outDir + '/' + fileName + ".solo_L1.tsv"
outFile = open(filePath, "w")
row = "#chrom" + "\t" + "bkp5prime" + "\t" + "contig5prime" + "\t" + "bkp3prime" + "\t" + "contig3prime" + "\t" + "rtClass" + "\t" + "insertionType" + "\t" + "mechanism" + "\t" + "filter" + "\t" + "score" + "\t" + "tsLen" + "\t" + "strand" + "\t" + "structure" + "\t" + "length" + "\n"
outFile.write(row)
outFile.write("\n".join(map(lambda x: str(x), soloL1List)))
outFile.write("\n")
outFile.close()

## Partnered L1
filePath = outDir + '/' + fileName + ".partnered_L1.tsv"
outFile = open(filePath, "w")
row = "#chrom" + "\t" + "bkp5prime" + "\t" + "contig5prime" + "\t" + "bkp3prime" + "\t" + "contig3prime" + "\t" + "rtClass" + "\t" + "insertionType" + "\t" + "mechanism" + "\t" + "filter" + "\t" + "score" + "\t" + "tsLen" + "\t" + "strand" + "\t" + "structure" + "\t" + "length" + "\n"
outFile.write(row)
outFile.write("\n".join(map(lambda x: str(x), partneredL1List)))
outFile.write("\n")
outFile.close()

## Orphan L1
filePath = outDir + '/' + fileName + ".orphan_L1.tsv"
outFile = open(filePath, "w")
row = "#chrom" + "\t" + "bkp5prime" + "\t" + "contig5prime" + "\t" + "bkp3prime" + "\t" + "contig3prime" + "\t" + "rtClass" + "\t" + "insertionType" + "\t" + "mechanism" + "\t" + "filter" + "\t" + "score" + "\t" + "tsLen" + "\t" + "strand" + "\t" + "structure" + "\t" + "length" + "\n"
outFile.write(row)
outFile.write("\n".join(map(lambda x: str(x), orphanL1List)))
outFile.write("\n")
outFile.close()

## Solo EI L1 
filePath = outDir + '/' + fileName + ".EI_L1.tsv"
outFile = open(filePath, "w")
row = "#chrom" + "\t" + "bkp5prime" + "\t" + "contig5prime" + "\t" + "bkp3prime" + "\t" + "contig3prime" + "\t" + "rtClass" + "\t" + "insertionType" + "\t" + "mechanism" + "\t" + "filter" + "\t" + "score" + "\t" + "tsLen" + "\t" + "strand" + "\t" + "structure" + "\t" + "length" + "\n"
outFile.write(row)
outFile.write("\n".join(map(lambda x: str(x), EIL1List)))
outFile.write("\n")
outFile.close()

## Solo Alu
filePath = outDir + '/' + fileName + ".solo_Alu.tsv"
outFile = open(filePath, "w")
row = "#chrom" + "\t" + "bkp5prime" + "\t" + "contig5prime" + "\t" + "bkp3prime" + "\t" + "contig3prime" + "\t" + "rtClass" + "\t" + "insertionType" + "\t" + "mechanism" + "\t" + "filter" + "\t" + "score" + "\t" + "tsLen" + "\t" + "strand" + "\t" + "structure" + "\t" + "length" + "\n"
outFile.write(row)
outFile.write("\n".join(map(lambda x: str(x), aluList)))
outFile.write("\n")
outFile.close()

## Finish ##
print
print "***** Finished! *****"
print




