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

## Get user's input ##
parser = argparse.ArgumentParser(description= "Make a table with the total number of MEI, number of L1, Alu and SVA per genotyped donor in a VCF")
parser.add_argument('inputVCF', help='Multi-sample VCF file containing genotyped MEI')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Read input multi-sample VCF and generate a VCF object
#############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Count the total number of MEI, number of L1, Alu and SVA per genotyped donor in the input VCF
#####################################################################################################

header("2. Count number of MEI per donor")

# Initialize dictionaries
nbTotalDict = {}
nbL1Dict = {}
nbAluDict = {}
nbSvaDict = {}

for donorId in donorIdList:

    nbTotalDict[donorId] = 0
    nbL1Dict[donorId] = 0
    nbAluDict[donorId] = 0
    nbSvaDict[donorId] = 0


# Iterate over the MEI in the VCF
for MEIObj  in VCFObj.lineList:
    MEItype = MEIObj.alt.split(":")[2][:-1]

    # Update counters on number of MEI per donor
    for donorId, genotype in MEIObj.genotypesDict.iteritems():

        # Individual bearing at least one copy of the current MEI
        if (genotype != "0") and (genotype != "0|0"):

            nbTotalDict[donorId] += 1

            # print "total: ", donorId, nbTotalDict[donorId]

            # a) LINE-1
            if (MEItype == "LINE1"):

                nbL1Dict[donorId] += 1
                # print "L1", MEItype, MEIObj.chrom, donorId, genotype, nbL1Dict[donorId]

            # b) ALU
            elif (MEItype == "ALU"):

                nbAluDict[donorId] += 1
                # print "ALU", MEItype, MEIObj.chrom, donorId, genotype, nbAluDict[donorId]

            # c) SVA
            else:

                nbSvaDict[donorId] += 1
                # print "SVA", MEItype, MEIObj.chrom, donorId, genotype, nbSvaDict[donorId]


#### 3. Make output table with number of MEI per donor
######################################################
header("3. Make output table with number of MEI per donor")

# Open output file
outFilePath = outDir + '/nbMEI_perDonor.txt'
outFile = open(outFilePath, 'a')

# Write header:
row = '#donorId' + "\t" + 'nbTotal' +  "\t" + 'nbL1' +  "\t" + 'nbAlu' +  "\t" + 'nbSva' + "\n"
outFile.write(row)

for donorId in donorIdList:

    nbTotal = nbTotalDict[donorId]
    nbL1 = nbL1Dict[donorId]
    nbAlu = nbAluDict[donorId]
    nbSva = nbSvaDict[donorId]

    row = donorId + "\t" + str(nbTotal) + "\t" + str(nbL1) +  "\t" + str(nbAlu) +  "\t" + str(nbSva) + "\n"
    outFile.write(row)

#### END
header("Finished")
