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

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: sampleId   tumorType vcfPath')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print




## PROCESS VCF FILES
#####################
outPath = outDir + '/' + fileName + '.tsv'

outFile = open(outPath, 'w')

## Write file header in the output file (TO DO)
row = "#sampleId" + "\t" + "tumorType" + "\t" + "nbTotal" + "\t" + "nbL1" + "\t" + "nbL1Solo" + "\t" + "nbL1TD" + "\t" + "nbL1DEL" + "\t" + "nbL1DUP" + "\t" + "nbAlu" + "\t" + "nbSVA" + "\t" + "nbERVK"  + "\t" + "nbPSD" + "\n"

outFile.write(row)

inputFile = open(inputPath, 'r')

# Per iteration, read a VCF
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")

    sampleId = line[0]  
    tumorType = line[1]
    vcfPath = line[2]

    print "Processing: ", sampleId, tumorType, vcfPath

    # Create VCF object
    VCFObj = formats.VCF()

    ## Raise error if the input VCF is not available
    if not os.path.isfile(vcfPath):

        print "[ERROR] Input file is not available"

    ## Input VCF available
    else:
        # Initialize counters:
        nbTotal = 0   
        nbL1 = 0    
        nbL1Solo = 0
        nbL1TD = 0
        nbL1DEL = 0
        nbL1DUP = 0
        nbAlu = 0
        nbSVA = 0
        nbERVK = 0
        nbPSD = 0

        # Read VCF and add information to VCF object
        VCFObj.read_VCF(vcfPath)
        
        ## For each MEI
        for MEIObj in VCFObj.lineList:  
  
            ## Select only those MEI that pass all the filters:
            if (MEIObj.filter == "PASS"):

                ## Update the counters 
                nbTotal += 1
           
                MEItype = MEIObj.infoDict["TYPE"] 

                # A) Processed pseudogene (PSD)
                if (MEItype == "PSD"):
                    nbPSD += 1

                # B) L1 mediated deletion
                elif ("GR" in MEIObj.infoDict) and (MEIObj.infoDict["GR"] == "DEL"):
                    nbL1DEL += 1

                # B) L1 mediated duplication
                elif ("GR" in MEIObj.infoDict) and (MEIObj.infoDict["GR"] == "DUP"):
                    nbL1DUP += 1

                # B) L1 transduction (orphan or partnered)
                elif ("GR" not in MEIObj.infoDict) and ((MEItype == "TD1") or (MEItype == "TD2")):
                    nbL1TD += 1
                    nbL1 += 1

                # C) Solo integration (L1, Alu, SVA or ERVK)
                elif ("GR" not in MEIObj.infoDict) and (MEItype == "TD0"):
    
                    MEIClass = MEIObj.infoDict["CLASS"]

                    ## a) L1:
                    if (MEIClass == "L1"):
                        nbL1Solo += 1
                        nbL1 += 1

                    ## b) Alu
                    elif (MEIClass == "Alu"):
                        nbAlu += 1
               
                    ## c) SVA
                    elif (MEIClass == "SVA"):
                        nbSVA += 1

                    ## e) ERVK
                    elif (MEIClass == "ERVK"):
                        nbERVK += 1 

        ## Write MEI counts into the output file
        row = sampleId + "\t" + tumorType + "\t" + str(nbTotal) + "\t" + str(nbL1) + "\t" + str(nbL1Solo) + "\t" + str(nbL1TD) + "\t" + str(nbL1DEL) + "\t" + str(nbL1DUP) + "\t" + str(nbAlu) + "\t" + str(nbSVA) + "\t" + str(nbERVK)  + "\t" + str(nbPSD)  + "\n"      
        outFile.write(row)

## End ##
print
print "***** Finished! *****"
print
