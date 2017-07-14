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
parser.add_argument('inputPath', help='Tabular text file containing one row per donor with the following consecutive fields: donorId   tumorType vcfPath')
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
outPath = outDir + fileName + '.tsv'

outFile = open(outPath, 'w')

## Write file header in the output file
row = "#donorId" + "\t" + "tumorType" + "\t" + "nbTotal" + "\t" + "nbSoloL1" + "\t" + "nbL1TD" + "\t" + "nbAlu" + "\t" + "nbSVA" + "\t" + "nbERVK"  + "\t" + "nbPSD" + "\t" + "nbL1DEL" + "\t" + "nbL1DUP" + "\n"     

outFile.write(row)

inputFile = open(inputPath, 'r')

# Per iteration, read a VCF
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")

    donorId = line[0]  
    tumorType = line[1]
    vcfPath = line[2]

    print "Processing: ", donorId, tumorType, vcfPath

    # Create VCF object
    VCFObj = formats.VCF()

    ## Raise error if the input VCF is not available
    if not os.path.isfile(vcfPath):

        print "[ERROR] Input file is not available"

    ## Input VCF available
    else:
        # Initialize counters:
        nbTotal = 0       
        nbSoloL1 = 0
        nbL1TD = 0
        nbAlu = 0
        nbSVA = 0
        nbERVK = 0
        nbPSD = 0
        nbL1DEL = 0
        nbL1DUP = 0

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
                    # print "PSEUDOGENE: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict   
       
                # B) L1 mediated deletion
                elif ("GR" in MEIObj.infoDict) and (MEIObj.infoDict["GR"] == "DEL"):
                    nbL1DEL += 1
                    # print "DELETION: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict   

                # B) L1 mediated duplication
                elif ("GR" in MEIObj.infoDict) and (MEIObj.infoDict["GR"] == "DUP"):
                    nbL1DUP += 1
                    # print "DUPLICATION: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict   

                # B) L1 transduction (orphan or partnered)
                elif ("GR" not in MEIObj.infoDict) and ((MEItype == "TD1") or (MEItype == "TD2")):
                    nbL1TD += 1
                    #print "TRANSDUCTION: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict   

                # C) Solo integration (L1, Alu, SVA or ERVK)
                elif ("GR" not in MEIObj.infoDict) and (MEItype == "TD0"):
    
                    MEIClass = MEIObj.infoDict["CLASS"]

                    ## a) L1:
                    if (MEIClass == "L1"):
                        nbSoloL1 += 1
                        #print "SOLO-L1: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict   

                    ## b) Alu
                    elif (MEIClass == "Alu"):
                        nbAlu += 1
                        #print "SOLO-Alu: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict  
               
                    ## c) SVA
                    elif (MEIClass == "SVA"):
                        nbSVA += 1
                        #print "SOLO-SVA: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict  

                    ## e) ERVK
                    elif (MEIClass == "ERVK"):
                        nbERVK += 1 

        ## Write MEI counts into the output file
        row = donorId + "\t" + tumorType + "\t" + str(nbTotal) + "\t" + str(nbSoloL1) + "\t" + str(nbL1TD) + "\t" + str(nbAlu) + "\t" + str(nbSVA) + "\t" + str(nbERVK)  + "\t" + str(nbPSD) + "\t" + str(nbL1DEL) + "\t" + str(nbL1DUP) + "\n"      
        outFile.write(row)

## End ##
print
print "***** Finished! *****"
print
