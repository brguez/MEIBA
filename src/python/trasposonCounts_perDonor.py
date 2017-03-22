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
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: tumorId   donorId   projectCode   histology vcfPath')
parser.add_argument('delPath', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
delPath = args.delPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "delPath: ", delPath
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## PROCESS DELETION COUNTS
############################
delFile = open(delPath, 'r')

delDict = {}

# Per iteration, read a VCF
for line in delFile:
    line = line.rstrip('\n')
    line = line.split("\t")
    donorId = line[0]
    projectCode = line[1]
    donorUniqueId = projectCode + "::" + donorId 
    nbDel = line[2]

    delDict[donorUniqueId] = nbDel

print "delDict: ", delDict

## PROCESS VCF FILES
#####################
outPath = outDir + '/transposonCounts_perDonor.tsv'

outFile = open(outPath, 'w')

## Write file header in the output file
row = "#donorId" + "\t" + "projectCode" + "\t" + "histology" + "\t" + "nbTotal" + "\t" + "nbSoloL1" + "\t" + "nbTD" + "\t" + "nbAlu" + "\t" + "nbSVA" + "\t" + "nbERVK"  + "\t" + "nbPSD" + "\t" + "nbDEL" + "\n"     

outFile.write(row)

inputFile = open(inputPath, 'r')

# Per iteration, read a VCF
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")

    donorId = line[0]
    projectCode = line[1]
    donorUniqueId = projectCode + "::" + donorId 
    projectCode = projectCode.split("-")[0]
    histology = line[2]
    vcfPath = line[3]

    print "Processing: ", donorId, projectCode, histology, vcfPath

    # Create VCF object
    VCFObj = formats.VCF()

    ## Raise error if the input VCF is not available
    if not os.path.isfile(vcfPath):

        print "[ERROR] Input file is not available"

    ## Input VCF available
    else:
        # Initialize counters:
        nbTotal = 0       
        nbPSD = 0
        nbSoloL1 = 0
        nbTD = 0
        nbAlu = 0
        nbSVA = 0
        nbERVK = 0
        nbDEL = delDict[donorUniqueId] if (donorUniqueId in delDict) else 0

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

                # B) L1 transduction (orphan or partnered)
                elif (MEItype == "TD1") or (MEItype == "TD2"):
                    nbTD += 1

                # C) Solo integration (L1, Alu, SVA or ERVK)
                elif (MEItype == "TD0"):
    
                    MEIClass = MEIObj.infoDict["CLASS"]

                    ## a) L1:
                    if (MEIClass == "L1"):
                        nbSoloL1 += 1

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
        row = donorId + "\t" + projectCode + "\t" + histology + "\t" + str(nbTotal) + "\t" + str(nbSoloL1) + "\t" + str(nbTD) + "\t" + str(nbAlu) + "\t" + str(nbSVA) + "\t" + str(nbERVK)  + "\t" + str(nbPSD) + "\t" + str(nbDEL) + "\n"      
        outFile.write(row)

## End ##
print
print "***** Finished! *****"
print
