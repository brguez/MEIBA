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
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: tumorSpecimenId   icgcDonorId vcfPath')
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

## PROCESS VCF FILES
#####################
outPath = outDir + '/PCAWG_somaticRetrotransposonCounts_perSpecimen_0.5.7.tsv'

outFile = open(outPath, 'w')

         
## Write file header in the output file
row = "icgc_donor_id" + "\t" + "dcc_project_code" + "\t" + "totalNbMEI" + "\t" + "nbL1" + "\t" + "nbAlu" + "\t" + "nbSVA" + "\t" + "nbERVK" + "\n"     

outFile.write(row)

inputFile = open(inputPath, 'r')

# Per iteration, read a VCF
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")

    tumorSpecimenId = line[0]
    icgcDonorId = line[1]
    tumorType = line[2]
    vcfPath = line[3]

    print "Processing: ", tumorSpecimenId, icgcDonorId, tumorType, vcfPath

    # Create VCF object
    VCFObj = formats.VCF()

    ## Raise error if the input VCF is not available
    if not os.path.isfile(vcfPath):

        print "[ERROR] Input file is not available"

    ## Input VCF available
    else:
        # Initialize counters:
        totalNbMEI = 0       
        nbL1 = 0
        nbAlu = 0
        nbSVA = 0
        nbERVK = 0
        
        # Read VCF and add information to VCF object
        VCFObj.read_VCF(vcfPath)
        
        ## For each MEI
        for MEIObj in VCFObj.lineList:  
            
            ## Select only those MEI that pass all the filters:
            if (MEIObj.filter == "PASS"):

                ## Update the counters 
                totalNbMEI += 1
                MEIClass = MEIObj.infoDict["CLASS"]

                ## a) L1:
                if (MEIClass == "L1"):
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
        row = tumorSpecimenId + "\t" + icgcDonorId + "\t" + tumorType + "\t" + str(totalNbMEI) + "\t" + str(nbL1) + "\t" + str(nbAlu) + "\t" + str(nbSVA) + "\t" + str(nbERVK) + "\n"      
        outFile.write(row)

## End ##
print
print "***** Finished! *****"
print
