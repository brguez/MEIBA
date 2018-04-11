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

## Write file header in the output file
row = "#sampleId" + "\t" + "tumorType" + "\t" + "fractionPolyA" + "\n"     

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
        nbPolyA = 0
       
        # Read VCF and add information to VCF object
        VCFObj.read_VCF(vcfPath)
                
        ## For each MEI
        for MEIObj in VCFObj.lineList:  
  
            ## Select only those MEI that pass all the filters:
            if (MEIObj.filter == "PASS"):

                ## Update the counters 
                nbTotal += 1

                # Count the number of MEI with 3' bkp
                if (int(MEIObj.infoDict["SCORE"]) > 3):

                    nbPolyA += 1

        ## Compute the proportion of MEI with 3' bkp for the given sample
        if nbTotal == 0:
            fractionPolyA = 0.0
        else:

            fractionPolyA = float(nbPolyA) / nbTotal

        ## Write MEI counts into the output file
        row = sampleId + "\t" + tumorType + "\t" + str(fractionPolyA) + "\n"      
        outFile.write(row)

## End ##
print
print "***** Finished! *****"
print
