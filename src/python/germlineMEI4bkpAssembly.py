#!/usr/bin/env python
#coding: utf-8

## Load modules/libraries
import sys
import argparse
import os
import errno

## Get user's input
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('insertions', help='')
parser.add_argument('sampleId', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
insertionsPath = args.insertions
normal_wgs_aliquot_id = args.sampleId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "insertions: ", insertionsPath
print "sampleId: ", normal_wgs_aliquot_id 
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print

## Create output file with insertions:
outPath = outDir + "/" + normal_wgs_aliquot_id + ".germline.td0.tsv"
outFile = open(outPath, "w" )

# Print header into the output file
header = "#chromPlus" + "\t" + "begPlus" + "\t" + "endPlus" + "\t" + "nbReadsPlus" + "\t" + "classPlus" + "\t" + "readListPlus" + "\t" + "chromMinus" + "\t" + "begMinus" + "\t" + "endMinus" + "\t" + "nbReadsMinus" + "\t" + "classMinus" + "\t" + "readListMinus" + "\t" + "insertionType" + "\t" + "chromSource" + "\t" + "begSource" + "\t" + "endSource" + "\t" + "strandSource" + "\t" + "tdBeg" + "\t" + "tdEnd" + "\t" + "tdRnaLen" + "\t" + "tdLen" + "\t" + "psdGene" + "\t" + "chromExonA" + "\t" + "begExonA" + "\t" + "endExonA" + "\t" + "chromExonB" + "\t" + "begExonB" + "\t" + "endExonB" + "\n"
outFile.write(header)

### Read insertions file 
insertions = open(insertionsPath, 'r')

# Read file line by line
for line in insertions:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
       
        chromPlus = fieldsList[0]
        begPlus = fieldsList[1] 
        endPlus = fieldsList[2]  
        nbReadsPlus = fieldsList[3]
        classPlus = fieldsList[4]
        readListPlus = fieldsList[5]
        chromMinus = fieldsList[6]
        begMinus = fieldsList[7] 
        endMinus = fieldsList[8]   
        nbReadsMinus = fieldsList[9] 
        classMinus = fieldsList[10]
        readListMinus = fieldsList[11]
            
        insertionType = "TD0" 
         
        ## Set transduction specific fields as not applicable (NA)
        chromSource = "NA" 
        begSource = "NA"
        endSource = "NA"
        strandSource = "NA"
        tdBeg = "NA"
        tdEnd = "NA"
        tdRnaLen = "NA"
        tdLen = "NA"

        ## Set pseudogene specific fields as not applicable (NA)
        psdGene = "NA"
        chromExonA = "NA"        
        begExonA = "NA"
        endExonA = "NA"
        chromExonB = "NA"
        begExonB = "NA"
        endExonB = "NA"

        # Print insertion information into the output file  
        row = chromPlus + "\t" + begPlus + "\t" + endPlus + "\t" + nbReadsPlus + "\t" + classPlus + "\t" + readListPlus + "\t" + chromMinus + "\t" + begMinus + "\t" + endMinus + "\t" + nbReadsMinus + "\t" + classMinus + "\t" + readListMinus + "\t" + insertionType + "\t" + chromSource + "\t" + begSource + "\t" + endSource + "\t" + strandSource + "\t" + tdBeg + "\t" + tdEnd + "\t" + tdRnaLen + "\t" + tdLen + "\t" + psdGene + "\t" + chromExonA + "\t" + begExonA + "\t" + endExonA + "\t" + chromExonB + "\t" + begExonB + "\t" + endExonB + "\n"

        outFile.write(row) 

print "***** Finished! *****"
print

