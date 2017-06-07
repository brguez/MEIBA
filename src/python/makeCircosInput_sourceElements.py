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
parser.add_argument('VCF', help='VCF with somatic MEI calls for a given sample')
parser.add_argument('sourceIds', help='List of source element identifiers to include transductions in the output file')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
VCF = args.VCF
sourceIds = args.sourceIds
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "sourceIds: ", sourceIds
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## PROCESS VCF FILES
#####################
outPath = outDir + '/' + fileName + '.tsv'

outFile = open(outPath, 'w')

# Source Ids list
sourceIdsList = sourceIds.split(',')

# Create VCF object
VCFObj = formats.VCF()

## Raise error if the input VCF is not available
if not os.path.isfile(VCF):

    print "[ERROR] Input file is not available"

## Input VCF available
else:

    # Read VCF and add information to VCF object
    VCFObj.read_VCF(VCF)
        
    tdCounter = 1

    ## For each MEI
    for MEIObj in VCFObj.lineList:  
            
        ## Select only those MEI that pass all the filters:
        if (MEIObj.filter == "PASS"):

            if ('SRCID' in MEIObj.infoDict) and (MEIObj.infoDict['SRCID'] in sourceIdsList):
                
                tdId = MEIObj.infoDict['SRCID'] + "_" + str(tdCounter)

                insertionChrom = 'hs' + MEIObj.chrom
                insertionBeg = str(MEIObj.pos)
                insertionEnd = str(MEIObj.pos + 1)
                
                sourceCoordList = MEIObj.infoDict['SRC'].split('_')
                sourceChrom = 'hs' + sourceCoordList[0]
                sourceBeg = sourceCoordList[1]
                sourceEnd = sourceCoordList[2]

    
                print tdId, insertionChrom, insertionBeg, insertionEnd
                print tdId, sourceChrom, sourceBeg, sourceEnd

                tdCounter += 1

                ## Write transduction coordinates into the output file
                row = tdId + '\t' + insertionChrom + '\t' + insertionBeg + '\t' + insertionEnd + '\n'
                outFile.write(row)

                row = tdId + '\t' + sourceChrom + '\t' + sourceBeg + '\t' + sourceEnd + '\n'
                outFile.write(row)
 

## End ##
print
print "***** Finished! *****"
print
