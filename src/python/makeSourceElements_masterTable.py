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
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


## This script takes as input a set of files contatining information about germline source elements and
# combines the info into a single table with the following fields:

# 1. source element id (cytoband labelled)
# 2. source element new coordinates
# 3. source element old coordinates
# 4. source element "present" or "absent" in the reference genome
# 5. source element strand
# 6. source element TSD length
# 7. source element novelty (novel, known)
# 8. source element activity status (none, low, moderate and high)

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('masterTable', help='')
parser.add_argument('status', help='')
parser.add_argument('known', help='')
parser.add_argument('sourceElementGt', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
masterTable = args.masterTable
status = args.status
known = args.known
sourceElementGt = args.sourceElementGt
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "masterTable: ", masterTable
print "status: ", status
print "known: ", known
print "inputVCF: ", sourceElementGt
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read source element activity file
##########################################
# Initialize a dictionary with the following structure:
# - dict: key(sourceIdNew) -> activityStatus
header("1. Read activity status file")
statusFile = open(status, 'r')

statusDict = {}

for line in statusFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')

        sourceIdNew = line[0]
        status = line[4]
        statusDict[sourceIdNew] = status

#### 2. Read source element novelty file
##########################################
# Initialize a dictionary with the following structure:
# - dict: key(chrom) -> [coordSrcX, coordSrcY, ... ]
header("2. Read source element novelty file")
knownFile = open(known, 'r')

knownDict = {}

for line in knownFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')
        chrom, pos = line

        if chrom not in knownDict:
            knownDict[chrom] = []
        
        knownDict[chrom].append(pos)        
      

#### 3. Read source element VCF to extract src elements strand, TSD length and 
##############################################################################
# if the element is absent or present in the reference genome
##############################################################
# Initialize a dictionary with the following structure:
# - dict1: key(sourceElementId) -> dict2: key1("strand") -> source element strand or orientation
#                                         key2("TSDlen") -> TSD length
#                                         key3("REF") -> "present" or "absent" in the reference genome
# sourceElementId: chr:beg-end

header("3. Read source element VCF to extract src elements strand and TSD length")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(sourceElementGt)

srcInfoDict = {}

## For each MEI:
for MEIObj in VCFObj.lineList:

    end = (MEIObj.infoDict["BKPB"] if "BKPB" in MEIObj.infoDict else "UNK")

    sourceElementId = MEIObj.chrom + ':' + str(MEIObj.pos) + '-' + str(end)
    #print "** source element ** ", sourceElementId

    ## Add source element strand and TSD length to the dictionary
    srcInfoDict[sourceElementId] = {}
    
    strand = (MEIObj.infoDict["STRAND"] if "STRAND" in MEIObj.infoDict else "UNK")
    TSDlen = (MEIObj.infoDict["TSLEN"] if "TSLEN" in MEIObj.infoDict else "UNK")
    refStatus = ("absent" if MEIObj.alt == "<MEI>" else "present")

    srcInfoDict[sourceElementId]["STRAND"] = strand
    srcInfoDict[sourceElementId]["TSLEN"] = TSDlen
    srcInfoDict[sourceElementId]["REF"] = refStatus
    
#### 4. Read source element master table 
#########################################
header("4. Read source element master table")

# Open output file
outFilePath = outDir + '/PCAWG_sourceElements_metadata.tsv'
outFile = open(outFilePath, 'w')

# Write header:
row = '#cytobandId' + "\t" + 'sourceIdNew' + "\t" + 'sourceIdOld' + "\t" + 'refStatus' + "\t" +  'strand' + "\t" + 'TSDlen' + "\t" + 'novelty' + "\t" + 'activityStatus' + "\n"
outFile.write(row)

# Read table line by line
masterTableFile = open(masterTable, 'r')
for line in masterTableFile:

    line = line.rstrip('\r\n')

    if not line.startswith("#"):
        line = line.split('\t')
        sourceIdOld, sourceIdNew, cytobandId = line
        activityStatus = statusDict[sourceIdNew]
        chrom, coords = sourceIdOld.split(":")
        beg, end = coords.split("-")

        ## Check source element novelty
        # a) Known source element (Reported in Tubio et al. (2014) Science)
        if (chrom in knownDict) and ((beg in knownDict[chrom]) or (end in knownDict[chrom])):              
            novelty = 'known'

        # b) Novel souce element
        else:
            novelty = 'novel'            
        
        ## Gather source element strand and TSD length
        ## If inconsistencies between source element identifier raise an error an skip
        # Problem only affects one element
        if sourceIdNew not in srcInfoDict:
            print "[ERROR] source element coordenate not found: ", sourceIdNew, sourceIdOld          
            continue
                    
        strand = srcInfoDict[sourceIdNew]["STRAND"] 
        TSDlen = srcInfoDict[sourceIdNew]["TSLEN"] 
        refStatus = srcInfoDict[sourceIdNew]["REF"] 

        row = cytobandId + "\t" + sourceIdNew + "\t" + sourceIdOld + "\t" + refStatus + "\t" + strand + "\t" + TSDlen + "\t" + novelty + "\t" + activityStatus + "\n"
        outFile.write(row)
        
####
header("Finished")
