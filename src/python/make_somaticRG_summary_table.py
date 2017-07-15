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
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: donorId   tumorType   vcfPath')
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
outPath = outDir + '/supplementaryTable3.tsv'

outFile = open(outPath, 'w')
  
## Write file header in the output file
row = '#icgc_donor_id' + "\t" + 'tumorType' + "\t" + 'chrom' + '\t' + 'beg' + '\t' + 'end' + '\t' + 'delLength' + '\t' + 'type' + '\t' + 'srcElementCoord' + '\t' + 'srcElementType' + '\t' + 'srcElementId' + '\n'
    
outFile.write(row)

inputFile = open(inputPath, 'r')

# Per iteration, read a VCF
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")

    donorId = line[0]    
    tumorType = line[1]
    vcfPath = line[2]

    #print "Processing: ", donorId, tumorType, vcfPath

    # Create VCF object
    VCFObj = formats.VCF()

    ## Raise error if the input VCF is not available
    if not os.path.isfile(vcfPath):

        print "[ERROR] Input file is not available"

    ## Input VCF available
    else:

        # Read VCF and add information to VCF object
        VCFObj.read_VCF(vcfPath)
        
        ## For each MEI
        for MEIObj in VCFObj.lineList:  
            
            ## Select only  L1-mediated deletions 
            if ("GR" in MEIObj.infoDict) and (MEIObj.infoDict["GR"] == "DEL"):
                chrom = MEIObj.chrom
                beg = MEIObj.pos
                end = MEIObj.infoDict['BKPB'] if 'BKPB' in MEIObj.infoDict else 'UNK' 
                rtType = MEIObj.infoDict['TYPE']  
                delLength = abs(int(MEIObj.infoDict['TSLEN'])) if 'TSLEN' in MEIObj.infoDict else 'UNK' 

                # a) Solo insertions
                if (rtType == "TD0"):
                    srcElementCoord = "NA"
                    srcElementType = "NA"
                    srcElementId = "NA"

                # b) Transduction
                else:   
                    srcElementCoord = MEIObj.infoDict['SRC']
                    srcElementType = MEIObj.infoDict['SRCTYPE']
                    srcElementId = MEIObj.infoDict['SRCID'] if 'SRCID' in MEIObj.infoDict else 'NA' 

                ## Write MEI into the output file
                row = donorId + "\t" + tumorType + "\t" + chrom + '\t' + str(beg) + '\t' + end + '\t' + str(delLength) + '\t' + rtType + '\t' + srcElementCoord + '\t' + srcElementType + '\t' + srcElementId + '\n'
                outFile.write(row)
 




## End ##
print
print "***** Finished! *****"
print
