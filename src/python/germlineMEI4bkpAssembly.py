#!/usr/bin/env python
#coding: utf-8

### Functions ###
def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"


def average_read_len(bam, sampleSize):
    """
        Compute the average read length based on a subset of reads. Reads selected following bam sorting order
        It is much faster than randomly subsampling and I would expect the output to be the same.  
    """

    bamFile = pysam.AlignmentFile(BAMPath, "rb")

    nbReads = 0
    readLengthList = []

    for alignment in bamFile.fetch():

        # Only select properly mapped read pairs
        if (alignment.is_unmapped == False) and (alignment.is_proper_pair == True):

            readLen = alignment.infer_read_length()
            readLengthList.append(readLen)        
            nbReads += 1

            if (nbReads == sampleSize):
                break

    avReadLen = int(numpy.mean(readLengthList))
    bamFile.close()

    return avReadLen



## Load modules/libraries
import sys
import argparse
import os
import errno
import pysam
import numpy
import time

## Get user's input
parser = argparse.ArgumentParser(description='Produce a correctly formated file for executing breakpoint assembly pipeline on the germline insertions from a given sample')
parser.add_argument('insertions', help='TraFiC germline retrotransposon insertions for a given sample')
parser.add_argument('bam', help='Sample BAM file')
parser.add_argument('sampleId', help='Sample id. Output file will be named accordingly')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
insertionsPath = args.insertions
BAMPath = args.bam
sampleId = args.sampleId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "insertions: ", insertionsPath
print "bam: ", BAMPath
print "sampleId: ", sampleId
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print


#### Infer average read length from bam file
#############################################
header("1. Infer average read length from bam file")

sampleSize = 10000

avReadLen = average_read_len(BAMPath, sampleSize)

print "AVERAGE_READ_LENGTH: ", avReadLen


#### Create output file with insertions
########################################
header("2. Create output file with insertions")

outPath = outDir + "/" + sampleId + ".germline.td0.tsv"
outFile = open(outPath, "w" )

# Print header into the output file
header = "#chromPlus" + "\t" + "begPlus" + "\t" + "endPlus" + "\t" + "nbReadsPlus" + "\t" + "classPlus" + "\t" + "readListPlus" + "\t" + "chromMinus" + "\t" + "begMinus" + "\t" + "endMinus" + "\t" + "nbReadsMinus" + "\t" + "classMinus" + "\t" + "readListMinus" + "\t" + "insertionType" + "\t" + "cytobandId" + "\t" + "sourceType" + "\t" + "chromSource" + "\t" + "begSource" + "\t" + "endSource" + "\t" + "strandSource" + "\t" + "tdBeg" + "\t" + "tdEnd" + "\t" + "tdRnaLen" + "\t" + "tdLen" + "\t" + "psdGene" + "\t" + "chromExonA" + "\t" + "begExonA" + "\t" + "endExonA" + "\t" + "chromExonB" + "\t" + "begExonB" + "\t" + "endExonB" + "\t" + "grType" + "\n"
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
        endPlus = str(int(fieldsList[2]) + avReadLen) 
        nbReadsPlus = fieldsList[3]
        classPlus = fieldsList[4]
        readListPlus = fieldsList[5]
        chromMinus = fieldsList[6]
        begMinus = fieldsList[7] 
        endMinus = str(int(fieldsList[8]) + avReadLen)   
        nbReadsMinus = fieldsList[9] 
        classMinus = fieldsList[10]
        readListMinus = fieldsList[11]
            
        insertionType = "TD0" 
         
        ## Set transduction specific fields as not applicable (NA)
        cytobandId = "NA"
        sourceType = "NA"
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

        ## Set rearrangement specific fields as not applicable (NA)
        grType = "NA" 

        # Print insertion information into the output file  
        row = chromPlus + "\t" + begPlus + "\t" + endPlus + "\t" + nbReadsPlus + "\t" + classPlus + "\t" + readListPlus + "\t" + chromMinus + "\t" + begMinus + "\t" + endMinus + "\t" + nbReadsMinus + "\t" + classMinus + "\t" + readListMinus + "\t" + insertionType + "\t" + cytobandId + "\t" + sourceType + "\t" + chromSource + "\t" + begSource + "\t" + endSource + "\t" + strandSource + "\t" + tdBeg + "\t" + tdEnd + "\t" + tdRnaLen + "\t" + tdLen + "\t" + psdGene + "\t" + chromExonA + "\t" + begExonA + "\t" + endExonA + "\t" + chromExonB + "\t" + begExonB + "\t" + endExonB + "\t" + grType + "\n"

        outFile.write(row) 

print "***** Finished! *****"
print

