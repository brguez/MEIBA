#!/usr/bin/env python
#coding: utf-8



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



def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string

def overlap(begA, endA, begB, endB):
    """
    Check if both ranges overlap. 2 criteria for defining overlap: 
    ## A) Begin of the range A within the range B         
    #       *beg* <---------range_A---------->                         
    # <---------range_B----------> 
                
    #    *beg* <-------range_A----->
    # <-------------range_B------------------>
    ## B) Begin of the range B within the range A     
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->
            
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>
    """    
       
    # a) Begin of the range A within the range B   
    if ((begA >= begB) and (begA <= endB)):
        overlap = True
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True

    # c) Ranges do not overlapping
    else:
        overlap = False

    return overlap


def findDiscordantPE(bamFile, chrom, begA, endA, begB, endB):
    """
    """  
    nbDiscordantPE = 0  
        
    # Extract alignments overlapping the genomic region
    iterator = bamFile.fetch(chrom, begA , endA)

    # Iterate over the alignments
    for alignment in iterator:
 
        properPair = alignment.is_proper_pair
        duplicate = alignment.is_duplicate
        matePos = alignment.next_reference_start

        if (properPair == False) and (duplicate == False) and (overlap(matePos, matePos, begB, endB)):
            nbDiscordantPE += 1

    return nbDiscordantPE


def worker(BED, donorIdBamPathList, discordantPEDict):
    """
    """  
    threadId = threading.currentThread().getName()
    info( threadId + ' launched')
    time.sleep(random.uniform(0, 1))

    ## Iterate over donor id/BAM files, performing discordantPE analysis on one donor at a time
    for line in donorIdBamPathList:

        # Extract donorId and BAM file path
        donorId = line[0]
        BAMPath = line[1]
  
        # Incorporate the donor into the dictionary: 
        discordantPEDict[donorId] = {}

        subHeader(threadId +" performing discordantPE analysis in donor " + donorId)
        time.sleep(random.uniform(0, 1))

        ## Open donor's BAM file for reading
        bamFile = pysam.AlignmentFile(BAMPath, "rb")

        BEDFile = open(BED, 'r')

        # Read file line by line
        for line in BEDFile:

            line = line.rstrip('\r\n')
    
            ## Discard header
            if not line.startswith("#"):
        
                fieldsList = line.split("\t")
                chrom = fieldsList[0]
                sourceId = chrom + ':' + fieldsList[1] + '-' + fieldsList[2]
                begA = int(fieldsList[1]) - 500
                endA = int(fieldsList[1]) + 500

                begB = int(fieldsList[2]) - 500
                endB = int(fieldsList[2]) + 500

                nbDiscordantPE = findDiscordantPE(bamFile, chrom, begA, endA, begB, endB)

                discordantPEDict[donorId][sourceId] = nbDiscordantPE

        subHeader("Finished " + donorId + " discordantPE analysis")
        bamFile.close()

    info( threadId + ' finished' ) 

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import numpy as np
import pandas as pd
from operator import itemgetter, attrgetter, methodcaller
import pysam
import threading
import random

# Global variables:
global debugBool ## debug logging mode. Boolean.

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('BED', help='')
parser.add_argument('BAMPaths', help='Tab separated file with two columns: 1) donor ids and 2) path to the corresponding donor normal BAM files.')
parser.add_argument('-t', '--threads', default=1, dest='threads', type=int, help='Number of threads. Default: 1' )
parser.add_argument('--debug', action='store_true', dest='debug', help='Enable debug mode. Display detailed log info about MEI genotyping')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
BED = args.BED
BAMPaths = args.BAMPaths
threads = args.threads
debugBool = args.debug
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "BED: ", BED
print "BAMPaths: ", BAMPaths
print "threads: ", threads
print "debug: ", debugBool
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

discordantPEDict = {}

# Make list of donorId/BAM files:
donorIdBamPathList = [line.split() for line in open(BAMPaths)]

# Split donors into X evenly sized chunks. X = number of threads
chunkSize = int(round(len(donorIdBamPathList) / float(threads)))

BAMChunks = [donorIdBamPathList[i:i + chunkSize] for i in xrange(0, len(donorIdBamPathList), chunkSize)]

threads = list()
counter = 1

# Generate a discordantPE thread per set of donors
for chunk in BAMChunks:

    print "chunk" + str(counter) + ": " + str(len(chunk)) + " donors for discordantPE analysis"

    threadName = "THREAD-" + str(counter)
    thread = threading.Thread(target=worker, args=(BED, chunk, discordantPEDict), name=threadName)
    threads.append(thread)

    counter += 1

# Launch threads
[t.start() for t in threads]

# Wait till threads have finished discordantPE analyses
[t.join() for t in threads]


### Convert dictionary into dataframe and write into tsv
info("Convert dictionary into dataframe and write into tsv")
discordantPEDf = pd.DataFrame(discordantPEDict)

## Save output into tsv
outFilePath = outDir + '/srcElements_discordantPE.tsv'
discordantPEDf.to_csv(outFilePath, sep='\t') 

print "discordantPEDict: ", discordantPEDict

header("Finished")        
