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

def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string


def overlap(begA, endA, begB, endB):

    """
    Check if both ranges overlap. 3 criteria for defining overlap: 

    ## A) Range A within range B
    # <-------------range_B----------------->
    #     <-------range_A------>

    ## B) Begin of the range A within the range B         
    #       *beg* <---------range_A---------->                         
    # <---------range_B----------> 
                
    #    *beg* <-------range_A----->
    # <-------------range_B------------------>

    ## C) Begin of the range B within the range A     
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->
         
    """    
       
    # a) A within B. beg within B and end within B
    if ((begA >= begB) and (begA <= endB)) and ((endA >= begB) and (endA <= endB)):
        overlap = True
        AwithinB = True

    # a) Begin of the range A within the range B   
    elif ((begA >= begB) and (begA <= endB)):        
        overlap = True
        AwithinB = False
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True
        AwithinB = False

    # c) Ranges do not overlapping
    else:
        overlap = False
        AwithinB = False

    return overlap, AwithinB


#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import argparse
import time
import sys
import os.path
import formats
import os.path
import numpy as np
import pybedtools 


## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('chromLen', help='')
parser.add_argument('annot', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
chromLen = args.chromLen
annot = args.annot
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "chromLen: ", chromLen
print "annot: ", annot
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Generate genomic intervals:
# 1Mb long

intervalsPath = outDir +'/intervals.txt'
intervalsFile = open(intervalsPath, 'w')

chromLenFile = open(chromLen, 'r')

windowLen = 1000000

for line in chromLenFile:
    line = line.rstrip('\n')
    lineList = line.split('\t')
    chrom = lineList[0]
    chromEnd = int(lineList[1])

    beg = 0

    for end in np.arange(windowLen, (chromEnd + 1), windowLen):

        row = str(chrom) + "\t" + str(beg) + "\t" + str(end) + "\n"
        intervalsFile.write(row)
        beg = end   


### 2. Intersect with annotation
intersectionPath = outDir +'/intersection.txt'
#command = 'bedtools intersect -wao -a ' + intervalsPath + ' -b ' + annot + ' > ' +  intersectionPath
# os.system(command)
# Executed outside the script. For some reason gives problems within the script... 
# Improper formated row in chr20 leading to truncated output without all the expected windows. 



## Sort intersection
sortedPath = outDir +'/intersection.sorted.txt'
command = 'sort -k1 -k2 -k5 -k6 -n ' + intersectionPath + ' > ' +  sortedPath

os.system(command)


### 3. Read intersected file and for each window compute the number of bases covered by genes 
sortedIntersectionFile = open(sortedPath, 'r')

nbBasesDict = {} 

for line in sortedIntersectionFile:
    line = line.rstrip('\n')
    lineList = line.split()
    
    #print "lineList: ", lineList
    chrom = lineList[0]
    beg = int(lineList[1])
    end = int(lineList[2])
    interval = chrom + "_" + str(beg) + "_" + str(end)
    begGnA = int(lineList[4])
    endGnA = int(lineList[5])
    
    # A) First instance of the interval
    if (interval not in nbBasesDict):
        nbBases = int(lineList[7])
        nbBasesDict[interval] = nbBases    

        # Update previous gene coordinates
        begGnB = begGnA
        endGnB = endGnA

    # B) Already reported intersection for this interval
    else:    
        overlappingGenes, AwithinB = overlap(begGnA, endGnA, begGnB, endGnB)   
 
        #print "OVERLAP: ", begGnA, endGnA, begGnB, endGnB, overlappingGenes, AwithinB

        # a) Not overlapping genes  
        if (overlappingGenes == False):
            nbBases = int(lineList[7])
            nbBasesDict[interval] = nbBases + nbBasesDict[interval]

            # Update previous gene coordinates
            begGnB = begGnA
            endGnB = endGnA
            
        # b) Overlapping genes
        else:

            ## Gene A not within Gene B interval
            if (AwithinB == False):
                overlapLen = endGnB - begGnA
                totalNbBases = int(lineList[7])  
                nbBases = totalNbBases - overlapLen  
                if nbBases < 0  : nbBases = 0
                nbBasesDict[interval] = nbBases + nbBasesDict[interval]

                # Update previous gene coordinates
                begGnB = begGnA
                endGnB = endGnA
            
            ## Gene A within Gene B interval
            # genomic region already covered by Gene B. 
            else:
                begGnB = begGnB
                endGnB = endGnB
    

    #print "NB-BASES: ", nbBasesDict[interval]

   
#print "nbBasesDict: ", nbBasesDict

### 4. Compute the percentage of bases covered by genes and print into file:
geneDensityPath = outDir +'/geneDensity_1Mb.txt'
geneDensityFile = open(geneDensityPath, 'w')

for interval in nbBasesDict:
    
    percBases = float(nbBasesDict[interval]) / windowLen
    chrom, beg, end = interval.split("_")
    row = str(chrom) + "\t" + str(beg) + "\t" + str(end) + "\t" + str(percBases) + "\n"
    geneDensityFile.write(row)

## End ##
print
print "***** Finished! *****"
print
