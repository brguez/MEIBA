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


def makeSourceRegions(sourceElements):
    """
    Create database of 5' source regions. The source region is the up to 13kb region that is transduced by the source element. Three posibilities depending on source L1 orientation:

    + orientation:

            >>>>>>>>>>>L1>>>>>>>>>>-------------TD_region--------------->

    - orientation:
    
            <-------------TD_region---------------<<<<<<<<<<L1<<<<<<<<<<<

    Unknown orientation (consider both sides):

            <-------------TD_region---------------??????????L1???????????-------------TD_region--------------->

    Make a dictionary with:
        key1 -> chromosome -> (list of source elements coordinates)
    """
    sourceElements = open(sourceElements, 'r')
    tdRegionsDict = {}

    # Read file line by line
    for line in sourceElements:
        line = line.rstrip('\r\n')
    
        ## Discard header
        if not line.startswith("#"):
            # 1p21.1    	1	102568909	102568923	L1	plus    	0	1	1
            # 1p13.1    	1	116980828	116980842	L1	minus	0	1	1
        
            fieldsList = line.split("\t")
            chrom = fieldsList[1]
            beg = int(fieldsList[2])
            end = int(fieldsList[3])
            orientation = fieldsList[5]

            msg = "Processing source element: " + chrom + '_' + str(beg) + '_' + str(end) + '_' + orientation
            log("TD-REGION", msg)

            ## Forward orientation
            if (orientation == "plus"):     
                tdRegionBeg = end
                tdRegionEnd = end + 13000    
                            
            ## Reverse orientation
            elif (orientation == "minus"): 

                tdRegionBeg = beg - 13000
                tdRegionEnd = beg   

            ## Unknown orientation
            else:        
                tdRegionBeg = beg - 13000
                tdRegionEnd = end + 13000 

            tdRegionCoord = str(tdRegionBeg) + "_" + str(tdRegionEnd)  

            msg = "Transduced region coordinates: " + tdRegionCoord
            log("TD-REGION", msg)

            ## First transduced region in a given chromosome
            if chrom not in tdRegionsDict:
                tdRegionsDict[chrom] = [ tdRegionCoord ]
        
            ## There are already transduced regions in the chromosome
            else:
                tdRegionsDict[chrom].append(tdRegionCoord)
    

        msg = "*********************"
        log("TD-REGION", msg)

    return tdRegionsDict

def overlap(begA, endA, begB, endB):

    """
    Check if both ranges overlap. 2 criteria for defining overlap: 

    ## A) Begin of the range A within the range B         
    #       *beg* <---------range_A---------->                         
    # <---------range_B----------> 
                
    #    *beg* <-------range_A----->
    # <-------------range_B------------------>

    ## B) Begin of the range B within the range A     
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->
            
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>
    """    
       
    # a) Begin of the range A within the range B   
    if ((begA >= begB) and (begA <= endB)):
        overlap = True
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True

    # c) Ranges do not overlapping
    else:
        overlap = False

    return overlap

#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import os.path
import formats
import os.path
from operator import attrgetter

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('VCF', help='VCF with processed pseudogene insertion calls for a given donor')
parser.add_argument('sourceElements', help='Bed file containing germline source elements')
parser.add_argument('donorId', help='Donor id. Used for naming output VCF')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
sourceElements = args.sourceElements
donorId = args.donorId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", inputVCF
print "sourceElements: ", sourceElements
print "donorId: ", donorId
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

outFilePath = outDir + '/' + donorId + ".somatic.PSD.filteredTD.vcf"

#### 1. Create database of source elements
tdRegionsDict = makeSourceRegions(sourceElements)

#### 2. Create somatic VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

#### 3. Filter out those PSD insertions that are likely transduction events (the mobilized exons are nearby a source element)

# Iterate over each somatic PSD in the VCF
for PSDObj in VCFObj.lineList:

    msg ="Processing PSD: " + PSDObj.chrom + "_" + str(PSDObj.pos) + " " + PSDObj.infoDict["SRCGENE"] + " " + PSDObj.infoDict["TDC"] 
    subHeader(msg)
    
    chrom, beg, end = PSDObj.infoDict["TDC"].split("_")
    
    ## There are transduced regions for the same chromosome as the mobilized pseudogene
    if chrom in tdRegionsDict:
 
        ## Iterate over the transduced regions list checking if their coordinates match with the mobilized exons called as pseudogenes
        for tdRegion in tdRegionsDict[chrom] :

            ## Define ranges for searching for overlap (extend ranges +- 500 to allow some desviation in the coordinates)
            # Assessed PSD range                
            chrom = chrom
            beg = int(beg) - 500
            end = int(end) + 500          

            # Transduced region range
            tdRegionList = tdRegion.split("_") 
            tdBeg = int(tdRegionList[0]) - 500
            tdEnd = int(tdRegionList[1]) + 500

            ## Assessed PSD is a transduction (mobilized region up to 13 kb down-stream of a source element)
            if (overlap(beg, end, tdBeg, tdEnd)):
                print "     PSD is a transduction!"
 
                ## Insertion passed all the filters   
                if (PSDObj.filter == "PASS"):
                    PSDObj.filter = "TD"
                   
                ## Insertion already filtered. Append TD as a filtering reason
                else:
                     PSDObj.filter = PSDObj.filter + ";TD"
                                         
#### 4. Make output VCF
## 4.1 Write header
VCFObj.write_header(outFilePath)

## 4.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print
print "***** Finished! *****"
print
