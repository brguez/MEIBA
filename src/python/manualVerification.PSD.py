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


def organizePSD(manualPSD):
    """
    Create database of manually curated PSD insertions 
    # Make a nested dictionary with:
        key1 -> donorUniqueId (projectCode::donorId) 
        key2 -> geneName -> value (list of PSD insertion objects)
    """
    manualPSD = open(manualPSD, 'r')
    manualPSDDict = {}

    # Read file line by line
    for line in manualPSD:
        line = line.rstrip('\r\n')

        ## Discard header
        if not line.startswith("#"):
        
            fieldsList = line.split("\t")
            chrom = fieldsList[0]
            beg = fieldsList[1]
            end = fieldsList[2]
            gene = fieldsList[3]
            projectCode = fieldsList[4]
            donorId = fieldsList[5]
            donorUniqueId = projectCode + "::" + donorId
            coords = chrom + "_" + beg + "_" + end
    
            ## First PSD in the given donor
            if donorUniqueId not in manualPSDDict:

                manualPSDDict[donorUniqueId] = {}
                manualPSDDict[donorUniqueId][gene] = {}
                manualPSDDict[donorUniqueId][gene] = [ coords ]

            ## Already PSD in the donor       
            else:
                
                ## First source gene in the donor
                if gene not in manualPSDDict[donorUniqueId]:
    
                    manualPSDDict[donorUniqueId][gene] = {}
                    manualPSDDict[donorUniqueId][gene] = [ coords ]
    
                ## Already PSD from the current source gene in the donor
                else:
                    manualPSDDict[donorUniqueId][gene].append(coords)
    
    return manualPSDDict

def insertionRange(MEIObj, overhang):
    """
    Define MEI range as follows:
                
                      range
       <------------------------------------>
       beg [bkpA - CIPOS - overhang]        end [(bkpB or bkpA) - CIPOS - overhang]  
    """
    bkpA = MEIObj.pos
    bkpB = int(MEIObj.infoDict["BKPB"]) if "BKPB" in MEIObj.infoDict else MEIObj.pos 
    CIPOS = int(MEIObj.infoDict["CIPOS"])
    begRange = bkpA - CIPOS - overhang
    endRange = bkpB + CIPOS + overhang

    return (begRange, endRange)

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
parser = argparse.ArgumentParser(description= "Set as PASS those PSD insertions that do not pass the filtering score but have been verified with IGV. Add the flag MANUAL to make clear that the call is the result of manual confirmation")
parser.add_argument('VCF', help='VCF with processed pseudogene insertion calls for a given donor')
parser.add_argument('manualPSD', help='Bed file containing manually verified pseudogene insertions through BAM inspection')
parser.add_argument('donorId', help='PCAWG submitter donor Id')
parser.add_argument('projectCode', help='PCAWG project code')

parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
manualPSD = args.manualPSD
donorId = args.donorId
projectCode = args.projectCode

outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", inputVCF
print "manualPSD: ", manualPSD
print "donorId: ", donorId
print "projectCode: ", projectCode
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

outFilePath = outDir + '/' + donorId + ".somatic.PSD.manual.vcf"

#### 1. Create database of manually curated PSD insertions 
# Make a nested dictionary with:
# key1 -> donorUniqueId (projectCode::donorId) 
# key2 -> geneName -> value (list of PSD insertion objects)

manualPSDDict = organizePSD(manualPSD)

#### 3. Create somatic VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

#### 4. Modify filtering status for those PSD manually verified
# Iterate over each somatic PSD in the VCF
donorUniqueId = projectCode + "::" + donorId

for PSDObj in VCFObj.lineList:

    chrom = PSDObj.chrom
    beg = PSDObj.pos
    end = int(PSDObj.infoDict["BKPB"]) if "BKPB" in PSDObj.infoDict else PSDObj.pos
    gene =  PSDObj.infoDict["SRCGENE"]

    msg ="Processing PSD: " + str(chrom) + "_" + str(beg) + "_" + str(end) + " " + gene + " " + PSDObj.infoDict["SCORE"] 
    subHeader(msg)

    ## There are manually verified PSD for this donor
    if donorUniqueId in manualPSDDict:

        ## There are manually verified PSD involving the same source gene
        if gene in manualPSDDict[donorUniqueId]:

            ## Iterate over the verified PSD list checking if their insertion coordinates match with the assessed PSD insertion
            for verifiedPSD in manualPSDDict[donorUniqueId][gene]:
                
                ## Define insertion ranges for searching for overlap (extend ranges +- 500 to allow some desviation in the coordinates)
                # Assessed PSD range                
                chrom = chrom
                beg = beg - 500
                end = end + 500          

                # Verified PSD range
                verifiedPSDlist = verifiedPSD.split("_") 
                print "     PSD-verified-coordinates: ", verifiedPSDlist
               
                vChrom = verifiedPSDlist[0]
                vBeg = int(verifiedPSDlist[1]) - 500
                vEnd = int(verifiedPSDlist[2]) + 500
     
                ## Assessed PSD has been manually verified (overlaps)
                if (chrom == vChrom) and (overlap(beg, end, vBeg, vEnd)):
                    print "     PSD have been manually verified!"
                    PSDObj.filter = "PASS"
                    PSDObj.infoDict["MANUAL"] = "1"

                    PSDObj.info = PSDObj.make_info()
                    
                    break


#### 4. Make output VCF
## 4.1 Write header
VCFObj.write_header(outFilePath)

## 4.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print
print "***** Finished! *****"
print
