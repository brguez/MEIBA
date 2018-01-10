#!/usr/bin/env python
#coding: utf-8

#### FUNCTIONS ####
def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, "****", string, "****"

def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string

def overlap(begA, endA, begB, endB, overhang):

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
    
    begA = max(0, (begA - overhang)) # Use 0 if it becomes negative 
    endA = endA + overhang
    begB = max(0, (begB - overhang)) # Use 0 if it becomes negative 
    endB = endB + overhang
 

    # A) Range A within range B
    #    *beg* <-------range_A----->   
    # <-------------range_B------------------>
    if ((begA >= begB) and (endA <= endB)):

        nbBases = endA - begA + 1
        overlap = "within"

    # B) Range B within range A
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>
    elif ((begB >= begA) and (endB <= endA)):

        nbBases = endB - begB + 1
        overlap = "within"

    # C) Range A partially overlapping range B
    #           *beg* <---------range_A----------> 
    #   <---------range_B---------->    
    elif ((begA > begB) and (begA <= endB)):
      
        nbBases = endB - begA + 1
        overlap = "partial"

    # D) Range B partially overlapping range A 
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->            
    elif ((begB > begA) and (begB <= endA)):

        nbBases = endA - begB + 1
        overlap = "partial"

    # E) Ranges do not overlapping
    else:
        nbBases = 0
        overlap = "not"

    return overlap, nbBases

#### CLASSES ####
class MEIdb():
    """
    """
    def __init__(self):
        """
        """
        self.MEIDict = {}

    def read_bed(self, inputBed):
        """
        """
        # Open germline MEI database input file
        inputBedFile = open(inputBed, 'r')

        ## Read first input file line by line
        for line in inputBedFile:
            line = line.rstrip('\n')

            ## Skip header
            if not line.startswith("#"):
                line = line.split('\t')

                ## Create germline MEI object:
                chrom = line[0]
                beg = line[1]
                end = line[2]
                MEIClass = line[3]

                MEIObj = MEI(chrom, beg, end, MEIClass, "1KGP")

                ## Create nested dictionaries as needed
                # A) First MEI of a given class
                if MEIObj.MEIClass not in self.MEIDict:
                    self.MEIDict[MEIObj.MEIClass] = {}
                    self.MEIDict[MEIObj.MEIClass][MEIObj.chrom] = {}

                # B) First MEI of a given class in the chromosome
                elif MEIObj.chrom not in self.MEIDict[MEIObj.MEIClass]:
                    self.MEIDict[MEIObj.MEIClass][MEIObj.chrom] = {}

                ## Add germline MEI object to the dictionary
                identifier = str(MEIObj.beg) + '-' + str(MEIObj.end)

                self.MEIDict[MEIObj.MEIClass][MEIObj.chrom][identifier] = MEIObj

        # Close germline MEI database file
        inputBedFile.close()


class MEI():
    """
    """
    def __init__(self, chrom, beg, end, MEIClass, database):
        """
        """
        self.chrom = chrom
        self.beg = int(beg)
        self.end = int(end)
        self.MEIClass = MEIClass
        self.database = database



#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import numpy as np
import scipy.stats as stats
from operator import itemgetter, attrgetter, methodcaller


## Get user's input ##
parser = argparse.ArgumentParser(description= "Compare PCAWG with 1KGP MEI and obtain different statistics")
parser.add_argument('VCF', help='PCAWG MEI')
parser.add_argument('oneKGPdb', help='BED file containing 1KGP MEI')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
oneKGPdb = args.oneKGPdb
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "oneKGPdb: ", oneKGPdb
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

###########
## Start ## 
###########

############################################
## 1. Create VCF object and read input VCF #
############################################
header("1. Create VCF object and read input VCF")
VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)

###########################################
## 2. Create germline MEI database object #
###########################################
header("2. Create germline MEI database object")
MEIdb1KGP = MEIdb()
MEIdb1KGP.read_bed(oneKGPdb)

######################################
## 3. Intersect PCAWG MEI with 1KGP ##
######################################
header("3. Intersect PCAWG MEI with 1KGP")

mergedMEIdbObj = MEIdb()

## For each PCAWG MEI
for MEIobjPCAWG in VCFObj.lineList:

    ## Select only those MEIs that pass the filters
    if (MEIobjPCAWG.filter == "PASS"):

        database = "PCAWG"

        chromPCAWG = MEIobjPCAWG.chrom
        begPCAWG = int(MEIobjPCAWG.pos)
        endPCAWG = int(MEIobjPCAWG.infoDict["BKPB"])
        classPCAWG = MEIobjPCAWG.infoDict["CLASS"]

        info("Checking if " + chromPCAWG + ":" + str(begPCAWG) + "-" +  str(endPCAWG) + " MEI is in 1KGP dataset..." )
                  
        ## There are 1KGP MEI of the same class and in the same chromosome
        if (classPCAWG in MEIdb1KGP.MEIDict) and (chromPCAWG in MEIdb1KGP.MEIDict[classPCAWG]):

            ## For each 1KGP MEI in the same chromosome check if it overlaps the current PCAWG MEI
            for MEIid1KGP, MEIobj1KGP in MEIdb1KGP.MEIDict[classPCAWG][chromPCAWG].iteritems():

                status, nbBases = overlap(begPCAWG, endPCAWG, MEIobj1KGP.beg, MEIobj1KGP.end, 10) 

                # PCAWG MEI overlaps a MEI in the 1KGP database 
                if (status != "not"):
                    database = "PCAWG,1KGP"
                    MEIobj1KGP.database = database
                    MEIdb1KGP.MEIDict[classPCAWG][chromPCAWG][MEIid1KGP] = MEIobj1KGP

                    #print "OVERLAP: ", begPCAWG, endPCAWG, MEIobj1KGP.beg, MEIobj1KGP.end
                    break

        ## Create MEI object
        MEIObj = MEI(chromPCAWG, begPCAWG, endPCAWG, classPCAWG, database)

        ## Create nested dictionaries as needed
        # A) First MEI of a given class
        if MEIObj.MEIClass not in mergedMEIdbObj.MEIDict:
            mergedMEIdbObj.MEIDict[MEIObj.MEIClass] = {}
            mergedMEIdbObj.MEIDict[MEIObj.MEIClass][MEIObj.chrom] = {}

        # B) First MEI of a given class in the chromosome
        elif MEIObj.chrom not in mergedMEIdbObj.MEIDict[MEIObj.MEIClass]:
            mergedMEIdbObj.MEIDict[MEIObj.MEIClass][MEIObj.chrom] = {}

        ## Add germline MEI object to the merged dictionary
        identifier = str(MEIObj.beg) + '-' + str(MEIObj.end)

        mergedMEIdbObj.MEIDict[MEIObj.MEIClass][MEIObj.chrom][identifier] = MEIObj        


#######################################################
## 4. Add MEI 1KGP specific to the merged dictionary ##  
#######################################################
header("4. Add MEI 1KGP specific to the merged dictionary")

number = 0
## For each class
for MEIClass in MEIdb1KGP.MEIDict:

    ## For each chromosome
    for MEIchrom in MEIdb1KGP.MEIDict[MEIClass]:

        ## For each 1KGP MEI of a given class in a given chromosome
        for MEIid1KGP, MEIobj1KGP in MEIdb1KGP.MEIDict[MEIClass][MEIchrom].iteritems():

            number += 1
        
            ## 1KGP specific MEI:
            if (MEIobj1KGP.database == "1KGP"):

                ## Create nested dictionaries as needed
                # A) First MEI of a given class
                if MEIClass not in mergedMEIdbObj.MEIDict:
                    mergedMEIdbObj.MEIDict[MEIClass] = {}
                    mergedMEIdbObj.MEIDict[MEIClass][MEIchrom] = {}

                # B) First MEI of a given class in the chromosome
                elif MEIchrom not in mergedMEIdbObj.MEIDict[MEIClass]:
                    mergedMEIdbObj.MEIDict[MEIClass][MEIchrom] = {}

                ## Add germline MEI object to the merged dictionary
                mergedMEIdbObj.MEIDict[MEIClass][MEIchrom][MEIid1KGP] = MEIobj1KGP        


################################################################
## 5. Write merged germline MEI (1KGP + PCAWG) into a bed file # 
################################################################
header("5. Write merged germline MEI (1KGP + PCAWG) into a bed file")

outFilePath = outDir + '/germlineMEI_db.bed'
outFile = open(outFilePath, 'w')

## For each class
for MEIClass in mergedMEIdbObj.MEIDict:

    ## For each chromosome
    for MEIchrom in mergedMEIdbObj.MEIDict[MEIClass]:

        ## For each MEI of a given class in a given chromosome
        for MEIobj in mergedMEIdbObj.MEIDict[MEIClass][MEIchrom].itervalues():
            row = MEIobj.chrom + '\t' + str(MEIobj.beg) + '\t' + str(MEIobj.end) + '\t' + MEIobj.MEIClass + '\t' + MEIobj.database + '\n'
            outFile.write(row)

outFile.close()

## End ##
print
print "***** Finished! *****"
print
