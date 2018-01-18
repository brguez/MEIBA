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

            print line
            ## Skip header
            if not line.startswith("#"):
                line = line.split('\t')

                ## Create germline MEI object:
                chrom = line[0]
                beg = line[1]
                end = line[2]
                MEIClass = line[3]
                database = line[4]

                MEIObj = MEI(chrom, beg, end, MEIClass, database)

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
from operator import itemgetter, attrgetter, methodcaller


## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='...')
parser.add_argument('germlineMEIDb', help='...')
parser.add_argument('sampleId', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
germlineMEIDb = args.germlineMEIDb
sampleId = args.sampleId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "germlineMEIDb: ", germlineMEIDb
print "sampleId: ", sampleId
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

############################################
## 1. Create VCF object and read input VCF #
############################################
header("1. Create VCF object and read input VCF")
VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)


###########################################
## 2. Create germline MEI database object #
###########################################
header("Read germline MEI database")
germlineMEIdb = MEIdb()
germlineMEIdb.read_bed(germlineMEIDb)


#########################################################################
## 3. Intersect input somatic candidate MEI with the germline database ##
#########################################################################
header("3. Intersect input somatic candidate MEI with the germline database")

## For each input MEI
for MEIobj in VCFObj.lineList:

    # Skip germline MEI annotation if pseudogene insertion or L1-mediated rearrangement. Only applicable to L1, Alu, SVA and ERVK insertions
    if (MEIobj.infoDict["TYPE"] == "PSD") or ("GR" in MEIobj.infoDict):
        print "[WARNING] Skip germline MEI annotation since " + MEIobj.infoDict["TYPE"] + " insertion"
        continue

    chrom = MEIobj.chrom
    beg = int(MEIobj.pos) - int(MEIobj.infoDict["CIPOS"]) 
    end = int(MEIobj.infoDict["BKPB"]) if "BKPB" in MEIobj.infoDict else int(MEIobj.pos) + int(MEIobj.infoDict["CIPOS"]) 
    iClass = MEIobj.infoDict["CLASS"]

    info("Checking if " + chrom + ":" + str(beg) + "-" +  str(end) + ":" + iClass + " MEI is in the germline database..." )
                  
    ## There are germline MEI of the same class and in the same chromosome
    if (iClass in germlineMEIdb.MEIDict) and (chrom in germlineMEIdb.MEIDict[iClass]):

        ## For each germline MEI in the same chromosome check if it overlaps the current somatic MEI
        for germlineMEIobj in germlineMEIdb.MEIDict[iClass][chrom].itervalues():

            status, nbBases = overlap(beg, end, germlineMEIobj.beg, germlineMEIobj.end, 10) 

            # PCAWG MEI overlaps a MEI in the 1KGP database 
            if (status != "not"):

                #print "OVERLAP: ", beg, end, germlineMEIobj.beg, germlineMEIobj.end, germlineMEIobj.database               
                MEIobj.infoDict["GERMDB"] = germlineMEIobj.database

                # Redefine info attribute with updated information
                MEIobj.info = MEIobj.make_info()
                break


########################
## 4. Make output VCF ##
########################
header("4. Make output VCF")
outFilePath = outDir + '/' + sampleId + ".germlineDbAnnot.vcf"

# 4.1 Write header
VCFObj.write_header(outFilePath)

# 4.2 Write variants
VCFObj.write_variants(outFilePath)


## End ##
print
print "***** Finished! *****"
print
