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
    print timeInfo, string,

#### CLASSES ####
class germlineMEIDb():
    """
            .....................

            Methods:
            -

    """
    def __init__(self):
        """

        """
        self.germlineMEIDict = {}

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
                pos = line[1]
                MEIClass = line[3]
                family = line[4]
                TSD = line[5]
                orientation = line[6]
                length = line [7]
                germlineMEIObj = germlineMEI(chrom, pos, MEIClass, family, TSD, orientation, length)

                ### Create a nested dictionary containing for each insertion class and chromosome a dictionary of
                ## germline MEI objects organized by coordinates
                # key1(MEIclass) -> value(dict2) -> key2(chrId) -> value(dict3) -> key3(coordinate) -> value(germlineMEIObj)

                ## Create nested dictionaries as needed
                # A) First MEI of a given class
                if germlineMEIObj.MEIClass not in self.germlineMEIDict:
                    self.germlineMEIDict[germlineMEIObj.MEIClass] = {}
                    self.germlineMEIDict[germlineMEIObj.MEIClass][germlineMEIObj.chrom] = {}

                # B) First MEI of a given class in the chromosome
                elif germlineMEIObj.chrom not in self.germlineMEIDict[germlineMEIObj.MEIClass]:
                    self.germlineMEIDict[germlineMEIObj.MEIClass][germlineMEIObj.chrom] = {}

                ## Add germline MEI object to the dictionary
                self.germlineMEIDict[germlineMEIObj.MEIClass][germlineMEIObj.chrom][germlineMEIObj.pos] = germlineMEIObj

        # Close germline MEI database file
        inputBedFile.close()

class germlineMEI():
    """
            .....................

            Methods:
            -

    """
    def __init__(self, chrom, pos, MEIClass, family, TSD, orientation, length):
        """

        """
        self.chrom = chrom
        self.pos = pos
        self.MEIClass = MEIClass
        self.family = family
        self.TSD = TSD
        self.orientation = orientation
        self.length = length


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
parser.add_argument('inputVCF', help='...')
parser.add_argument('germlineMEIBed', help='...')
parser.add_argument('sampleId', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
germlineMEIBed = args.germlineMEIBed
sampleId = args.sampleId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "germlineMEIBed: ", germlineMEIBed
print "sampleId: ", sampleId
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## 2. Create germline MEI database object
germlineMEIDbObj = germlineMEIDb()
germlineMEIDbObj.read_bed(germlineMEIBed)

## 3. For each MEI in the input VCF check if already included in the database

for VCFlineObj in VCFObj.lineList:

    # Skip repeat novelty annotation if pseudogene insertion or L1-mediated deletion. Repeat annotation only applicable to L1, Alu, SVA and ERVK? insertions
    if (VCFlineObj.infoDict["TYPE"] == "PSD") or (VCFlineObj.infoDict["TYPE"] == "DEL"):
        print "[WARNING] Skip repeat novelty annotation since " + VCFlineObj.infoDict["TYPE"] + " insertion"
        continue

    info("Checking if " + VCFlineObj.chrom + ":" + str(VCFlineObj.pos) + " MEI is in 1000 genomes germline database..." )

    # A) There are MEI in the germline database in the same chromosome as the input MEI
    if VCFlineObj.chrom in germlineMEIDbObj.germlineMEIDict[VCFlineObj.infoDict["CLASS"]].keys():
        
	## Make numpy array with the chromosomal positions for the germline MEI on the same chromosome in the 1000 genomes database
        MEIDbPosArr = np.array(sorted(germlineMEIDbObj.germlineMEIDict[VCFlineObj.infoDict["CLASS"]][VCFlineObj.chrom].keys(), key=int), dtype='int')

        ## Identify those MEI in the germline database overlapping the first insertion breakpoint
        # Determine bkpA beg and end range to search for overlap with 1000 genomes MEI database
        begBkpA = VCFlineObj.pos - int(VCFlineObj.infoDict["CIPOS"]) - 10
        endBkpA = VCFlineObj.pos + int(VCFlineObj.infoDict["CIPOS"]) + 10

        # Select those MEI in the database within bkpA beg and end range
        filteredArrBkpA = MEIDbPosArr[(MEIDbPosArr >= begBkpA) & (MEIDbPosArr <= endBkpA)]

        ## A) Target site info available -> Identify those MEI in the germline database overlapping the second insertion breakpoint
        if 'TSLEN' in VCFlineObj.infoDict:

            # a) Consistent TSD
            if (VCFlineObj.infoDict["TSLEN"] != "inconsistent"):
                # Determine second breakpoint coordinates and bkpB beg and end range to search for overlap with 1000 genomes MEI database
                bkpB = VCFlineObj.pos + abs(int(VCFlineObj.infoDict["TSLEN"]))
                begBkpB = bkpB - int(VCFlineObj.infoDict["CIPOS"]) - 10
                endBkpB = bkpB + int(VCFlineObj.infoDict["CIPOS"]) + 10

                # Select those MEI in the database within bkpA beg and end range
                filteredArrBkpB = MEIDbPosArr[(MEIDbPosArr >= begBkpB) & (MEIDbPosArr <= endBkpB)]

            # b) Inconsistent TSD
            else:
                # Empty array
                filteredArrBkpB = np.array([ ])

        ## B) Target site info not available
        else:
            # Empty array
            filteredArrBkpB = np.array([ ])

    # B) There are not MEI in the germline database in the same chromosome as the input MEI
    else:
        # Empty arrays
	    filteredArrBkpA = np.array([])
    filteredArrBkpB = np.array([])

    ## Make a single array containing a not redundant list of 1000 genomes MEI coordinates overlapping the first and/or second insertion breakpoints.
    filteredArr = np.array(np.unique(np.concatenate((filteredArrBkpA, filteredArrBkpB), axis=0)), dtype='int')

    ## A) MEI already reported in the germline database
    if (len(filteredArr) != 0) :
        print 'yes'
        VCFlineObj.infoDict["GERMDB"] = "1KGENOMES"

        # Redefine info attribute with updated information
        VCFlineObj.info = VCFlineObj.make_info()


    ## B) MEI not reported in the germline database
    else:
        print 'not'



## 4. Make output VCF
outFilePath = outDir + '/' + sampleId + ".germlineDbAnnot.vcf"

# 4.1 Write header
VCFObj.write_header(outFilePath)

# 4.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print
print "***** Finished! *****"
print
