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


def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


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


def organizeMEI(MEIList):
        """
        Organize MEI into a nested dictionary containing for each insertion class and chromosome 
        a dictionary of germline MEI objects organized by coordinates. Structure: 
        
        key1(MEIclass) -> value(dict2) -> key2(chrId) -> value(MEIObjList)
        MEIObjList sorted in increasing coordinates ordering
        """

        MEIDict = {}
 
        ### Build dictionary
        # Per MEI object 
        for MEIObj in MEIList:

            # a) First MEI of a given class    
            if MEIObj.infoDict["CLASS"] not in MEIDict:
            
                MEIDict[MEIObj.infoDict["CLASS"]] = {}
                MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom] = [MEIObj]
            
            # b) First MEI of a given class in the chromosome
            elif MEIObj.chrom not in MEIDict[MEIObj.infoDict["CLASS"]]:
                MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom] = [MEIObj]

            # c) There are already MEI of this class in the chromosome
            else:
                MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom].append(MEIObj)
                
        # Sort MEI list in increasing coordinates ordering
        for MEIclass in MEIDict:
            for chrom in MEIDict[MEIclass]:

                MEIlist = MEIDict[MEIclass][chrom]
                MEIlistSorted = sorted(MEIlist, key=lambda MEIObj: MEIObj.pos)
                MEIDict[MEIclass][chrom] = MEIlistSorted
            
        return MEIDict


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('VCF', help='multi-sample VCF file containing genotyped MEI')
parser.add_argument('BED', help='BED file containing a row per 1KGP variant and its MAF in European donors')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
BED = args.BED
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "BED: ", BED
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"


## Start ## 

#### 1. Read PCAWG VCF 
########################

header("1. Read PCAWG VCF")

## ReadVCF
VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(VCF)


##### 2. Final number of Alu, L1, SVA and ERVK
###############################################

header("2. Final number of Alu, L1, SVA and ERVK")

nbMEIDict = {}
nbMEIDict["total"] = 0
nbMEIDict["Alu"] = 0
nbMEIDict["L1"] = 0
nbMEIDict["SVA"] = 0
nbMEIDict["ERVK"] = 0

## For each PCAWG MEI 
for MEIObj in VCFObj.lineList:
    
    ## Select only those MEI that passes all the filters    
    if (MEIObj.filter == "PASS"):
        MEIclass = MEIObj.infoDict["CLASS"]
        nbMEIDict["total"] += 1
        nbMEIDict[MEIclass] += 1 

for MEIclass in nbMEIDict:
    nbMEI = nbMEIDict[MEIclass]
    print "nb_MEI: ", MEIclass, nbMEI

##### 3. Median size of SV sites.
################################
# Median size of the inserted mobile elements (Alu, L1, SVA). For ERV we do not have this info

header("3. Median size of SV sites")

sizesDict = {}
sizesDict["Alu"] = []
sizesDict["L1"] = []
sizesDict["SVA"] = []

sizesDelDict = {}
sizesDelDict["Alu"] = []
sizesDelDict["L1"] = []
sizesDelDict["SVA"] = []

## For each PCAWG MEI 
for MEIObj in VCFObj.lineList:

    ## Select only those MEI that passes all the filters and their length info is available:
    if (MEIObj.filter == "PASS") and ("LEN" in MEIObj.infoDict):
        MEIclass = MEIObj.infoDict["CLASS"]
        MEILen = int(MEIObj.infoDict["LEN"])
        sizesDict[MEIclass].append(MEILen)
        
        ## Select 5' deleted variants:
        if (MEIObj.infoDict["STRUCT"] == "DEL"):
            sizesDelDict[MEIclass].append(MEILen)
   
medianSizeDict = {}
    
for MEIclass in sizesDict:
    medianSize = np.median(sizesDict[MEIclass])
    medianSizeDict[MEIclass] = medianSize

    print "Median_size: ", MEIclass, int(round(medianSize))


medianSizeDelDict = {}
    
for MEIclass in sizesDelDict:
    medianSize = np.median(sizesDelDict[MEIclass])
    medianSizeDelDict[MEIclass] = medianSize

    print "Median_size_del: ", MEIclass, int(round(medianSize))



# 4. Median kbp per donor
###############################
## How many kilobase-pairs are affected by the respective variant type per person (median)? 
# I will assume that kilobase-pairs affected would be the same as how many base pairs are different with respect to the reference genome due to a given variant. 

### I see two possible criteria to compute this:
## 1) Consider both the TSD and the MEI length
# For a given MEI. The number of bases affected would be: MEI_length + abs(TSD_length) (it could be duplication or deletion)
# The number of kbp affected in a given donor would be: (MEI_length_1 + TSD_length_1)+(MEI_length_2 + TSD_length_2)+...+(MEI_length_N + TSD_length_N)/1000
# Limitation: we do not have the MEI_length for every MEI (missing for ERVK and for those elements 5' inverted...). What should we do in these cases? infra-estimate not counting the length for these cases? use a consesus length for the element?

## 2) Consider only the TSD length (I will use this criteria since I think it is more )
# For a given MEI. The number of bases affected would be: abs(TSD_length) (it could be duplication or deletion)
# The number of kbp affected in a given donor would be: (TSD_length_1)+(TSD_length_2)+...+(TSD_length_N)/1000
# Limitation: we are infra-estimating the kilobase-pairs affected...

## Note that applies both for 1) and 2). If individual have two copies of the variant multiply by 2 the number of bases affected. 

header("4. Median kbp per donor")

# Initialize empty dictionaries:
nbKbpPerDonorDict = {}

nbKbpPerDonorDict["Alu"] = {}
nbKbpPerDonorDict["L1"] = {}
nbKbpPerDonorDict["SVA"] = {}

nbKbpPerDonorDict2 = {}
nbKbpPerDonorDict2["Alu"] = {}
nbKbpPerDonorDict2["L1"] = {}
nbKbpPerDonorDict2["SVA"] = {}

for donorId in donorIdList:

    nbKbpPerDonorDict["Alu"][donorId] = 0
    nbKbpPerDonorDict["L1"][donorId] = 0
    nbKbpPerDonorDict["SVA"][donorId] = 0

    nbKbpPerDonorDict2["Alu"][donorId] = 0
    nbKbpPerDonorDict2["L1"][donorId] = 0
    nbKbpPerDonorDict2["SVA"][donorId] = 0

## For each MEI
for MEIObj in VCFObj.lineList:

    ## Select only those MEI that passes all the filters and discard ERVK as their length info is not available:
    if (MEIObj.filter == "PASS") and (MEIObj.infoDict["CLASS"] != "ERVK"):
        
        MEIClass = MEIObj.infoDict["CLASS"]
        TSlen = abs(float(MEIObj.infoDict["TSLEN"]))
        MEIlen = float(MEIObj.infoDict["LEN"] if "LEN" in MEIObj.infoDict else float(medianSizeDelDict[MEIClass])) 
    

        ## For each donor
        for donorId, genotypeField in MEIObj.genotypesDict.iteritems():
    
            genotypeFieldList = genotypeField.split(":")
            genotype = genotypeFieldList[0]

            ## Determine number of kbp affected by the current variant
            # a) Homozygous alternative        
            if (genotype == "1/1"): 
                #nbKbpMEI = 2*TSlen/1000
                #nbKbpMEI2 = (2*(TSlen+MEIlen))/1000
                nbKbpMEI = TSlen/1000
                nbKbpMEI2 = (TSlen+MEIlen)/1000
                
            # b) Heterozygous or haploid carrier (for male variants in the X or Y chromosomes outside the PAR region)
            elif (genotype == "0/1") or (genotype == "1"):
                nbKbpMEI = TSlen/1000
                nbKbpMEI2 = (TSlen+MEIlen)/1000

            # c) Donor do not carrying the variant or missing genotype
            else:
                nbKbpMEI = float(0)
                nbKbpMEI2 = float(0)
        
            ## Update the counters incorporating the number of kbp affected by the current variant
            nbKbpPerDonorDict[MEIClass][donorId] += nbKbpMEI
            nbKbpPerDonorDict2[MEIClass][donorId] += nbKbpMEI2

for MEIclass in nbKbpPerDonorDict:

    median = np.median(nbKbpPerDonorDict[MEIclass].values())
        
    print "nbKbp_median: ", MEIclass, median

for MEIclass in nbKbpPerDonorDict2:

    median = np.median(nbKbpPerDonorDict2[MEIclass].values())
        
    print "nbKbp_median2: ", MEIclass, median


##### 5. Sensitivity estimate 
################################
# Approach: Proportion (in %) of 1000 Genomes Project Phase 3 variant sites with >1% MAF (common) from European Populations that are recovered in PCAWG

header("5. Sensitivity estimate")

### 5.1 Organize MEI into a dictionary
MEIDict = organizeMEI(VCFObj.lineList)

### 5.2 Read 1KGP BED file and check for each variant if is reported in PCAWG
BEDFile = open(BED, 'r')

snDict = {}
snDict["Alu"] = {}
snDict["L1"] = {}
snDict["SVA"] = {}
snDict["Alu"]["totalNbCommon"] = 0
snDict["L1"]["totalNbCommon"] = 0
snDict["SVA"]["totalNbCommon"] = 0
snDict["Alu"]["nbCommonInPCAWG"] = 0
snDict["L1"]["nbCommonInPCAWG"] = 0
snDict["SVA"]["nbCommonInPCAWG"] = 0

# For each 1KGP MEI
for line in BEDFile: 

    line = line.rstrip('\r\n')

    # Skip header
    if not line.startswith("#"):
        line = line.split('\t')

        chrom = line[0]
        beg1KGP = int(line[1]) - 10
        end1KGP = int(line[2]) + 10
        MEIclass = line[3]
        EURMAF = float(line[4])
        
        ## Only consider MEI in autosomal chromosomes with >1% MAF (common variants)
        if (chrom != "X") and (EURMAF > 1):        

            snDict[MEIclass]["totalNbCommon"] += 1 

            ## For each PCAWG MEI 
            for MEIObj in MEIDict[MEIclass][chrom]:

                ## Select only those MEI that passes all the filters    
                if (MEIObj.filter == "PASS"):

                    begPCAWG, endPCAWG = insertionRange(MEIObj, 10)

                    ## 1KGP MEI reported in PCAWG
                    if (overlap(beg1KGP, end1KGP, begPCAWG, endPCAWG)):
                        snDict[MEIclass]["nbCommonInPCAWG"] += 1 
                        break

for MEIclass in snDict:

    totalNbCommon = float(snDict[MEIclass]["totalNbCommon"])
    commonInPCAWG = float(snDict[MEIclass]["nbCommonInPCAWG"])
    sn = float(commonInPCAWG)/float(totalNbCommon)*100
    print "sn_estimation: ", MEIclass, round(sn, 2)


#### END
header("Finished")


