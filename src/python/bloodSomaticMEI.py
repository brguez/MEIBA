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
    print timeInfo, string,

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import scipy.stats as stats
import numpy as np
from matplotlib import pyplot as plt

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCFnormal', help='multi-sample VCF file containing genotyped germline MEI in normal genomes')
parser.add_argument('VCFtumor', help='multi-sample VCF file containing genotyped germline MEI in tumor genomes')
parser.add_argument('--max-zygosity', default=1, type=float, dest='maxZygosity', help='Maximum zygosity for considering a MEI as somatic. Default: 1, so disabled.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCFnormal = args.VCFnormal
VCFtumor =  args.VCFtumor
maxZygosity = args.maxZygosity
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCFnormal: ", VCFnormal
print "VCFtumor: ", VCFtumor
print "maxZygosity: ", maxZygosity
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read input VCFs and generate VCF objects
###################################################
# Important requirement: the two VCF must have the same MEIs sorted also in the same order
header("1. Process input VCFs ")


## Normal genome multi-sample VCF
VCFObjNormal = formats.VCF()
donorIdListNormal = VCFObjNormal.read_VCF_multiSample(VCFnormal)

## Tumor genome multi-sample VCF
VCFObjTumor = formats.VCF()
donorIdListTumor = VCFObjTumor.read_VCF_multiSample(VCFtumor)

#### 2. Identify MEI that are blood/normal specific singletons
###############################################################
# These cases are potential blood somatic events...
header("2. Identify MEI that are blood/normal specific singletons")

# Make VCF object with all the candidate blood somatic MEI
VCFObjBloodSomatic = formats.VCF()

### For each MEI
nbMEI = len(VCFObjNormal.lineList)
for index in range(0, nbMEI):

    ### Generate MEI objects containing
    # a) MEI genotypes in normal genomes  
    MEIObjNormal = VCFObjNormal.lineList[index]

    info("Checking if " + MEIObjNormal.chrom + ":" + str(MEIObjNormal.pos) + " MEI is a blood/normal specific singleton..." )

    # b) MEI genotypes in tumor genomes
    MEIObjTumor = VCFObjTumor.lineList[index]

    ### Sanity check to make sure we are comparing the tumor and normal genotype for the same MEI
    if (MEIObjNormal.chrom != MEIObjTumor.chrom) or (MEIObjNormal.pos != MEIObjTumor.pos):
        print '\n[ERROR] VCF do not sorted in the same order'
        sys.exit(2)

    ### Count the number of donors carrying the MEI
    # a) Normal genome genotypes
    nbCarriersNormal = 0

    # For each donor
    carriersIdNormalList = []

    for donorId in donorIdListNormal:
        genotypePlusReadsNormal = MEIObjNormal.genotypesDict[donorId]
        genotypeNormal = genotypePlusReadsNormal.split(':')[0]

        # MEI carrier
        if (genotypeNormal == '1/1') or (genotypeNormal == '0/1'):
            nbCarriersNormal += 1
            carriersIdNormalList.append(donorId)

    # b) Tumor genome genotypes
    nbCarriersTumor = 0

    # For each donor
    carriersIdTumorList = []

    for donorId in donorIdListTumor:
        genotypePlusReadsTumor = MEIObjTumor.genotypesDict[donorId]
        genotypeTumor = genotypePlusReadsTumor.split(':')[0]

        # MEI carrier
        if (genotypeTumor == '1/1') or (genotypeTumor == '0/1'):

            nbCarriersTumor += 1
            carriersIdTumorList.append(donorId)

    ## Blood-specific singleton
    if (nbCarriersNormal == 1) and (nbCarriersTumor == 0):

        carrierId = carriersIdNormalList[0]
        genotypeField = MEIObjNormal.genotypesDict[carrierId]
        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]
        nbReadsMEI = float(genotypeFieldList[1])
        totalNbReads = float(genotypeFieldList[2])
        zygosity =  nbReadsMEI / totalNbReads

        ## And zygosity/VAF <= maxZygosity -> candidate Blood somatic MEI
        if (zygosity <= maxZygosity):
            VCFObjBloodSomatic.addLine(MEIObjNormal)

    # b) Not blood-specific singleton
    else:
        print 'not'


#### 3. Make summary table with candidate Blood somatic MEI
#############################################################
header("3. Make summary table with candidate Blood somatic MEI")

# Open output file
fileName = 'candidate_blood_somatic_MEI.maxZygosity' + str(maxZygosity) + '.tsv'
outFilePath = outDir + '/' + fileName
outFile = open(outFilePath, 'w')

# Write header
header = "#CHROM" + "\t" + "POS" + "\t" + "CLASS" + "\t" + "ZYGOSITY" + "\t" + "DONOR_ID" + "\n"
outFile.write(header)

# For each MEI blood-specific singleton
for MEIObj in VCFObjBloodSomatic.lineList:

    # For each genotyped donor
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        # MEI carrier (there should be just one carrier per MEI since they are singletons)
        if (genotype == '1/1') or (genotype == '0/1'):

            nbReadsMEI = float(genotypeFieldList[1])
            totalNbReads = float(genotypeFieldList[2])


            # Compute genotyping VAF:
            if totalNbReads == 0:
                zygosity = 0
            else:
                zygosity =  nbReadsMEI / totalNbReads

            ## Write row into output table:
            row = MEIObj.chrom + "\t" + str(MEIObj.pos) + "\t" + MEIObj.infoDict['CLASS'] + "\t" + str(zygosity) + "\t" + donorId + "\n"
            outFile.write(row)
