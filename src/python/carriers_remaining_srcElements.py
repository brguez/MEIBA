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
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
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


##### 2. 
###############################################

header("2. ")


## For each MEI
for MEIObj in VCFObj.lineList:

    print MEIObj.chrom, MEIObj.pos, MEIObj.infoDict["BKPB"]

    ## For each donor
    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():
    
        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        ## Determine number of kbp affected by the current variant
        # a) Homozygous alternative        
        if (genotype == "1/1"): 
            print "HOMOZYGOUS_DONOR: ", donorId, genotypeField
                
        # b) Heterozygous or haploid carrier (for male variants in the X or Y chromosomes outside the PAR region)
        elif (genotype == "0/1") or (genotype == "1"):
            print "HETEROZYGOUS_DONOR: ", donorId, genotypeField
        # c) Donor do not carrying the variant or missing genotype
        #else:
        
    print "*********"

#### END
header("Finished")


