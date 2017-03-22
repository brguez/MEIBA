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
    print timeInfo, string

####### CLASSES #######
class cohort():
    """
        .....................

        Methods:
        -
    """
    def __init__(self):
        """
        """
        self.VCFdict = {}


    def read_VCFs(self, inputPath):
        """
        """
        inputFile = open(inputPath, 'r')

        info("Read input VCFs ")

        # Per iteration, read a VCF, generate a VCF object and add it to the cohort
        for line in inputFile:
            line = line.rstrip('\n')
            line = line.split("\t")

            donorId = line[0]
            projectCode = line[1].split("-")[0]
            VCFfile = line[2]

            #print "tiooo: ", donorId, projectCode, VCFfile

            # Create VCF object
            VCFObj = formats.VCF()

            info("Reading " + VCFfile + "...")

            # Input VCF available
            if os.path.isfile(VCFfile):

                # Read VCF and add information to VCF object
                VCFObj.read_VCF(VCFfile)

                # Initialize the donor list for a given project if needed 
                if projectCode not in self.VCFdict:

                    self.VCFdict[projectCode] = []   

                # Add donor VCF to cohort
                self.VCFdict[projectCode].append(VCFObj)

            else:
                print "[ERROR] Input file does not exist"



#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
from operator import itemgetter, attrgetter, methodcaller
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per donor with the following consecutive fields: projectCode donorId vcf_path')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ##

### 1. Initialize cohort object
cohortObj = cohort()

### 2. Read VCF files, create VCF objects and organize them
cohortObj.read_VCFs(inputPath)

### 3. Make a dictionary containing per tumor type the total number of retrotransposition events and the retrotransposition rate
# Retrotransposition rate defined as the average number of retrotransposition events for the donors in a given tumor

rtTumorTypeDict = {}
totalNbEvents = 0

## For each tumor type
for projectCode in cohortObj.VCFdict:

    ## Initialize category counts    
    rtTumorTypeDict[projectCode] = {}
    nbEvents = 0
    nbDonors = 0

    ## For each donor 
    for VCFObj in cohortObj.VCFdict[projectCode]:

        nbDonors += 1

        # For each MEI
        for MEIObj in VCFObj.lineList:  
            
            ## Select only those MEI that pass all the filters:
            if (MEIObj.filter == "PASS"):

                totalNbEvents += 1 
                nbEvents += 1 

    # Retrotransposition rate (average number of events per donor):
    rtRate = float(nbEvents) / float(nbDonors)     

    # Save into dictionary:
    rtTumorTypeDict[projectCode]["nbDonors"] =  nbDonors
    rtTumorTypeDict[projectCode]["nbEvents"] =  nbEvents
    rtTumorTypeDict[projectCode]["rtRate"] = rtRate


print  "totalNbEvents: ", totalNbEvents
print "rtTumorTypeDict: ", rtTumorTypeDict

### 4. Make dataframe with the info gathered in 3
#                   nbDonors    nbEvents    rtRate
# ProjectCode1      X1           Y1           Z1     
# ProjectCode2   X2           Y2           Z2
# ...
# ProjectCodeN

rtTumorTypeDataFrame =  pd.DataFrame(rtTumorTypeDict).transpose()

# Save output into tsv
outFilePath = outDir + '/nbDonors_nbEvents_rtRate_perTumorType.tsv'
rtTumorTypeDataFrame.to_csv(outFilePath, sep='\t') 

## End ##
print
print "***** Finished! *****"
print
