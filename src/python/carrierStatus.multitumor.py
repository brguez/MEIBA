#!/usr/bin/env python
#coding: utf-8


#### FUNCTIONS ####
def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"


#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import os.path
import formats
import os.path
from operator import attrgetter
import pandas as pd

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='...')
parser.add_argument('samples', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
samples = args.samples
donorId = args.donorId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "samples: ", samples
print "donorId: ", donorId
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read input VCF
#######################
header("1. Read input VCF")

VCF = VCF.rstrip('\n')        
VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)  


#### 2. Make dictionary contatining for each MEI and sample if it carries or not the MEI
#############################################################################################
header("2. Make dictionary contatining for each MEI and sample if it carries or not the MEI")

allSampleList = filter(None, samples.split(','))   ## List with all the samples for a given donor

MEIdict = {}

## For each MEI
for MEIObj in VCFObj.lineList:

    ## Select only MEI that passes all the filters:    
    if (MEIObj.filter == "PASS"):
            
        MEI = MEIObj.infoDict["CLASS"] + "_" + MEIObj.infoDict["TYPE"] + "_" + MEIObj.chrom + "_" + str(MEIObj.pos)
        
        carrierSampleList = MEIObj.genotype.split(':')[2].split(';') # Make list of samples carrying the variant   

        MEIdict[MEI] = {}

        ## For each sample check if carries the MEI
        for sampleId in allSampleList:
    
            # A) Carries the MEI
            if (sampleId in carrierSampleList):
                MEIdict[MEI][sampleId] = 1
            
            # B) Not carries the MEI             
            else:
                MEIdict[MEI][sampleId] = 0

#### 3. Convert dictionary into dataframe and compute the number of samples carrying each variant
##################################################################################################
header("3. Convert dictionary into dataframe and compute the number of samples carrying each variant")

## Convert to dataframe
MEIdf = pd.DataFrame(MEIdict).transpose()

## Compute and add as first column the number of samples carrying each variant
MEIdf.insert(0, "nbSamples", MEIdf.sum(axis=1))

## Add as first column the donor id
MEIdf.insert(0, "donorId", donorId)

# Save dataframe into tsv
outFilePath = outDir + '/' + donorId + '.carrierStatus.tsv'
MEIdf.to_csv(outFilePath, sep='\t') 


## End ##
print
print "***** Finished! *****"
print
