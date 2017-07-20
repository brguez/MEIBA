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


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
from operator import itemgetter, attrgetter, methodcaller
import pandas as pd

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per donor with the following consecutive fields: donorId vcf_path')
parser.add_argument('donorMetadata', help='PCAWG donors metadata')
parser.add_argument('srcElements', help='L1 source elements metadata')
parser.add_argument('outFileName', help='Output file name identifier')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
donorMetadata = args.donorMetadata
srcElements = args.srcElements
outFileName = args.outFileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "donorMetadata: ", donorMetadata
print "srcElements: ", srcElements
print "outFileName: ", outFileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ##

#### 1. Read metadata file and create dictionary with germline source element activity
info("1. Read metadata file and create dictionary with germline source element activity")

donorMetadataFile = open(donorMetadata, 'r')
germlineSrcActivityDict = {}

## For donor's metadata:
for line in donorMetadataFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
        donorId = line[1]
        exclusionStatus = line[3]

        if (exclusionStatus == "Whitelist"):

            germlineSrcActivityDict[donorId] = {}

#### 2. Read germline source elements file and initialize dictionary with source element activity
info("2. Read germline source elements file and initialize dictionary with source element activity")

srcElementsFile = open(srcElements, 'r')

# For each source element
for line in srcElementsFile:
    
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
        srcId = line[0]
         
        ## For each donor
        for donorId in germlineSrcActivityDict.keys():

            germlineSrcActivityDict[donorId][srcId] = 0
        

#### 2. Read donor's VCF files containing the somatic MEI and fill the dictionary
info("2. Read input VCFs")

inputFile = open(inputPath, 'r')
somaticSrcActivityDict = {}

# Per iteration process a donor's VCF
for line in inputFile:
    
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
    
        donorId = line[0]
        VCFfile = line[2]

        # Create VCF object
        VCFObj = formats.VCF()

        info("Reading " + VCFfile + "...")

        # Input VCF available
        if os.path.isfile(VCFfile):

            ## If not excluded donor:
            if donorId in germlineSrcActivityDict:
    
                # Read donor's VCF and add information to VCF object
                VCFObj.read_VCF(VCFfile)
    
                ## For each somatic MEI in the VCf
                for MEIObj in VCFObj.lineList:
    
                    # Select transductions that passess all the filters. Exclude L1-mediated rearrangements mediated by transductions:
                    if ("GR" not in MEIObj.infoDict) and (MEIObj.filter == "PASS") and ((MEIObj.infoDict["TYPE"] == "TD1") or (MEIObj.infoDict["TYPE"] == "TD2")):

                        # A) germline source element                     
                        if (MEIObj.infoDict["SRCTYPE"] == "GERMLINE"):
 
                            srcId = MEIObj.infoDict["SRCID"] # Use the cytoband as the source element identifier                       
                            germlineSrcActivityDict[donorId][srcId] += 1  
                    
                        # B) somatic germline element
                        else:
                            srcId = MEIObj.infoDict["SRC"]
                            print "SOMATIC: ", srcId, MEIObj.infoDict["SRCTYPE"]                               

                            # a) somatic source element not reported yet
                            if srcId not in somaticSrcActivityDict:
                                somaticSrcActivityDict[srcId] = {}
                                somaticSrcActivityDict[srcId][donorId] = 1

                            # b) somatic source element not reported yet in a given donor
                            elif donorId not in somaticSrcActivityDict[srcId]:
                                somaticSrcActivityDict[srcId][donorId] = 1
                    
                            # c) somatic source element already reported in the donor
                            else:
                                somaticSrcActivityDict[srcId][donorId] += 1
                           
### 3. Convert germline dictionary into dataframe and write into tsv
info("3. Convert dictionary into dataframe and write into tsv")
germlineSrcActivityDf = pd.DataFrame(germlineSrcActivityDict)

print "germlineSrcActivityDf: ", germlineSrcActivityDf

sys.exit(1)

## Save output into tsv
outFilePath = outDir + '/germline_' +  outFileName +'.tsv'
germlineSrcActivityDf.to_csv(outFilePath, sep='\t') 

### 4. Generate tsv with the somatic source elements activity:
## tsv with 5 Fields:
# srcId, nbSamples, maxActivity, sampleList, activityList

info("4. Generate tsv with the somatic source elements activity")

outFilePath = outDir + '/somatic_' +  outFileName +'.tsv'
outFile = open(outFilePath, 'w')

header = "\t" + 'nbSamples' + "\t" + 'maxActicity' + "\t" + 'sampleListStr' + "\t" + 'activityListStr' + "\n"
outFile.write(header)

# For each somatic source element
for srcId in somaticSrcActivityDict:
    sampleList = somaticSrcActivityDict[srcId].keys()
    activityList = somaticSrcActivityDict[srcId].values()

    nbSamples = len(sampleList)
    maxActicity = max(activityList)

    sampleListStr = ",".join(sampleList)
    activityListStr = ",".join(map(str, activityList))
    
    # Write row into the tsv
    row = srcId + "\t" + str(nbSamples) + "\t" + str(maxActicity) + "\t" + sampleListStr + "\t" + activityListStr + "\n"
    outFile.write(row)
   
outFile.close()


## End ##
print
print "***** Finished! *****"
print
