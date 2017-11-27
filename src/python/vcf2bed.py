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


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import time
import formats

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('vcf', help='VCF containing MEI')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('--overHang', default=500, dest='overHang', type=int, help='Bed coordinates will be defined as (beg-overHang) and (end+overHang). Default: 500' )
parser.add_argument('--filter', action='store_true', dest='filter', help='Flag to filter out insertions that do have PASS filtering status.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.vcf
fileName = args.fileName
overHang = args.overHang
filterBool = args.filter
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "fileName: ", fileName
print "overHang: ", overHang
print "filterBool: ", filterBool
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Read input VCF and generate a VCF object
#################################################
header("1. Process VCF as input")

VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

#### 2. Write each MEI as a BED data line
##########################################
header("2. Write each MEI as a BED data line")

# Open output file
outFilePath = outDir + '/' + fileName + '.bed'
outFile = open(outFilePath, 'w')

# Write header:
row = '#chrom' + "\t" + 'beg' +  "\t" + 'end' + "\t" + 'iClass' + "\t" + 'iType' + "\t" + 'srcId' + '\t' + "GR" + "\n"
outFile.write(row)

## For each MEI
for MEIObj in VCFObj.lineList:

    ## If filter enabled select only those MEI that passes all the filters
    # If not enabled select all the MEI
    if (filterBool == False) or (filterBool == True and MEIObj.filter == "PASS"):

        chrom = MEIObj.chrom
        beg = str(MEIObj.pos - overHang)
        end = str(int(MEIObj.infoDict['BKPB']) + overHang) if 'BKPB' in MEIObj.infoDict else str(MEIObj.pos + overHang)
        iClass = MEIObj.infoDict["CLASS"]
        iType = MEIObj.infoDict["TYPE"]
        srcId = MEIObj.infoDict['SRCID'] if 'SRCID' in MEIObj.infoDict else 'NA' 
        GR = MEIObj.infoDict['GR'] if 'GR' in MEIObj.infoDict else 'NA' 
        row = chrom + '\t' + beg + '\t' + end + '\t' + iClass + '\t' + iType + '\t' + srcId + '\t' + GR + '\n'
        outFile.write(row)
       

#### End
header("FINISH!!")


