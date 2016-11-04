#!/usr/bin/env python
#coding: utf-8

def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='...')
parser.add_argument('driverDb', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
driverDb = args.driverDb
donorId = args.donorId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "driver-database: ", driverDb
print "donorId: ", donorId
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

outFilePath = outDir + '/' + donorId + ".driverAnnot.vcf"

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## 2. Create dictionary containing driver genes
cancerDriversDict = {}

# Open driver database input file
driverDbFile = open(driverDb, 'r')

## Read first input file line by line
for line in driverDbFile:
    line = line.rstrip('\n')

    ## Skip header
    if not line.startswith("#"):
        line = line.split('\t')

        ## Save information into the dictionary:
        gnSymbol = line[0]
        role = line[1]
        cosmic = line[2]
        cpg = line[3]

        # Create driver annotation dictionary
        annotDict = {}
        annotDict["ROLE"] = role
        annotDict["COSMIC"] = cosmic
        annotDict["CPG"] = cpg

        # Add driver annotation dictionary as value
        cancerDriversDict[gnSymbol] = annotDict

## 3. Add driver annotation info to VCF file
for VCFlineObj in VCFObj.lineList:

        # Insertion overlapping a Gene:
    if "GENE" in VCFlineObj.infoDict:
        overlappingGn = VCFlineObj.infoDict["GENE"]

        ## Gene in cancer driver genes database:
        if overlappingGn in cancerDriversDict:

            ## Add cancer driver annotation to the VCF info field
            VCFlineObj.infoDict.update(cancerDriversDict[overlappingGn])

            VCFlineObj.info = VCFlineObj.make_info()

## 3. Make output VCF

# 3.1 Write header
VCFObj.write_header(outFilePath)

# 3.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print
print "***** Finished! *****"
print
