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
import formats
import time
import scipy.stats as stats
import numpy as np
from matplotlib import pyplot as plt

## Get user's input ##
parser = argparse.ArgumentParser(description= "Filter out genotyped MEI according to two criteria: allele count and % missing genotypes")
parser.add_argument('inputVCF', help='multi-sample VCF file containing genotyped MEI')
parser.add_argument('metadata', help='PCAWG donors metadata file.')
parser.add_argument('fileName', help='Identifier to name the output file.')
parser.add_argument('--minAlleleCount', default=1, dest='minAlleleCount', type=int, help='Minimum allele count. Default: 1' )
parser.add_argument('--maxPercMissingGt', default=5, dest='maxPercMissingGt', type=int, help='Maximum percentage of missing genotypes. Default: 5' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
metadata = args.metadata
fileName =  args.fileName
minAlleleCount = args.minAlleleCount
maxPercMissingGt = args.maxPercMissingGt
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "metadata: ", metadata
print "fileName: ", fileName
print "minAlleleCount: ", minAlleleCount
print "maxPercMissingGt: ", maxPercMissingGt
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 0. Create dictionary with donor id gender equivalencies
#############################################################
header("0. Create dictionary with donor id gender equivalencies")

metadata = open(metadata, 'r')
genderDict = {}

# Read file line by line
for line in metadata:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        donorId = fieldsList[0]
        gender = fieldsList[5]

        genderDict[donorId] = gender

#### 1. Read input multi-sample VCF and generate a VCF object
###############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)


#### 2. Apply allele count and missing genotypes filters
############################################################
header("2. Apply allele count and missing genotypes filters")

### Make table containing for each germline insertions its VAF
## Two columns:
# insertionsId(chrom:pos)   VAF

# Open output file
outFilePath = outDir + '/' + fileName + '.tsv'
outFile = open(outFilePath, 'w')

# Write header:
row = '#MEI' + "\t" + 'alleleCount' + "\t" + 'nbMissingGt' + "\t" + 'percMissing' + "\t" + "\n"
outFile.write(row)

#header("3. Compute parameters")

totalNbGt = 0
totalNbMissingGt = 0

for MEIObj in VCFObj.lineList:

    ## Compute MEI allele count, genotyping VAF and update counters to heterozygous and homozygous MEI per donor
    # MEI allele count: number of chromosomes in the population harvouring the MEI
    # Genotyping VAF: Ratio (number_reads_supporting_MEI)/(total_nb_reads_covering_MEI)
    chrom = MEIObj.chrom
    alleleCount = 0
    nbGt = 0
    nbMissingGt = 0

    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]
        gender = genderDict[donorId]

        ## Update counters and store VAF values
        # a) Homozygous alternative        
        if (genotype == "1/1"): 
            nbGt += 1
            totalNbGt += 1
            alleleCount += 2

        # b) Heterozygous or haploid carrier (for male variants in the X or Y chromosomes outside the PAR region)
        elif (genotype == "0/1") or (genotype == "1"):
            nbGt += 1
            totalNbGt += 1
            alleleCount += 1

        # c) Homozygous reference or haploid not carrier (for male variants in the X or Y chromosomes outside the PAR region)        
        elif (genotype == "0/0") or (genotype == "0"):
            nbGt += 1
            totalNbGt += 1

        # d) Missing genotype (diploid (./.) or haploid (.))    
        else:

            ## Variant outside the Y chromosome or in the Y but in males (not consider variants in the Y for females as they do not have)
            if (chrom != "Y") or ((chrom == "Y") and (gender == "male")):
                nbGt += 1
                totalNbGt += 1 
                totalNbMissingGt += 1 
                nbMissingGt += 1

    percMissing = float(nbMissingGt)/float(nbGt)*100

    # a) Both filters failed (Allele count < minAlleleCount and % missing genotypes > maxPercMissingGt)
    if (alleleCount < minAlleleCount) and (percMissing > maxPercMissingGt):
        MEIObj.filter = "COUNT;MISSGT"
    
    # b) Allele count < minAlleleCount
    elif (alleleCount < minAlleleCount):
        MEIObj.filter = "COUNT"

    # c) Percentage of missing genotypes > 5
    elif (percMissing > maxPercMissingGt):
        MEIObj.filter = "MISSGT"

    # d) Passed all the filters 
    else:
        MEIObj.filter = "PASS"

    ## Add row to the table
    row = MEIObj.chrom + ":" + str(MEIObj.pos) + ":" + MEIObj.infoDict['CLASS'] + "\t" + str(alleleCount) + "\t" + str(nbMissingGt) + "\t" + str(percMissing) + "\n"
    outFile.write(row)

percMissingTotal = float(totalNbMissingGt)/float(totalNbGt)*100

print "STATS: ", totalNbGt, totalNbMissingGt, percMissingTotal 

#### 3. Write output VCF file with filtered VCF
header("3. Write output VCF file")

outFilePath = outDir + '/' + fileName + '.vcf'

# 1. Write header
VCFObj.write_header(outFilePath)

# 2. Write variants
VCFObj.write_variants_multiSample(donorIdList, outFilePath)



####
header("Finished")
