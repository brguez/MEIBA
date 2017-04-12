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
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputVCF', help='multi-sample VCF file containing genotyped MEI')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

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
outFilePath = outDir + '/MEI_alleleCount_nbMissingGt_percMissing.tsv'
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
    alleleCount = 0
    nbGt = 0    
    nbMissingGt = 0
    
#    print "MEI: ", MEIObj.chrom, MEIObj.pos, MEIObj.infoDict["CLASS"]

    nbGt = len(MEIObj.genotypesDict)

    for donorId, genotypeField in MEIObj.genotypesDict.iteritems():

        genotypeFieldList = genotypeField.split(":")
        genotype = genotypeFieldList[0]

        #print donorId, genotype

        ## Update counters and store VAF values
        # a) Homozygous alternative        
        if (genotype == "1/1"): 
            totalNbGt += 1
            alleleCount += 2

        # b) Heterozygous or haploid carrier (for male variants in the X or Y chromosomes outside the PAR region)
        elif (genotype == "0/1") or (genotype == "1"):
            totalNbGt += 1
            alleleCount += 1

        # c) Homozygous reference or haploid not carrier (for male variants in the X or Y chromosomes outside the PAR region)        
        elif (genotype == "0/0") or (genotype == "0"):
            totalNbGt += 1

        # d) Missing genotype (diploid (./.) or haploid (.))    
        else:
            totalNbGt += 1 
            totalNbMissingGt += 1 
            nbMissingGt += 1

    percMissing = float(nbMissingGt)/float(nbGt)*100

    # a) Both filters failed
    if (alleleCount == 0) and (percMissing > 5):
        MEIObj.filter = "COUNT;MISSGT"
    
    # b) Allele count of 0
    elif (alleleCount == 0):
        MEIObj.filter = "COUNT"

    # c) Percentage of missing genotypes > 5
    elif (percMissing > 5):
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

outFilePath = outDir + '/' + "germline_MEI_genotyped_normal_PCAWG.filtered.vcf"

# 1. Write header
VCFObj.write_header(outFilePath)

# 2. Write variants
VCFObj.write_variants_multiSample(donorIdList, outFilePath)



####
header("Finished")
