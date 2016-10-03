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

## Get user's input ## 
parser = argparse.ArgumentParser(description= "Make table with number of homozygous and heterozygous donors per MEI from a multi-sample VCF file")
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
#############################################################
header("1. Process multi-sample VCF as input")

VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(inputVCF)

#### 2. Count the number of homozygous and heterozygous donors per MEI
#########################################################################
header("2. Count the number of homozygous and heterozygous donors per MEI")

# Open output file
outFilePath = outDir + '/nbDonors_perMEI_statistics.txt' 
outFile = open(outFilePath, 'a')

# Write header:
row = '#MEI' + "\t" + 'totalNb' + "\t" + 'nbHeteroz' + "\t" + 'nbHomoz' + "\n"
outFile.write(row)


for MEIObj in VCFObj.lineList:
	
	nbTotalDonors = 0
	nbHeterozygousDonors = 0 
	nbHomozygousDonors =  0

	MEIid =   MEIObj.chrom + '_'  + str(MEIObj.pos)

 	for donorId, genotypeField in MEIObj.genotypesDict.iteritems():		

		genotypeFieldList = genotypeField.split(":")
		genotype = genotypeFieldList[0]

		# a) Heterozygous
		if (genotype == "0/1"):
			nbTotalDonors += 1
			nbHeterozygousDonors += 1
	
		# b) Homozygous		
		elif (genotype == "1/1"):
			nbTotalDonors += 1
			nbHomozygousDonors += 1

	# Write into output file 	
	row = MEIid + "\t" + str(nbTotalDonors) + "\t" + str(nbHeterozygousDonors) + "\t" + str(nbHomozygousDonors) + "\n" 
	outFile.write(row)

#### Finish:
header("Finished")
