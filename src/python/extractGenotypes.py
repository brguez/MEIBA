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
parser = argparse.ArgumentParser(description= "Make a text file containing genotyping information (donorId, genotype) per target MEI")
parser.add_argument('inputVCF', help='Multi-sample VCF file containing genotyped MEI')
parser.add_argument('targetMEI', help='Text file in tabular format containing one line per MEI of interest and three columns: chromosome, position and MEI class')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
targetMEI =  args.targetMEI
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "targetMEI: ", targetMEI
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


#### 2. Make a text file with genotyping information for each target MEI
########################################################################
header("2. Make a text file with genotyping information for each target MEI")
targetMEI = open(targetMEI, 'r')


## For each target MEI
for line in targetMEI:
	
	# Skip header	
	if not line.startswith("#"):
		
		# Make variables with needed info
		line = line.rstrip('\n')
	        line = line.split('\t')
		chrom = line[0] 
		pos = line[1]
		category = line[2]
		gene = line[3]
	
		# Open output file
		outFilePath = outDir + '/' + 'chr' + chrom + '_' + pos + '_' + category + '_genotypes.txt' 
		outFile = open(outFilePath, 'a')

		# Write header:
		row = '#donorId' + "\t" + 'genotype' + "\n"
		outFile.write(row)

		# Iterate over the MEI in the VCF
		for MEIobj  in VCFObj.lineList:
			
			# Current MEI is the target MEI
			if (chrom == str(MEIobj.chrom)) and (pos == str(MEIobj.pos)) and (category == MEIobj.infoDict["CLASS"]):
	 
				# Print into the output file a row per donor with its donorId and genotype for the current MEI
				for donorId in donorIdList:
					genotype = MEIobj.genotypesDict[donorId].split(':')[0]

					row = donorId + "\t" + genotype + "\n"
					outFile.write(row)
				
				# Stop iterating once target MEI found 
				break

	    
#### END
header("Finished")
