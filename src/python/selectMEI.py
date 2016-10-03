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
parser = argparse.ArgumentParser(description= "Select MEI from VCF based on a list of target MEIs. Report VCF with those MEI from initial VCF included in the list of target MEIs")
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

#### 2. Select MEI from VCF based on a list of target MEIs
############################################################
header("2. Select MEI from VCF based on a list of target MEIs")
targetMEI = open(targetMEI, 'r')

outVCFObj = formats.VCF()

counter = 1

## For each target MEI
for line in targetMEI:

	msg = "Assessing if there is a MEI in the input VCF overlapping target MEI: " + str(counter)
	info(msg)	
	
	# Skip header	
	if not line.startswith("#"):
		
		# Make variables with needed info
		line = line.rstrip('\n')
	        line = line.split('\t')
		chrom = line[0] 
		pos = int(line[1])
		category = line[2]

		# Iterate over the MEI in the VCF
		for MEIObj  in VCFObj.lineList:
			MEItype = MEIObj.alt.split(":")[2][:-1]
			
			# Use our MEI type naming convention (1000 genomes uses a different one..) 
			if (MEItype == "LINE1"):
				MEItype = "L1"

			elif (MEItype == "ALU"):
				MEItype = "Alu"

			else:
				MEItype = "SVA"

			# Current MEI in the same chromosome and of the same type as the target MEI
			if (chrom == str(MEIObj.chrom)) and (category == MEItype):
				
				# Define range +-20 bp to search for overlap between target and MEI in VCF breakpoints			
				beg = pos - 20
				end = pos + 20 

				# Current MEI overlapping target MEI breakpoint range
				if (int(MEIObj.pos) >= beg) and (int(MEIObj.pos) <= end):
					
					# Add overlapping MEI object
					outVCFObj.addLine(MEIObj)

					# Stop iterating when match found
					break

				# Current MEI with a breakpoint with a position higher than end range -> stop iterating				
				if (int(MEIObj.pos) >= end):
					break
	counter+=1

#### 3. Make output multisample VCF file 
############################################
header("3. Produce multi-sample VCF with the selected MEI as ouput ")

fileName = 'selectedMEI.genotyped.vcf'
outFilePath = outDir + '/' + fileName

# 3.1 Write header
outVCFObj.header = VCFObj.header
outVCFObj.write_header(outFilePath)

# 3.2 Write variants

outVCFObj.write_variants_multiSample(donorIdList, outFilePath)


#### END
header("Finished")
