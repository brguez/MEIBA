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
parser = argparse.ArgumentParser(description= "Select MEI from VCF based on a list of target MEIs in bed format. Report VCF with those MEI from initial VCF included in the bed.")
parser.add_argument('inputVCF', help='Uni-sample or multi-sample VCF file containing MEI')
parser.add_argument('inputBed', help='Bed file containing one line per MEI of interest and four columns: chromosome, beg, end and MEI class.')
parser.add_argument('sampleId', help='Identifier to name output file.')
parser.add_argument('-multiSample', action='store_true')
parser.add_argument('--overhang', default=0, type=int, dest='overhang', help='Maximum overhang for searching overlap between MEI VCF and Bed coordinates. Default: 0 base pairs.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='Output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.inputVCF
inputBed =  args.inputBed
sampleId =  args.sampleId
multiSample = args.multiSample
overhang = args.overhang
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputVCF: ", inputVCF
print "inputBed: ", inputBed
print "sampleId: ", sampleId
print "multiSample: ", multiSample
print "overhang: ", overhang
print 
print "***** Executing ", scriptName, ".... *****"
print 

## Start ## 

#### 1. Read input VCF and generate a VCF object
#############################################################
header("1. Process input VCF" )

VCFObj = formats.VCF()

# a) One sample VCF file
if (multiSample == False):
	print "read-onesample"
	VCFObj.read_VCF(inputVCF)

# b) Multi-sample VCF file
else:
	print "read-multisample"
	donorIdList = VCFObj.read_VCF_multiSample(inputVCF)


#### 2. Select MEI from VCF based on a list of target MEIs
############################################################
header("2. Select MEI from VCF based on a list of target MEIs")
inputBed = open(inputBed, 'r')

outVCFObj = formats.VCF()

counter = 1

## For each target MEI
for line in inputBed:

	msg = "Assessing if there is a MEI in the input VCF overlapping target MEI: " + str(counter)
	info(msg)	
	
	# Skip header	
	if not line.startswith("#"):
		
		# Make variables with needed info
		line = line.rstrip('\n')
	        line = line.split('\t')
		chrom = line[0] 
		beg = int(line[1]) - overhang
		end = int(line[2]) + overhang
		category = line[3]

		print 'bed_MEI:', chrom, beg, end, category

		# Iterate over the MEI in the VCF
		for MEIObj  in VCFObj.lineList:
			MEItype = MEIObj.infoDict['CLASS']
 
			# Current MEI in the same chromosome and of the same type as the target MEI
			if (chrom == str(MEIObj.chrom)) and (category == MEItype):
			
				# Current MEI overlapping target MEI breakpoint range
				if (int(MEIObj.pos) >= beg) and (int(MEIObj.pos) <= end):
		
					print 'vcf_MEI: ', MEIObj.chrom, MEIObj.pos, MEItype				
		
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

fileName = sampleId + '.vcf'
outFilePath = outDir + '/' + fileName

## 3.1 Write header
outVCFObj.header = VCFObj.header
outVCFObj.write_header(outFilePath)

## 3.2 Write variants
# a) One sample VCF file
if (multiSample == False):
	print "write-onesample"
	outVCFObj.write_variants(outFilePath)

# b) Multi-sample VCF file
else:
	print "write-multisample"
	outVCFObj.write_variants_multiSample(donorIdList, outFilePath)


#### END
header("Finished")
