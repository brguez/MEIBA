#!/usr/bin/env python
#coding: utf-8 

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
    
def log(label, string):
    """ 
        Display labelled information 
    """ 
    print "[" + label + "]", string


def genotyper(bamFile, VCFlineObj, HETVAF, HOMVAF):
	'''
	
	'''
	### 1. Count the number of reads supporting the reference and MEI polymorphism
	nbReadsRef = 0
	nbReadsMEI = 0

	# Define genomic region +- (CIPOS + 5bp) around insertion breakpoint 
	# + 1 for the end interval position, reason: 
	# Convert from 1-based in VCF to 0-based half-open intervals in BAM [X,Y)
	# The position X is part of the interval, but Y is not.
	chrom = str(VCFlineObj.chrom)
	beg = int(VCFlineObj.pos) - (int(VCFlineObj.infoDict["CIPOS"]) + 5)  
	end = int(VCFlineObj.pos) + (int(VCFlineObj.infoDict["CIPOS"]) + 5) + 1

	subHeader("Genotype " + VCFlineObj.infoDict["CLASS"] + ":" + VCFlineObj.chrom + "_" + str(VCFlineObj.pos) + " with a CIPOS of " + VCFlineObj.infoDict["CIPOS"])

	# Extract alignments within the genomic region
	iterator = bamFile.fetch(chrom, beg, end)
	
	# Iterate over the alignments
	for alignment in iterator:

		# Discard unmapped reads
		if (alignment.is_unmapped == False):
			log("GENOTYPE", "Alignment info (bkpACoord, start, CIGAR): " + str(VCFlineObj.pos) + " " + str(alignment.reference_start) + " " + alignment.cigarstring)
	
			# Assess if the alignment supports the reference allele. Conditions:
			# a) Read overlap the insertion site with an overhang of 20 or more nucleotides on each side
			# b) Each segment properly aligned according to the aligner (bitwise flag 0x2)
			# c) Not PCR nor optical duplicated (bitwise flag 0x400)

			overhang = 20 
			overlap = alignment.get_overlap(VCFlineObj.pos - overhang, VCFlineObj.pos + overhang)
			properPair = alignment.is_proper_pair
			duplicate = alignment.is_duplicate
	
			log("GENOTYPE", "Reference allele info (overhang, properPair, duplicate): " + str(overlap) + " " + str(properPair) + " " + str(duplicate))
		
			# A) Read supporting the reference allele:
			if (overlap >= 40) and (properPair == True) and (duplicate == False): 
				nbReadsRef += 1			
				log("GENOTYPE", "Alignment supports REFERENCE")
			
			# B) Read do not supporting the reference allele   
			else:
				# Assess if the alignment supports the MEI polymorphism. Conditions:
				# a) Only beginning of the read soft (Operation=4) or hard clipped (Operation=5) 
				# b) Alignment start position within +- (CIPOS + 5bp) compared to the insertion breakpoint 
				# c) Not PCR nor optical duplicated (bitwise flag 0x400)
	
				firstOperation = alignment.cigartuples[0][0] 
				lastOperation = alignment.cigartuples[-1][0] 
		
				log("GENOTYPE", "MEI allele info (firstOperation, lastOperation, duplicate): " + str(firstOperation) + " " + str(lastOperation) + " " + str(duplicate))

				## Assess clipping
				# a) Clipping at the beginning of the read while not at the end
				if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation != 4) and (lastOperation != 5)): 
				
					## Assess breakpoint position
					beg = int(VCFlineObj.pos) - (int(VCFlineObj.infoDict["CIPOS"]) + 5) 
					end = int(VCFlineObj.pos) + (int(VCFlineObj.infoDict["CIPOS"]) + 5) 
	

					# bkp within range
					if (beg <= alignment.reference_start) and (alignment.reference_start <= end):
						bkp = True	
					else:
						bkp = False
	
				# b) Not clipped 
				else:
					bkp = False
					
				# Read supporting the MEI polymorphism:
				if (bkp == True) and (duplicate == False): 
					nbReadsMEI += 1			
					log("GENOTYPE", "Alignment supports MEI")
				 			
				
			print "--------"
	
	### 2. Compute the VAF
	totatNbReads = nbReadsMEI + nbReadsRef

	if totatNbReads == 0:
		vaf = 0
	else:
		vaf =  float(nbReadsMEI) / totatNbReads
	
	### 3. Determine donor genotype for the current variant
	# Note: Required at least three MEI supporting reads for considering the variant
		
	# A) Homozygous for MEI
	if (vaf >= HOMVAF) and (nbReadsMEI >= 3):
		genotype = '1/1'
        	
	# B) Heterozygous                
	elif (vaf >= HETVAF) and (nbReadsMEI >= 3):
		genotype = '0/1' 

	# C) Homozygous for reference	
	else:
		genotype = '0/0' 

	info("MEI genotype (genotype, nbReadsMEI, totatNbReads, VAF): " + genotype + " " + str(nbReadsMEI) + " " + str(totatNbReads) + " " + str(vaf) + "\n")	
	
	return (genotype, nbReadsMEI, totatNbReads, vaf)


#### CLASSES ####
#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import numpy as np
from operator import itemgetter, attrgetter, methodcaller
import pysam

## Get user's input ## 
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='...')
parser.add_argument('BAMlist', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
BAMlist = args.BAMlist
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

# Constants
HETVAF = 0.10  # Min VAF threshold for heterozygous (0/1)
HOMVAF = 0.9   # Min VAF threshold for homozygous alternative (1/1)

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "BAMlist: ", BAMlist
print "min-vaf-heterozygous: ", HETVAF
print "min-vaf-homozygous: ", HOMVAF
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 

## Start ## 

numDonors = sum(1 for line in open(BAMlist))
	
outFilePath = outDir + '/' + "genotyped.MEI.PCAWG." + str(numDonors) + ".vcf"

#### 1. Create VCF object and read input VCF

header("1. Process input VCF with polymorphic MEI for genotyping")

VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)


#### 2. Genotype donors
header("2. Genotype donors")

## Iterate over donor BAM files, genotyping one donor at each time
for line in open(BAMlist):

	BAMPath = line.rstrip('\n')
	
	base = os.path.basename(BAMPath)
	donorId = os.path.splitext(base)[0]

	subHeader("Donor: " + donorId)
	
	## Open donor's BAM file for reading
	bamFile = pysam.AlignmentFile(BAMPath, "rb")

	## Genotype each MEI polymorphism for the current donor
	for VCFlineObj in VCFObj.lineList:
		
		genotype, nbReadsMEI, totatNbReads, vaf = genotyper(bamFile, VCFlineObj, HETVAF, HOMVAF)

	bamFile.close()

sys.exit()
	
