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
    print timeInfo, string,

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

# Constants
inputVCF = "/Users/brodriguez/Research/Scripts/Bash/TEIBA/test/0.3.3/testing.vcf"
inputBAM = "/Users/brodriguez/Research/Scripts/Bash/TEIBA/test/0.3.3/testing.bam"

HETVAF = 0.10  # Min VAF threshold for heterozygous (0/1)
HOMVAF = 0.9   # Min VAF threshold for homozygous alternative (1/1)
    
## Start ## 

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## 2. Open BAM file for reading:

bamfile = pysam.AlignmentFile(inputBAM, "rb")
print bamfile

## 3. Genotype each MEI polymorphism 

for VCFlineObj in VCFObj.lineList:
	
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

	print "insertion_coord", VCFlineObj.chrom, VCFlineObj.pos, VCFlineObj.infoDict["CIPOS"] 

	# Extract alignments within the genomic region
	iterator = bamfile.fetch(chrom, beg, end)
	
	# Iterate over the alignments
	for alignment in iterator:
		print "alignment: ", alignment, alignment.reference_start

		# Assess if the alignment supports the reference allele. Conditions:
		# a) Read overlap the insertion site with an overhang of 20 or more nucleotides on each side
		# b) Each segment properly aligned according to the aligner (bitwise flag 0x2)
		# c) Not PCR nor optical duplicated (bitwise flag 0x400)
		# d) Not secondary alignment (bitwise flag 0x100)

		overhang = 20 
		overlap = alignment.get_overlap(VCFlineObj.pos - overhang, VCFlineObj.pos + overhang)
		properPair = alignment.is_proper_pair
		duplicate = alignment.is_duplicate
		secondary = alignment.is_secondary

		print "proterties: ", overhang, properPair, duplicate, secondary

		# A) Read supporting the reference allele:
		if (overlap >= 40) and (properPair == True) and (duplicate == False) and (secondary == False): 
			nbReadsRef += 1			
			print "REFERENCE"
		
		# B) Read do not supporting the reference allele   
		else:
			# Assess if the alignment supports the MEI polymorphism. Conditions:
			# a) Only beginning of the read soft (Operation=4) or hard clipped (Operation=5) 
			# b) Alignment start position within +- (CIPOS + 5bp) compared to the insertion breakpoint 
			# c) Not PCR nor optical duplicated (bitwise flag 0x400)

			firstOperation = alignment.cigartuples[0][0] 
			lastOperation = alignment.cigartuples[-1][0] 

			print "operations: ", firstOperation, lastOperation

			## Assess clipping
			# a) Clipping at the beginning of the read while not at the end
			if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation != 4) and (lastOperation != 5)): 
		
				## Assess breakpoint position
				beg = int(VCFlineObj.pos) - (int(VCFlineObj.infoDict["CIPOS"]) + 5) 
				end = int(VCFlineObj.pos) + (int(VCFlineObj.infoDict["CIPOS"]) + 5) 

				print "interval: ", beg, end, alignment.reference_start

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
				print "MEI"	
			 			
			
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
		genotype = '1/1:' + str(nbReadsMEI) + ':' + str(totatNbReads)
        
	# B) Heterozygous                
	elif (vaf >= HETVAF) and (nbReadsMEI >= 3):
		genotype = '0/1:' + str(nbReadsMEI) + ':' + str(totatNbReads)

	# C) Homozygous for reference	
	else:
		genotype = '0/0:' + str(nbReadsMEI) + ':' + str(totatNbReads)
  
	print "results: ", nbReadsMEI, nbReadsRef, totatNbReads, vaf, genotype
	

	print "*****************"



