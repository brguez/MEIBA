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


def worker(VCF, BAMset):
	'''
	'''

	threadId = threading.currentThread().getName()
	info( threadId + ' launched')
	time.sleep(random.uniform(0, 1))

	## Iterate over donor BAM files, genotyping one donor at a time
	for line in BAMset:
	
		BAMPath = line.rstrip('\n')
	
		base = os.path.basename(BAMPath)
		donorId = os.path.splitext(base)[0]

		subHeader(threadId +" genotyping donor " + donorId)
		time.sleep(random.uniform(0, 1))

		## Open donor's BAM file for reading
		bamFile = pysam.AlignmentFile(BAMPath, "rb")

		## Genotype each MEI polymorphism for the current donor
		for VCFlineObj in VCFObj.lineList:

			subHeader("Genotype " + VCFlineObj.infoDict["CLASS"] + ":" + VCFlineObj.chrom + "_" + str(VCFlineObj.pos) + " with a CIPOS of " + VCFlineObj.infoDict["CIPOS"])			
			chrom = str(VCFlineObj.chrom)

			## A) Estimate VAF for MEI breakpoint A
			# ------------------__TSD__*********TE******AAA__TSD__------------- Donor genome
			#		         bkpB		     bkpA
			# -------------------------------__TSD__--------------------------- Reference genome
			# 			      bkpA    bkpB
			info("Breakpoint A")
			bkpPos = int(VCFlineObj.pos)
			nbReadsMEI_A, totatNbReads_A, vaf_A = computeBkpVaf(bamFile, chrom, bkpPos, "A")

			## B) Estimate VAF for MEI breakpoint B 
			info("Breakpoint B")
			bkpPos = int(VCFlineObj.pos) + int(VCFlineObj.infoDict["TSLEN"])
			nbReadsMEI_B, totatNbReads_B, vaf_B = computeBkpVaf(bamFile, chrom, bkpPos, "B")
			
			## Print results
			print "TMP_VAF-A", nbReadsMEI_A, totatNbReads_A, vaf_A
			print "TMP_VAF-B", nbReadsMEI_B, totatNbReads_B, vaf_B

		subHeader("Finished " + donorId + " genotyping")
		bamFile.close()
	
	info( threadId + ' finished')
	

def computeBkpVaf(bamFile, chrom, bkpPos, bkpCat):
	'''
	
	'''
	### 1.Count the number of reads supporting the reference and MEI polymorphism
	nbReadsRef = 0
	nbReadsMEI = 0

	# Define genomic region +- 5bp around insertion breakpoint 
	# + 1 for the end interval position, reason: 
	# Convert from 1-based in VCF to 0-based half-open intervals in BAM [X,Y)
	# The position X is part of the interval, but Y is not.
	beg = bkpPos - 5  
	end = bkpPos + 6

	# Extract alignments overlapping the genomic region
	iterator = bamFile.fetch(chrom, beg, end)
	
	# Iterate over the alignments
	for alignment in iterator:

		# Discard unmapped reads
		if (alignment.is_unmapped == False):
			log("VAF-" + bkpCat, "Alignment info (bkpCoord, start, end, CIGAR): " + str(bkpPos) + " " + str(alignment.reference_start) + " " + str(alignment.reference_end) + " " + alignment.cigarstring)
	
			# Assess if the alignment supports the reference allele. Conditions:
			# a) Read overlap the insertion site with an overhang of 20 or more nucleotides on each side
			#    Note regarding overhang: it should be equal to the minimum anchor length for soft and hard clipped reads supporting the insertion. 
			#    Otherwise we are underestimating VAF values, since we are applying a more stringent criteria for considering reads supporting 
			#    the alternative than the reference allele...
			# b) Each segment properly aligned according to the aligner (bitwise flag 0x2)
			# c) Not PCR nor optical duplicated (bitwise flag 0x400)

			overhang = 20 
			overlap = alignment.get_overlap(bkpPos - overhang, bkpPos + overhang)
			properPair = alignment.is_proper_pair
			duplicate = alignment.is_duplicate
	
			log("VAF-" + bkpCat, "Reference allele info (overhang, properPair, duplicate): " + str(overlap) + " " + str(properPair) + " " + str(duplicate))
		
			# A) Read supporting the reference allele:
			if (overlap >= 40) and (properPair == True) and (duplicate == False): 
				nbReadsRef += 1			
				log("VAF-" + bkpCat, "Alignment supports REFERENCE")
			
			# B) Read do not supporting the reference allele   
			else:
				# Assess if the alignment supports the MEI polymorphism. Conditions:
				# a) bkpA: Only beginning of the read soft (Operation=4) or hard clipped (Operation=5) 
				#    bkpB: Only end of the read soft (Operation=4) or hard clipped (Operation=5)
				# b) Alignment start position within +- 5bp compared to the insertion breakpoint 
				# c) Not PCR nor optical duplicated (bitwise flag 0x400)
	
				firstOperation = alignment.cigartuples[0][0] 
				lastOperation = alignment.cigartuples[-1][0] 
		
				log("VAF-" + bkpCat, "MEI allele info (firstOperation, lastOperation, duplicate): " + str(firstOperation) + " " + str(lastOperation) + " " + str(duplicate))

				### Assess clipping
		
				## Define region to search for clipped reads	
				beg = bkpPos - 5 
				end = bkpPos + 5

				# A) Breakpoint A. Clipping at the beginning of the read while not at the end
				if (bkpCat == "A") and ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation != 4) and (lastOperation != 5)): 
				
					# a) clipping breakpoint within range (+-5 bkp in VCF)
					# beg <---------> end
					if (beg <= alignment.reference_start) and (alignment.reference_start <= end):
						bkpBool = True	

					# b) bkp not within range
					else:
						bkpBool = False

				# B) Breakpoint B. Clipping at the end of the read while not at the beginning
				elif (bkpCat == "B") and ((lastOperation == 4) or (lastOperation == 5)) and ((firstOperation != 4) and (firstOperation != 5)):

					# a) clipping breakpoint within range
					## beg <---------> end
					if (beg <= alignment.reference_end) and (alignment.reference_end <= end):
						bkpBool = True	

					# b) not within range					
					else:
						bkpBool = False

				# C) Not clipped or not properly clipped 
				else:
					bkpBool = False
	
				## Read supporting the MEI polymorphism:
				if (bkpBool == True) and (duplicate == False): 
					nbReadsMEI += 1			
					log("VAF-" + bkpCat, "Alignment supports MEI")		

	### 2. Compute the VAF
	totatNbReads = nbReadsMEI + nbReadsRef

	if totatNbReads == 0:
		vaf = 0
	else:
		vaf =  float(nbReadsMEI) / totatNbReads

	log("VAF-" + bkpCat, "VAF (nbReadsMEI, totatNbReads, VAF): " + str(nbReadsMEI) + " " + str(totatNbReads) + " " + str(vaf))

	return (nbReadsMEI, totatNbReads, vaf)

# -------------------------------------------------------
### 3. Determine donor genotype for the current variant
# Note: Required at least three MEI supporting reads for considering the variant
		
# A) Homozygous for MEI
#if (vaf >= HOMVAF) and (nbReadsMEI >= 3):
#genotype = '1/1'
        	
# B) Heterozygous                
#elif (vaf >= HETVAF) and (nbReadsMEI >= 3):
#genotype = '0/1' 

# C) Homozygous for reference	
#else:
#genotype = '0/0' 

#info("MEI genotype (genotype, nbReadsMEI, totatNbReads, VAF): " + genotype + " " + str(nbReadsMEI) + " " + str(totatNbReads) + " " + str(vaf) + "\n")			
# return (genotype, nbReadsMEI, totatNbReads, vaf)
# --------------------------------------------------

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
import threading
import random

## Get user's input ## 
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='VCF with the MEI polymorphism to genotype across the set of donors')
parser.add_argument('BAMPaths', help='Text file with the path to the donor BAM files. One BAM per row and donor')
parser.add_argument('-t', '--threads', default=1, dest='threads', type=int, help='Number of threads. Default: 1.' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
BAMPaths = args.BAMPaths
threads = args.threads
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

# Constants
HETVAF = 0.10  # Min VAF threshold for heterozygous (0/1)
HOMVAF = 0.9   # Min VAF threshold for homozygous alternative (1/1)

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "BAMPaths: ", BAMPaths
print "min-vaf-heterozygous: ", HETVAF
print "min-vaf-homozygous: ", HOMVAF
print "threads: ", threads
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 

## Start ## 


#### 1. Create VCF object and read input VCF

header("1. Process input VCF with polymorphic MEI for genotyping")

VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)

#### 2. Genotype donors
header("2. Genotype donors")

# Make list of input BAM files:
BAMList = [line.rstrip('\n') for line in open(BAMPaths)]

# Split donors into X evenly sized chunks. X = number of threads
chunkSize = int(round(len(BAMList) / float(threads)))

BAMChunks = [BAMList[i:i + chunkSize] for i in xrange(0, len(BAMList), chunkSize)]

threads = list()
counter = 1 

# Generate a genotyping thread per set of donors
for chunk in BAMChunks:

	print "chunk" + str(counter) + " : " + str(len(chunk)) + " donors to genotype" 
	
	threadName = "THREAD-" + str(counter)
	thread = threading.Thread(target=worker, args=(VCF, chunk), name=threadName)
	threads.append(thread)

	counter += 1

# Launch threads
[t.start() for t in threads]
	
# Wait till threads have finished
[t.join() for t in threads]

header("Finished")
