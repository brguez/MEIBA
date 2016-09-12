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


def worker(VCFObj, donorIdBamPathList, homVaf, hetVaf):
	'''
	'''

	threadId = threading.currentThread().getName()
	info( threadId + ' launched')
	time.sleep(random.uniform(0, 1))

	## Iterate over donor id/BAM files, genotyping one donor at a time
	for line in donorIdBamPathList:
	
		# Extract donorId and BAM file path 
		donorId = line[0]
		BAMPath = line[1]

		subHeader(threadId +" genotyping donor " + donorId)
		time.sleep(random.uniform(0, 1))

		## Open donor's BAM file for reading
		bamFile = pysam.AlignmentFile(BAMPath, "rb")

		## Genotype each MEI polymorphism for the current donor
		for VCFlineObj in VCFObj.lineList:
			genotype = genotypeMEI(bamFile, VCFlineObj, homVaf, hetVaf)
			VCFlineObj.genotypesDict[donorId] = genotype 

		subHeader("Finished " + donorId + " genotyping")
		bamFile.close()
	
	info( threadId + ' finished')
	
def genotypeMEI(bamFile, VCFlineObj, homVaf, hetVaf):
	'''
	'''

	subHeader("Genotype " + VCFlineObj.infoDict["CLASS"] + ":" + VCFlineObj.chrom + "_" + str(VCFlineObj.pos) + " with a CIPOS of " + VCFlineObj.infoDict["CIPOS"])			
	chrom = str(VCFlineObj.chrom)

	### 1. Compute VAF 
	## 1.1) Estimate VAF for MEI breakpoint A
	# ------------------__TSD__*********TE******AAA__TSD__------------- Donor genome
	#		         bkpB		     bkpA
	# -------------------------------__TSD__--------------------------- Reference genome
	# 			      bkpA    bkpB
	msg = "Breakpoint A VAF"
	if debugBool == True: info(msg)  
	bkpPos = int(VCFlineObj.pos)
	nbReadsMEI_A, totalNbReads_A, vaf_A = computeBkpVaf(bamFile, chrom, bkpPos, "A")

	## 1.2) Estimate VAF for MEI breakpoint B 
	msg = "Breakpoint B VAF"
	if debugBool == True: info(msg)	
	bkpPos = int(VCFlineObj.pos) + int(VCFlineObj.infoDict["TSLEN"])
	nbReadsMEI_B, totalNbReads_B, vaf_B = computeBkpVaf(bamFile, chrom, bkpPos, "B")
		
	## 1.3) Compute average VAF for the MEI
	msg = "Average VAF"
	if debugBool == True: info(msg)
	
	minNbClipped = 5

	## a) Both breakpoints supported by at least X clipped-reads 
	if (nbReadsMEI_A >= minNbClipped) and (nbReadsMEI_B >= minNbClipped): 
		nbReadsMEI = int(nbReadsMEI_A + nbReadsMEI_B) / 2 
		totalNbReads = int(totalNbReads_A + totalNbReads_B) / 2
		vaf = float(nbReadsMEI) / totalNbReads

	## b) Breakpoint A supported by at least X clipped-reads
	elif (nbReadsMEI_A >= minNbClipped):
		nbReadsMEI = nbReadsMEI_A 
		totalNbReads = totalNbReads_A
		vaf = vaf_A

	## c) Breakpoint B supported by at least X clipped-reads
	elif (nbReadsMEI_B >= minNbClipped):
		nbReadsMEI = nbReadsMEI_B 
		totalNbReads = totalNbReads_B
		vaf = vaf_B
 
	## d) No breakpoint supported by at least X clipped-reads 
	else:
		nbReadsMEI = int(nbReadsMEI_A + nbReadsMEI_B) / 2 
		totalNbReads =	int(totalNbReads_A + totalNbReads_B) / 2
		vaf = 0

	### 2. Determine donor genotype for the current variant
	# a) Homozygous for MEI
	if (vaf >= homVaf):
		genotype = '1/1'
        	
	# b) Heterozygous for MEI             
	elif (vaf >= hetVaf):
		genotype = '0/1' 

	# c) Homozygous for reference	
	else:
		genotype = '0/0' 

	msg = "MEI genotype (genotype, nbReadsMEI, totalNbReads, VAF): " + genotype + " " + str(nbReadsMEI) + " " + str(totalNbReads) + " " + str(vaf) + "\n"
	if debugBool == True: info(msg)

	## Make genotype field and return output
	genotypeField = genotype + ':' + str(nbReadsMEI) + ':' + str(totalNbReads)
	
	return (genotypeField)


def computeBkpVaf(bamFile, chrom, bkpPos, bkpCat):
	'''
	'''
	tag = "VAF-" + bkpCat

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
			msg = "Alignment info (bkpCoord, start, end, CIGAR): " + str(bkpPos) + " " + str(alignment.reference_start) + " " + str(alignment.reference_end) + " " + alignment.cigarstring
			if (debugBool == True): log(tag, msg)  					

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
	
			msg = "Reference allele info (overhang, properPair, duplicate): " + str(overlap) + " " + str(properPair) + " " + str(duplicate)
			if (debugBool == True): log(tag, msg) 	

			# A) Read supporting the reference allele:
			if (overlap >= 40) and (properPair == True) and (duplicate == False): 
				nbReadsRef += 1	
				msg = "Alignment supports REFERENCE"
				if (debugBool == True): log(tag, msg) 	

			# B) Read do not supporting the reference allele   
			else:
				# Assess if the alignment supports the MEI polymorphism. Conditions:
				# a) bkpA: Only beginning of the read soft (Operation=4) or hard clipped (Operation=5) 
				#    bkpB: Only end of the read soft (Operation=4) or hard clipped (Operation=5)
				# b) Alignment start position within +- 3bp compared to the insertion breakpoint 
				# c) Not PCR nor optical duplicated (bitwise flag 0x400)
	
				firstOperation = alignment.cigartuples[0][0] 
				lastOperation = alignment.cigartuples[-1][0] 
			
				msg = "MEI allele info (firstOperation, lastOperation, duplicate): " + str(firstOperation) + " " + str(lastOperation) + " " + str(duplicate)
				if (debugBool == True): log(tag, msg) 	

				### Assess clipping
		
				## Define region to search for clipped reads	
				beg = bkpPos - 3
				end = bkpPos + 3

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

					msg = "Alignment supports MEI"
					if (debugBool == True): log(tag, msg) 
		
		msg = "-------------------"
		if (debugBool == True): print msg 

	### 2. Compute the VAF
	totalNbReads = nbReadsMEI + nbReadsRef

	if totalNbReads == 0:
		vaf = 0
	else:
		vaf =  float(nbReadsMEI) / totalNbReads

	msg = "VAF (nbReadsMEI, totalNbReads, VAF): " + str(nbReadsMEI) + " " + str(totalNbReads) + " " + str(vaf)
	if (debugBool == True): log(tag, msg) 

	return (nbReadsMEI, totalNbReads, vaf)


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

# Global variables:
global debugBool ## debug logging mode. Boolean.

## Get user's input ## 
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='VCF with the MEI polymorphism to genotype across the set of donors')
parser.add_argument('BAMPaths', help='Tab separated text file with two columns: 1) donor ids and 2) path to the corresponding donor normal BAM files. One donor per row.')
parser.add_argument('sampleId', help='Identifier to name output file.')
parser.add_argument('--hetVaf', default=0.10, dest='hetVaf', type=float, help='Min VAF threshold for heterozygous (0/1). Default: 0.1' )
parser.add_argument('--homVaf', default=0.9, dest='homVaf', type=float, help='Min VAF threshold for homozygous alternative (1/1). Default: 0.9' )
parser.add_argument('-t', '--threads', default=1, dest='threads', type=int, help='Number of threads. Default: 1' )
parser.add_argument('--debug', action='store_true', dest='debug', help='Enable debug mode. Display detailed log info about MEI genotyping')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
BAMPaths = args.BAMPaths
sampleId =  args.sampleId
hetVaf = args.hetVaf
homVaf = args.homVaf
threads = args.threads
debugBool = args.debug
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "BAMPaths: ", BAMPaths
print "sampleId: ", sampleId
print "min-vaf-heterozygous: ", hetVaf
print "min-vaf-homozygous: ", homVaf
print "threads: ", threads
print "debug: ", debugBool
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

# Make list of donorId/BAM files:
donorIdBamPathList = [line.split() for line in open(BAMPaths)]

# Split donors into X evenly sized chunks. X = number of threads
chunkSize = int(round(len(donorIdBamPathList) / float(threads)))

BAMChunks = [donorIdBamPathList[i:i + chunkSize] for i in xrange(0, len(donorIdBamPathList), chunkSize)]

threads = list()
counter = 1 

# Generate a genotyping thread per set of donors
for chunk in BAMChunks:

	print "chunk" + str(counter) + ": " + str(len(chunk)) + " donors to genotype" 
	
	threadName = "THREAD-" + str(counter)
	thread = threading.Thread(target=worker, args=(VCFObj, chunk, homVaf, hetVaf), name=threadName)
	threads.append(thread)

	counter += 1

# Launch threads
[t.start() for t in threads]
	
# Wait till threads have finished genotyping
[t.join() for t in threads]

#### 3. Make output multisample VCF file
header("3. Produce multi-sample VCF file as ouput")

fileName = sampleId + '.vcf'
outFilePath = outDir + '/' + fileName

# 3.1 Write header
VCFObj.write_header(outFilePath)

# 3.2 Write variants
donorIdList = [line[0] for line in donorIdBamPathList] 

VCFObj.write_variants_multiSample(donorIdList, outFilePath)


header("Finished")
