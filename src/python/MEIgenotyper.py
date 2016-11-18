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


def worker(VCFObj, donorIdBamPathList, maxDist, minClipped, minMapQual, minREF, hetVaf, homVaf):
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
            genotype = genotypeMEI(bamFile, VCFlineObj, maxDist, minClipped, minMapQual, minREF, hetVaf, homVaf)
            VCFlineObj.genotypesDict[donorId] = genotype

        subHeader("Finished " + donorId + " genotyping")
        bamFile.close()

    info( threadId + ' finished')

def genotypeMEI(bamFile, VCFlineObj, maxDist, minClipped, minMapQual, minREF, hetVaf, homVaf):
    '''
    '''

    subHeader("Genotype " + VCFlineObj.infoDict["CLASS"] + ":" + VCFlineObj.chrom + "_" + str(VCFlineObj.pos) + " with a CIPOS of " + VCFlineObj.infoDict["CIPOS"])
    chrom = str(VCFlineObj.chrom)

    ### 1. Compute VAF
    ## 1.1) Estimate VAF for MEI breakpoint A
    msg = "Breakpoint A VAF"
    if debugBool == True: info(msg)
    bkpPos = int(VCFlineObj.pos)
 
    ## a) Insertion absent in reference genome
    # ------------------__TSD__*********TE******AAA__TSD__------------- Donor genome
    #                          bkpB                bkpA
    # -------------------------------__TSD__--------------------------- Reference genome
    #                                bkpA   bkpB
    # * Clipping at the beginning of the read
    if (VCFlineObj.alt == "<MEI>"):

        nbReadsALT_A, totalNbReads_A, vaf_A = computeBkpVaf(bamFile, chrom, bkpPos, "Beg", maxDist, minMapQual)

    ## b) Insertion in reference genome and absent in donor genome
    # -------------------------------__TSD__--------------------------- Donor genome
    #                                bkpB   bkpA                           
    # ------------------__TSD__*********TE******AAA__TSD__------------- Reference genome genome
    #                          bkpA                bkpB
    # * Clipping at the ending of the read
    elif (VCFlineObj.ref == "<MEI>"):
        
        nbReadsALT_A, totalNbReads_A, vaf_A = computeBkpVaf(bamFile, chrom, bkpPos, "End", maxDist, minMapQual)

    ## c) Raise error...  
    else:
        msg="Incorrectly formated VCF line"
        info(msg)

    ## 1.2) Estimate VAF for MEI breakpoint B if possible
    # A) Second breakpoint information available
    if ('BKPB' in VCFlineObj.infoDict):

        msg = "Breakpoint B VAF"
        if debugBool == True: info(msg)
        bkpPos = int(VCFlineObj.infoDict["BKPB"])

        ## a) Insertion absent in reference genome
        # * Clipping at the ending of the read
        if (VCFlineObj.alt == "<MEI>"):
        
            nbReadsALT_B, totalNbReads_B, vaf_B = computeBkpVaf(bamFile, chrom, bkpPos, "End", maxDist, minMapQual)

        ## b) Insertion in reference genome and absent in donor genome
        # * Clipping at the beginning of the read
        elif (VCFlineObj.ref == "<MEI>"):

            nbReadsALT_B, totalNbReads_B, vaf_B = computeBkpVaf(bamFile, chrom, bkpPos, "Beg", maxDist, minMapQual)

        ## c) Raise error...  
        else:
            msg="Incorrectly formated VCF line"
            info(msg)

        ### Compute average VAF for the MEI (bkpA and bkpB)
        msg = "Average VAF"
        if debugBool == True: info(msg)
        nbReadsALT, totalNbReads, vaf = computeAverageVaf(nbReadsALT_A, totalNbReads_A, vaf_A, nbReadsALT_B, totalNbReads_B, vaf_B, minClipped)

    # B) Breakpoint B information not available (use bkp A information)
    else:
        nbReadsALT = nbReadsALT_A
        totalNbReads = totalNbReads_A
        vaf = vaf_A
        
    ### 2. Determine donor genotype for the current variant
    nbReadsREF = totalNbReads - nbReadsALT

    # a) Homozygous for alternative allele
    if (vaf >= homVaf):
        genotype = '1/1'

    # b) Heterozygous for alternative
    elif (vaf >= hetVaf):
        genotype = '0/1'
        
    # c) Homozygous for reference (at least 5 reference supporting reads required by default)
    elif (nbReadsREF >= minREF):
        genotype = '0/0'

    # d) Missing genotype (Not enough number of reads supporting alternative nor reference alleles)
    else:
        genotype = './.'        

    msg = "MEI genotype (genotype, nbReadsALT, totalNbReads, VAF): " + genotype + " " + str(nbReadsALT) + " " + str(totalNbReads) + " " + str(vaf) + "\n"
    if debugBool == True: info(msg)

    ## Make genotype field and return output
    genotypeField = genotype + ':' + str(nbReadsALT) + ':' + str(totalNbReads)

    return (genotypeField)


def computeBkpVaf(bamFile, chrom, bkpPos, clippingPos, maxDist, minMapQual):
    '''
    '''

    tag = "VAF-" + clippingPos

    ### 1.Count the number of reads supporting the reference and alternative allele
    nbReadsREF = 0
    nbReadsALT = 0

    # Define genomic region +- X bp (3 default) around insertion breakpoint
    # + 1 for the end interval position, reason:
    # Convert from 1-based in VCF to 0-based half-open intervals in BAM [X,Y)
    # The position X is part of the interval, but Y is not.
    beg = bkpPos - maxDist
    end = bkpPos + (maxDist + 1)

    # Extract alignments overlapping the genomic region
    iterator = bamFile.fetch(chrom, beg, end)

    # Iterate over the alignments
    for alignment in iterator:

        ### Discard unmapped reads and reads with mapping quality lower than X (20 default) 
        if (alignment.is_unmapped == False) and (alignment.mapping_quality >= minMapQual):

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
                nbReadsREF += 1
                msg = "Alignment supports REF"
                if (debugBool == True): log(tag, msg)

            # B) Read do not supporting the reference allele
            else:
                # Assess if the alignment supports the alternative allele. Conditions:
                # a) bkpA: Only beginning of the read soft (Operation=4) or hard clipped (Operation=5)
                #    bkpB: Only end of the read soft (Operation=4) or hard clipped (Operation=5)
                # b) Alignment start position within +- 3bp compared to the insertion breakpoint
                # c) Not PCR nor optical duplicated (bitwise flag 0x400)

                firstOperation = alignment.cigartuples[0][0]
                lastOperation = alignment.cigartuples[-1][0]

                msg = "MEI allele info (firstOperation, lastOperation, duplicate, mapping_quality): " + str(firstOperation) + " " + str(lastOperation) + " " + str(duplicate) + " " + str(alignment.mapping_quality) 
            
                if (debugBool == True): log(tag, msg)

                ### Assess clipping

                ## Define genomic region +- X bp (3 default) around insertion breakpoint to search for clipped reads
                beg = bkpPos - maxDist
                end = bkpPos + maxDist

                # A) Clipping at the beginning of the read while not at the end
                # *******--------- (clipped bases: *)
                
                if (clippingPos == "Beg") and ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation != 4) and (lastOperation != 5)):

                    # a) clipping breakpoint within range 
                    # beg <---------> end
                    if (beg <= alignment.reference_start) and (alignment.reference_start <= end):
                        bkpBool = True

                    # b) bkp not within range
                    else:
                        bkpBool = False

                # B) Clipping at the end of the read while not at the beginning
                # ---------******* (clipped bases: *)
                elif (clippingPos == "End") and ((lastOperation == 4) or (lastOperation == 5)) and ((firstOperation != 4) and (firstOperation != 5)):

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

                ## Read supporting the alternative allele:
                if (bkpBool == True) and (duplicate == False):
                    nbReadsALT += 1

                    msg = "Alignment supports ALT"
                    if (debugBool == True): log(tag, msg)

        ### Unmapped or bad quality read filtered out
        else:
            msg = "Read filtered out (unmapped, mapping_quality): " + str(alignment.is_unmapped) + " " + str(alignment.mapping_quality) 
            if (debugBool == True): print msg

        msg = "-------------------"
        if (debugBool == True): print msg

    ### 2. Compute the VAF
    totalNbReads = nbReadsALT + nbReadsREF

    if totalNbReads == 0:
        vaf = 0
    else:
        vaf =  float(nbReadsALT) / totalNbReads

    msg = "VAF (nbReadsALT, totalNbReads, VAF): " + str(nbReadsALT) + " " + str(totalNbReads) + " " + str(vaf)
    if (debugBool == True): log(tag, msg)

    return (nbReadsALT, totalNbReads, vaf)


def computeAverageVaf(nbReadsALT_A, totalNbReads_A, vaf_A, nbReadsALT_B, totalNbReads_B, vaf_B, minClipped):
    '''
    '''

    ## a) Both breakpoints supported by at least X clipped-reads (default: 5)
    if (nbReadsALT_A >= minClipped) and (nbReadsALT_B >= minClipped):
        nbReadsALT = int(nbReadsALT_A + nbReadsALT_B) / 2
        totalNbReads = int(totalNbReads_A + totalNbReads_B) / 2
        vaf = float(nbReadsALT) / totalNbReads

    ## b) Breakpoint A supported by at least X clipped-reads (default: 5)
    elif (nbReadsALT_A >= minClipped):
        nbReadsALT = nbReadsALT_A
        totalNbReads = totalNbReads_A
        vaf = vaf_A

    ## c) Breakpoint B supported by at least X clipped-reads (default: 5)
    elif (nbReadsALT_B >= minClipped):
        nbReadsALT = nbReadsALT_B
        totalNbReads = totalNbReads_B
        vaf = vaf_B

    ## d) No breakpoint supported by at least X clipped-reads (default: 5)
    else:
        nbReadsALT = int(nbReadsALT_A + nbReadsALT_B) / 2
        totalNbReads =  int(totalNbReads_A + totalNbReads_B) / 2
        vaf = 0

    return nbReadsALT, totalNbReads, vaf


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
parser = argparse.ArgumentParser(description= "Genotype a set of germline mobile element insertions (MEI) across a group of donors. Produce a multi-sample VCF as output")
parser.add_argument('VCF', help='VCF with the MEI polymorphism to genotype across the set of donors')
parser.add_argument('BAMPaths', help='Tab separated text file with two columns: 1) donor ids and 2) path to the corresponding donor normal BAM files. One donor per row.')
parser.add_argument('outFileName', help='Identifier to name the output file.')
parser.add_argument('--maxDist', default=3, dest='maxDist', type=int, help='Maximum distance to breakpoint when searching for clipped reads supporting the alternative allele. Default: 3' )
parser.add_argument('--minClipped', default=5, dest='minClipped', type=int, help='Minimum number of clipped reads supporting the alternative allele. Default: 5' )
parser.add_argument('--minREF', default=5, dest='minREF', type=int, help='Minimum number of reads supporting the reference allele. Default: 5' )
parser.add_argument('--minMapQual', default=0, dest='minMapQual', type=int, help='Minimum mapping quality for considering a clipped read as supporting the alternative allele. Default: 0 (disabled)' )
parser.add_argument('--hetVaf', default=0.10, dest='hetVaf', type=float, help='Min VAF threshold for heterozygous (0/1). Default: 0.1' )
parser.add_argument('--homVaf', default=0.9, dest='homVaf', type=float, help='Min VAF threshold for homozygous alternative (1/1). Default: 0.9' )
parser.add_argument('-t', '--threads', default=1, dest='threads', type=int, help='Number of threads. Default: 1' )
parser.add_argument('--debug', action='store_true', dest='debug', help='Enable debug mode. Display detailed log info about MEI genotyping')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
BAMPaths = args.BAMPaths
outFileName =  args.outFileName
maxDist = args.maxDist
minClipped = args.minClipped
minREF = args.minREF
minMapQual = args.minMapQual
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
print "outFileName: ", outFileName
print "max-distance-bkp: ", maxDist
print "min-nb-clipped: ", minClipped
print "min-nb-reference: ", minREF
print "min-mapping-quality: ", minMapQual
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
    thread = threading.Thread(target=worker, args=(VCFObj, donorIdBamPathList, maxDist, minClipped, minMapQual, minREF, hetVaf, homVaf), name=threadName)
    threads.append(thread)

    counter += 1

# Launch threads
[t.start() for t in threads]

# Wait till threads have finished genotyping
[t.join() for t in threads]

#### 3. Make output multisample VCF file
header("3. Produce multi-sample VCF file as ouput")

fileName = outFileName + '.vcf'
outFilePath = outDir + '/' + fileName

# 3.1 Write header
VCFObj.write_header(outFilePath)

# 3.2 Write variants
donorIdList = [line[0] for line in donorIdBamPathList]

VCFObj.write_variants_multiSample(donorIdList, outFilePath)

header("Finished")        
