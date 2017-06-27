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

####### CLASSES #######
class cohort():
    """
        .....................

        Methods:
        -
    """
    def __init__(self):
        """
        """
        self.VCFdict = {}


    def read_VCFs(self, inputPath):
        """
        """
        inputFile = open(inputPath, 'r')

        info("Read input VCFs ")

        # Per iteration, read a VCF, generate a VCF object and add it to the cohort
        for line in inputFile:
            line = line.rstrip('\n')
            line = line.split("\t")

            donorId = line[0]
            projectCode = line[1].split("-")[0]
            VCFfile = line[2]

            #print "tiooo: ", donorId, projectCode, VCFfile

            # Create VCF object
            VCFObj = formats.VCF()

            info("Reading " + VCFfile + "...")

            # Input VCF available
            if os.path.isfile(VCFfile):

                # Read VCF and add information to VCF object
                VCFObj.read_VCF(VCFfile)

                # Initialize the donor list for a given project if needed 
                if projectCode not in self.VCFdict:

                    self.VCFdict[projectCode] = []   

                # Add donor VCF to cohort
                self.VCFdict[projectCode].append(VCFObj)

            else:
                print "[ERROR] Input file does not exist"

class fasta():
    """
    .... class.
    .....
    Methods:
    - fasta_reader
    """

    def __init__(self, fastaFile):
        """
            Initialize fasta object.
            Input:
            1)
            Output:
            -
        """
        self.fastaDict = self.fasta_reader(fastaFile)

    #### FUNCTIONS ####
    def write_fasta(self, outFilePath):
        """
            Input:
            1)
            Output:
            1)
        """
        fastaDict = {}
        outFile = open(outFilePath, 'w')
        
        for seqId in self.fastaDict:
            header = ">" + seqId + '\n'
            outFile.write(header)
            
            seq = self.fastaDict[seqId] + '\n'
            outFile.write(seq)

    def fasta_reader(self, fastaFile):
        """
            Input:
            1)
            Output:
            1)
        """
        fastaDict = {}

        subHeader("Fasta reader")

        fh = open(fastaFile)
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()

            # drop the info
            header = header.split(" ")[0]

            info("Reading " + header + "...")
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            fastaDict[header] = seq

        return fastaDict

    #### FUNCTIONS ####
    def rev_complement(self):
        """
            Make the reverse complementary of a fasta file
        """
           
        revFastaDict = {}

        for seqId in self.fastaDict:
    
            seq = self.fastaDict[seqId]
            baseComplementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
            seq = seq.upper()
            revSeq = seq[::-1] # Make reverse sequence
            letters = list(revSeq)
            letters = [baseComplementDict[base] for base in letters]
            revComplementSeq = ''.join(letters) # Make complement of reverse sequence

            revFastaDict[seqId] = revComplementSeq
        return revFastaDict

    def complement(self):
        """
            Make the complementary of a fasta file
        """
           
        complFastaDict = {}

        for seqId in self.fastaDict:
    
            seq = self.fastaDict[seqId]
            baseComplementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
            seq = seq.upper()
            letters = list(seq)
            letters = [baseComplementDict[base] for base in letters]
            complComplementSeq = ''.join(letters) # Make complement of the sequence

            complFastaDict[seqId] = complComplementSeq
        return complFastaDict

    def rev(self):
        """
            Make the reverse of a fasta file
        """
           
        revFastaDict = {}

        for seqId in self.fastaDict:
    
            seq = self.fastaDict[seqId]
            seq = seq.upper()
            revSeq = seq[::-1] # Make reverse sequence

            revFastaDict[seqId] = revSeq
        return revFastaDict


#### MAIN ####

## Import modules ##
import argparse
import sys
import os
import os.path
import formats
import time
from operator import itemgetter, attrgetter, methodcaller
from itertools import groupby


## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per donor with the following consecutive fields: tumorType donorId vcf_path')
parser.add_argument('genome', help='Reference genome in fasta format')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
genome = args.genome
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "genome: ", genome
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ##

### 1. Initialize cohort object
#################################
cohortObj = cohort()

### 2. Read VCF files, create VCF objects and organize them
#############################################################
cohortObj.read_VCFs(inputPath)

### 3. Extract sequences
#########################
out5PrimePlusPath = outDir + '/seqs_5primeBkp_plus.fa'
out5PrimeMinusPath = outDir + '/seqs_5primeBkp_minus.fa'
out3PrimePlusPath = outDir + '/seqs_3primeBkp_plus.fa'
out3PrimeMinusPath = outDir + '/seqs_3primeBkp_minus.fa'

## For each tumor type
for tumorType in cohortObj.VCFdict:

    ## For each donor 
    for VCFObj in cohortObj.VCFdict[tumorType]:

        # For each MEI
        for MEIObj in VCFObj.lineList:  
            
            ## Select only those MEI that pass all the filters 
            if (MEIObj.filter == "PASS"):


                ####### Insertion with at least the 5' bkp reconstructed (polyA extreme)
                if (int(MEIObj.infoDict["SCORE"]) == 3) or (int(MEIObj.infoDict["SCORE"]) == 5):

                    #### Select sequence around 5' insertion breakpoint
                    # A) Insertion in plus orientation 
                    #  (POS)>>>>>>>>>>>>>>AAAAAAAA-------
                    if (MEIObj.infoDict["STRAND"] == '+'):
                    
                        beg = int(MEIObj.pos) - 9
                        end = int(MEIObj.pos) + 10                   
                        targetInterval = MEIObj.chrom + ':' + str(beg) + '-' + str(end)
                
                        # Extract sequence and print into fasta   
                        command="samtools faidx " + genome + ' ' + targetInterval + ' >> ' + out5PrimePlusPath
                        os.system(command)

                    # B) Insertion in minus orientation 
                    # AAAAAAAA<<<<<<<<<<<<<(BKPB)-------
                    else:
                        beg = int(MEIObj.infoDict["BKPB"]) - 9
                        end = int(MEIObj.infoDict["BKPB"]) + 10                   
                        targetInterval = MEIObj.chrom + ':' + str(beg) + '-' + str(end)
            
                        # Extract sequence and print into fasta                    
                        command="samtools faidx " + genome + ' ' + targetInterval + ' >> ' + out5PrimeMinusPath
                        os.system(command)


                ####### Insertion with at least the 3' bkp reconstructed (polyA extreme)
                if (int(MEIObj.infoDict["SCORE"])>3):

                    #### Select sequence around 3' insertion breakpoint
                    # A) Insertion in plus orientation 
                    #  >>>>>>>>>>>>>>AAAAAAAA(BKPB)-------
                    if (MEIObj.infoDict["STRAND"] == '+'):
                    
                        beg = int(MEIObj.infoDict["BKPB"]) - 9
                        end = int(MEIObj.infoDict["BKPB"]) + 10                   
                        targetInterval = MEIObj.chrom + ':' + str(beg) + '-' + str(end)
                
                        # Extract sequence and print into fasta   
                        command="samtools faidx " + genome + ' ' + targetInterval + ' >> ' + out3PrimePlusPath
                        os.system(command)
    
                    # B) Insertion in minus orientation 
                    # (POS)AAAAAAAA<<<<<<<<<<<<<-------
                    else:
                        beg = int(MEIObj.pos) - 9
                        end = int(MEIObj.pos) + 10                   
                        targetInterval = MEIObj.chrom + ':' + str(beg) + '-' + str(end)
                
                        # Extract sequence and print into fasta                    
                        command="samtools faidx " + genome + ' ' + targetInterval + ' >> ' + out3PrimeMinusPath
                        os.system(command)
    

## Reformat 3' breakpoint sequences
######################################
## Reverse complement plus fasta:
outPlusRevPath = outDir + '/seqs_3primeBkp_plus_compl.fa'
fastaObj = fasta(out3PrimePlusPath)
fastaObj.fastaDict = fastaObj.complement()
fastaObj.write_fasta(outPlusRevPath)

## Reverse minus fasta:
outMinusRevPath = outDir + '/seqs_3primeBkp_minus_rev.fa'
fastaObj = fasta(out3PrimeMinusPath)
fastaObj.fastaDict = fastaObj.rev()
fastaObj.write_fasta(outMinusRevPath)



## End ##
print
print "***** Finished! *****"
print
