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
parser.add_argument('VCFfile', help='VCF containing structural variants')
parser.add_argument('genome', help='Reference genome in fasta format')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
VCFfile = args.VCFfile
genome = args.genome
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCFfile: ", VCFfile
print "genome: ", genome
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print



## Start ##

### 1.
#################################
   

info("Reading " + VCFfile + "...")

# Input VCF available
if os.path.isfile(VCFfile):

    # Create VCF object
    VCFObj = formats.VCF()

    # Read VCF and add information to VCF object
    VCFObj.read_VCF(VCFfile)

else:
    print "[ERROR] Input file does not exist"


### 3. 
#########################
bkpSeqPath = outDir + '/bkpSeqs.fa'

# For each MEI
for bkpObj in VCFObj.lineList:  

    if "CIPOS" in bkpObj.infoDict:
        CIPOS = bkpObj.infoDict["CIPOS"].split(",")
        intervalLen = abs(int(CIPOS[0])) + abs(int(CIPOS[1]))
        
    else:
        CIPOS = "0,0"
        intervalLen = 0

    #if (bkpObj.filter == "PASS") and (intervalLen <= 10):
    if (bkpObj.filter == "PASS") and (intervalLen <= 0):
        beg = int(bkpObj.pos) - 9
        end = int(bkpObj.pos) + 10                   
        targetInterval = bkpObj.chrom + ':' + str(beg) + '-' + str(end)
                
        # Extract sequence and print into fasta                    
        command="samtools faidx " + genome + ' ' + targetInterval + ' >> ' + bkpSeqPath
        os.system(command)

## End ##
print
print "***** Finished! *****"
print
