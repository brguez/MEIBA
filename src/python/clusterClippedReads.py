#!/usr/bin/env python
#coding: utf-8


# Class cluster:
#    clippedReadDict = {}
# Keys(readId) -> {}
#        key(alignmentObj)
#        key(seq)

     
#### CLASSES ####
class fasta():
    """
    .... class.

    .....

    Methods:
    - fasta_reader

    """

    def __init__(self):
        """
            Initialize fasta object.

            Input:
            1)

            Output:
            -
        """
        self.fastaDict = {}

    #### FUNCTIONS ####
    def fasta_reader(self, fastaFile):
        """
        """
        fastaDict = {}

        subHeader("Fasta reader")

        fh = open(fastaFile)
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()

            # drop the info
            header = header.split(" ")[0]

            info("Reading " + header + "...")
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            fastaDict[header] = seq

        self.fastaDict = fastaDict

    def write_fasta(self, outFilePath):
        """
        """
        outFile = open(outFilePath, "w" )

        for header, seq in self.fastaDict.iteritems():
            header = ">" + header

            outFile.write("%s\n" % header)
            outFile.write("%s\n" % seq)

        # Close output fasta file
        outFile.close()

class cluster():
    """
    """
    def __init__(self, alignmentObj, clippedSide):
        """
        """
        self.chrom = alignmentObj.reference_name      
        self.clippedSide = clippedSide
        self.bkpPos = alignmentObj.reference_start if clippedSide == "beg" else alignmentObj.reference_end
        self.clippedReadDict = {}       
        self.consensusSeq = ""

    def addClippedRead(self, alignmentObj):
        """
        """
        mate = '/1' if alignmentObj.is_read1 else '/2'
        readId = alignmentObj.query_name + mate
            
        self.bkpPos = alignmentObj.reference_start if self.clippedSide == "beg" else alignmentObj.reference_end        

        self.clippedReadDict[readId] = {}
        self.clippedReadDict[readId]["alignmentObj"] = alignmentObj

    def nbReads(self):
        """
        """
        return len(self.clippedReadDict)

    def addReadSeqs(self, fastaObj):
        """
        """
        for readId in self.clippedReadDict.keys():

            alignmentObj = self.clippedReadDict[readId]["alignmentObj"]
            
            ## Make the reverse complementary of reads aligned on the reverse strand
            if (alignmentObj.is_reverse == True):

                readSeq = rev_complement(fastaObj.fastaDict[readId])

                #print "SEQ: ", fastaObj.fastaDict[readId]
                #print "REV-SEQ: ", readSeq
            else:
                readSeq = fastaObj.fastaDict[readId]
            
            self.clippedReadDict[readId]["seq"]= readSeq

    def makeConsensusSeq(self):
        """
        """
                
        readIdAlignmentLenTupleList = [(readId, self.clippedReadDict[readId]["alignmentObj"].reference_length) for readId in self.clippedReadDict.keys()]
        
        # A) Single clipped read composing the cluster
        if len(readIdAlignmentLenTupleList) == 1:
            readId = readIdAlignmentLenTupleList[0][0]
            self.consensusSeq = self.clippedReadDict[readId]["seq"] 
            #print "CONSENSUS-1: ", self.consensusSeq

        # B) Multiple clipped reads composing the cluster
        else:  
            
            readIdAlignmentLenTupleSortedList =  sorted(readIdAlignmentLenTupleList, key=lambda x: x[1])
            lastElementIndex = len(readIdAlignmentLenTupleSortedList) - 1 

            minAlignmentSeq = self.clippedReadDict[readIdAlignmentLenTupleSortedList[0][0]]["seq"]  
            maxAlignmentSeq = self.clippedReadDict[readIdAlignmentLenTupleSortedList[lastElementIndex][0]]["seq"] 

            minAlignmentObj = self.clippedReadDict[readIdAlignmentLenTupleSortedList[0][0]]["alignmentObj"]  
            maxAlignmentObj = self.clippedReadDict[readIdAlignmentLenTupleSortedList[lastElementIndex][0]]["alignmentObj"] 

            ## a) Cluster composed by reads clipped at the end      
            #  maxAlignment
            # <-------------###
            #    <----------#####
            #         <-----##########
            #            <--##############
            #                minAlignment
            if (self.clippedSide == "end"):
                
                seqA = maxAlignmentSeq[:maxAlignmentObj.reference_length - 1]
                seqB = minAlignmentSeq[minAlignmentObj.reference_length:]
                self.consensusSeq = seqA + seqB

                #print "maxAlignment: ", maxAlignmentObj.is_reverse, maxAlignmentObj.reference_length, maxAlignmentSeq, seqA
                #print "minAlignment: ", minAlignmentObj.is_reverse, minAlignmentObj.reference_length, minAlignmentSeq, seqB
                #print "consensusSeq: ", self.consensusSeq

            ## b) Cluster composed by reads clipped at the beg      
            #                 maxAlignment
            #             ###-------------->
            #           #####---------->
            #      ##########----->
            #  ##############-->
            #   minAlignment              
            else:
                readLenA = len(minAlignmentSeq)
                clippedLenA = readLenA - minAlignmentObj.reference_length 
                
                readLenB = len(maxAlignmentSeq)
                clippedLenB = readLenB - maxAlignmentObj.reference_length 

                seqA = minAlignmentSeq[:clippedLenA - 1]
                seqB = maxAlignmentSeq[clippedLenB:]               
                self.consensusSeq = seqA + seqB

                #print "TIOO: ", readLenA, clippedLenA, readLenB, clippedLenB
                #print "minAlignment: ", minAlignmentObj.is_reverse, minAlignmentObj.reference_length, minAlignmentSeq, seqA
                #print "maxAlignment: ", maxAlignmentObj.is_reverse, maxAlignmentObj.reference_length, maxAlignmentSeq, seqB
                #print "consensusSeq: ", self.consensusSeq

           
#### FUNCTIONS ####
def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string

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


def rev_complement(seq):
    """
        Make the reverse complementary of a dna sequence

        Input:
        1) seq. DNA sequence

        Output:
        1) revComplementSeq. Reverse complementary of input DNA sequence
    """
    baseComplementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    seq = seq.upper()
    revSeq = seq[::-1] # Make reverse sequence
    letters = list(revSeq)
    letters = [baseComplementDict[base] for base in letters]
    revComplementSeq = ''.join(letters) # Make complement of reverse sequence

    return revComplementSeq

def overlap(begA, endA, begB, endB):

    """
    Check if both ranges overlap. 2 criteria for defining overlap: 
    ## A) Begin of the range A within the range B         
    #       *beg* <---------range_A---------->                         
    # <---------range_B----------> 
                
    #    *beg* <-------range_A----->
    # <-------------range_B------------------>
    ## B) Begin of the range B within the range A     
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->
            
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>
    """    
       
    # a) Begin of the range A within the range B   
    if ((begA >= begB) and (begA <= endB)):
        overlap = True
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True

    # c) Ranges do not overlapping
    else:
        overlap = False

    return overlap


def getClipped(chrom, beg, end, bamFile):
    '''
    '''
    #print "** pickClipped function **"
    
    clippedBegList = []
    clippedEndList = []
    
    # Extract alignments in the interval
    iterator = bamFile.fetch(chrom, beg, end)

    # Iterate over the alignments
    for alignmentObj in iterator:

        ### Discard unmapped reads 
        if (alignmentObj.is_unmapped == False):

            firstOperation = alignmentObj.cigartuples[0][0]
            lastOperation = alignmentObj.cigartuples[-1][0]

            #### Cleck if soft-clipped read
            # Note: soft (Operation=4) or hard clipped (Operation=5)
            # Discard reads clipped both in the beginning and ending

            ## a) Clipping at the beginning of the read while not clipping at all at the end
            # *******--------- (clipped bases: *)     
            if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation != 4) and (lastOperation != 5)):

                #print "BEG-SOFT: ", alignmentObj, alignmentObj.query_name, alignmentObj.query_sequence, alignmentObj.reference_name, alignmentObj.reference_start       
                clippedBegList.append(alignmentObj)
                
            ## b) Clipping at the end of the read while not clipping at all at the beginning
            # ---------******* (clipped bases: *)     
            elif ((lastOperation == 4) or (lastOperation == 5)) and ((firstOperation != 4) and (firstOperation != 5)):
                #print "END-SOFT: ", alignmentObj, alignmentObj.query_name, alignmentObj.query_sequence, alignmentObj.reference_name, alignmentObj.reference_start

                clippedEndList.append(alignmentObj)      
  
    return clippedBegList, clippedEndList


def clusterCLipped(clippedList, clippedSide):
    '''
    '''
    #print "** clusterCLipped  function **"

    ### Sort the list of clipped reads in increasing coordinates order
    if (clippedSide == "beg"):    
        clippedSortedList = sorted(clippedList, key=lambda alignmentObj: alignmentObj.reference_start, reverse=False)
    else:
        clippedSortedList = sorted(clippedList, key=lambda alignmentObj: alignmentObj.reference_end, reverse=False)

    ### Make clipped read clusters:
    clusterList = []

    ## For each clipped read alignment 
    for alignmentObj in clippedSortedList:

        # A) No cluster in the list -> Create first cluster
        if not clusterList:

            #msg = "Initialize first cluster"
            #log("CLUSTER", msg) 

            clusterObj = cluster(alignmentObj, clippedSide) 
            clusterObj.addClippedRead(alignmentObj)
            clusterList.append(clusterObj)

            #print "TEST: ", clusterObj.chrom, clusterObj.clippedSide, clusterObj.bkpPos, clusterObj.clippedReadDict

        # B) There is already at least one cluster in the list -> Check if current clipped read within the latest cluster
        else:
                    
            #msg = "Check if clipped read within latest cluster"
            #log("CLUSTER", msg) 

            ## Define bkp position:
            bkpPos = alignmentObj.reference_start if clippedSide == "beg" else alignmentObj.reference_end        

            ## Define cluster range for searching for overlap
            lastClusterObj = clusterList[-1]     
            begClusterRange = lastClusterObj.bkpPos 
            endClusterRange = lastClusterObj.bkpPos + 3
            #endClusterRange = lastClusterObj.bkpPos + 1
    
            #### Check if clipped read within cluster range
            overlapping = overlap(bkpPos, bkpPos, begClusterRange, endClusterRange) 

            #msg = "cluster_range,clipped_range: " + str(begClusterRange) + " " + str(endClusterRange) + " " + str(bkpPos) + " " + str(bkpPos) + " " + str(overlapping)
            #log("CLUSTER", msg) 
    
            ## a) Overlapping ranges, so clipped read within previous cluster interval -> add read to the cluster                                
            if overlapping:

                #msg = "clipped read within cluster -> add to the cluster"
                #log("CLUSTER", msg) 
                lastClusterObj.addClippedRead(alignmentObj)
                     
            ## b) Clipped read outside previous cluster interval -> create new cluster and add it into the list
            else:
            
                #msg = "Clipped read outside the cluster -> create new cluster "
                #log("CLUSTER", msg) 
                clusterObj = cluster(alignmentObj, clippedSide) 
                clusterObj.addClippedRead(alignmentObj)
                clusterList.append(clusterObj)
                 
            #msg = "Number of clipped reads within cluster: ", len(clusterObj.clippedReadDict)
            #log("CLUSTER", msg) 
            #msg = "----------------------"
            #log("CLUSTER", msg) 

    return clusterList

        
def filterClusters(clusterList):
    '''
    '''
    print "** filterClusters  function **"

    for clusterObj in clusterList:
        print clusterObj, len(clusterObj.clippedReadDict), clusterObj.nbReads()



#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import time
from operator import itemgetter, attrgetter, methodcaller
import pysam
import itertools 
import subprocess

# Global variables:
global debugBool ## debug logging mode. Boolean.

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('insertions', help='')
parser.add_argument('bam', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
insertionsPath = args.insertions
bam = args.bam
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "insertionsPath: ", insertionsPath
print "bam: ", bam
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

## 1) 
##############################

## Open input files
insertions = open(insertionsPath, 'r')

## Open donor's BAM file for reading
bamFile = pysam.AlignmentFile(bam, "rb")

clustersDict = {}

## Read insertions file line by line
for line in insertions:

    ## Ignore comment lines (e.g. header)
    if line.startswith('#'):
        continue
        
    line = line.rstrip('\n')
    fieldsList = line.split("\t")

    ## Insertion line with the expected number of columns
    if (int(len(fieldsList)) == 31):
        chrPlus = fieldsList[0]
        begPlus = fieldsList[1]
        endPlus = fieldsList[2] 
        nbReadsPlus = fieldsList[3]
        familyPlus = fieldsList[4]
        readPairListPlus = fieldsList[5].split(",")
        chrMinus = fieldsList[6]
        begMinus = fieldsList[7]
        endMinus = fieldsList[8]
        nbReadsMinus = fieldsList[9]
        familyMinus = fieldsList[10]
        readPairListMinus = fieldsList[11].split(",")
        insertionType = fieldsList[12]
        rgType = fieldsList[30]

        print "###### INSERTION: ", chrPlus, begPlus, endPlus, chrMinus, begMinus, endMinus, rgType

        ## Define an insertion id (insertion coordinates defined by the end
        # of + cluster and beg of - cluster)
        insertionId = familyPlus + ":" + insertionType + ":" + chrPlus + "_" + endPlus + "_" + begMinus 


        ### 1. Search for clipped reads
        ## Positive cluster
        chrom = chrPlus

        if (rgType == "DUP"):
            beg = int(begPlus) - 200
            end = int(begPlus) + 200
        else:
            beg = int(endPlus) + 100 - 200
            end = int(endPlus) + 100 + 200

        print "range_+: ", chrom, beg, end
        clippedBegPlusList, clippedEndPlusList = getClipped(chrom, beg, end, bamFile)

        #print "clipped-Reads-Plus: ", len(clippedBegPlusList), len(clippedEndPlusList), clippedBegPlusList, clippedEndPlusList

        ## Negative cluster
        chrom = chrMinus

        if (rgType == "DUP"):
            beg = int(endMinus) + 100 - 200
            end = int(endMinus) + 100 + 200

        else:
            beg = int(begMinus) - 200
            end = int(begMinus) + 200
        
        print "range_-: ", chrom, beg, end
        clippedBegMinusList, clippedEndMinusList = getClipped(chrom, beg, end, bamFile)

        #print "clipped-Reads-Minus: ", len(clippedBegMinusList), len(clippedEndMinusList), clippedBegMinusList, clippedEndMinusList

        ### 2. Merge clipped read lists:
        clippedBegList = list(set(clippedBegPlusList + clippedBegMinusList))
        clippedEndList = list(set(clippedEndPlusList + clippedEndMinusList))

        ### 3. Cluster clipped reads:
        ## CLipping at the beginn
        clusterBegList = clusterCLipped(clippedBegList, "beg")

        print "clusterBegList: ", len(clusterBegList), clusterBegList

        ## Clipping at the end
        clusterEndList = clusterCLipped(clippedEndList, "end")
         
        print "clusterEndList: ", len(clusterEndList), clusterEndList

        ### 4. Filter clusters:
        #filterClusters(clusterBegList)
        
        ### 5. Add the 2 cluster lists to the dictionary:
        clustersDict[insertionId] = {}
        clustersDict[insertionId]["beg"] = clusterBegList
        clustersDict[insertionId]["end"] = clusterEndList

#print "clustersDict: ", clustersDict

bamFile.close()

## 2) Make fasta containing the reads supporting the clusters of clipped reads
################################################################################

## 1. Make file containing the list of read pair ids supporting the clipped clusters
allReadPairIdList = []

for insertionId in clustersDict:
   
    clusterBegList = clustersDict[insertionId]["beg"] 
    clusterEndList = clustersDict[insertionId]["end"]   

    for clusterObj in clusterBegList:
        readPairIdList = [readId.split("/")[0] for readId in list(clusterObj.clippedReadDict.keys())]
        #readIdList = [readId for readId in list(clusterObj.clippedReadDict.keys())]
        allReadPairIdList = allReadPairIdList + readPairIdList

    for clusterObj in clusterEndList:
        readPairIdList = [readId.split("/")[0] for readId in list(clusterObj.clippedReadDict.keys())]
        #readIdList = [readId for readId in list(clusterObj.clippedReadDict.keys())]
        allReadPairIdList = allReadPairIdList + readPairIdList    

allReadPairIdList = list(set(allReadPairIdList))


## Write read pair ids in an output file:
readPairsPath = outDir +'/allReadPairs.txt'
readPairsFile = open(readPairsPath, 'w')

for readPairId in allReadPairIdList:
    row = readPairId + "\n"
    readPairsFile.write(row)

## Important to close! otherwhise next step won't work properly...
readPairsFile.close()

## 2. Extract clipped read sequences with picard and generate fasta containing all the reads supporting the clusters
# PICARD = "java -jar /Users/brodriguez/Research/Apps/Picard/2.12.1/picard.jar"
# PICARD="java -Xms10G -Xmx10G -jar /software/CGP/external-apps/picard-tools-1.80/lib/FilterSamReads.jar" # Sanger
PICARD = "java -jar /Users/brodriguez/Research/Apps/Picard/2.12.1/picard.jar FilterSamReads" # Laptop

readPairsFasta = outDir + '/allReadPairs.fa'

command = PICARD + ' I=' + bam + ' O=/dev/stdout READ_LIST_FILE=' + readPairsPath + ' FILTER=includeReadList WRITE_READS_FILES=false VALIDATION_STRINGENCY=SILENT QUIET=true | samtools fasta - > '  + readPairsFasta
print command
os.system(command) # returns the exit status

## 3) Add to the reads supporting the clusters its complete sequence from fasta and 
####################################################################################
# generate consensus sequence
##############################
fastaObj = fasta()
fastaObj.fasta_reader(readPairsFasta)

for insertionId in clustersDict:
    #print "********** ", insertionId, " *************"
    clusterBegList = clustersDict[insertionId]["beg"] 
    clusterEndList = clustersDict[insertionId]["end"]   

    #print "--- clusterBeg ---"
    for clusterObj in clusterBegList:
        clusterObj.addReadSeqs(fastaObj) 
        #print "TIO: ", clusterObj.clippedReadDict
        clusterObj.makeConsensusSeq()
        #print "......"

    #print "--- clusterEnd ---"
    for clusterObj in clusterEndList:
        clusterObj.addReadSeqs(fastaObj) 
        #print "TIO: ", clusterObj.clippedReadDict
        clusterObj.makeConsensusSeq()
        #print "......"

## 4) For each insertion generate a fasta containing the consensus sequences for each cluster
##############################################################################################
for insertionId in clustersDict:
    print "********** ", insertionId, " *************"
    
    fastaDict = {}
    clusterList = clustersDict[insertionId]["beg"] + clustersDict[insertionId]["end"]

    for clusterObj in clusterList:
        header = "cluster" + "_" + clusterObj.chrom + "_" + str(clusterObj.bkpPos) + "_" + clusterObj.clippedSide + "_" + str(clusterObj.nbReads())
        seq = clusterObj.consensusSeq
        fastaDict[header] = seq

    fastaObj = fasta()
    fastaObj.fastaDict = fastaDict
    print insertionId, fastaObj.fastaDict

    fileName = insertionId + ".fa"
    outFilePath = outDir + "/" + fileName
    fastaObj.write_fasta(outFilePath)


print "***** Finished! *****"
print

