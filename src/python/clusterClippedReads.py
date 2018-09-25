#!/usr/bin/env python
#coding: utf-8


#### CLASSES ####
class fasta():
    """
    """

    def __init__(self):
        """
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

        operation = alignmentObj.cigartuples[0][0] if self.clippedSide == "beg" else alignmentObj.cigartuples[-1][0]
        clipType = "soft" if operation == 4 else "hard"

        self.clippedReadDict[readId] = {}
        self.clippedReadDict[readId]["alignmentObj"] = alignmentObj
        self.clippedReadDict[readId]["clipType"] = clipType

    def nbReads(self):
        """
        """
        return len(self.clippedReadDict)

    def readIdList(self):
        """
        """
        return list(self.clippedReadDict.keys())

    def addReadSeqs(self, fastaObj):
        """
        """
        for readId in self.clippedReadDict.keys():

            alignmentObj = self.clippedReadDict[readId]["alignmentObj"]
            
            ## Make the reverse complementary of reads aligned on the reverse strand
            if (alignmentObj.is_reverse == True):

                readSeq = rev_complement(fastaObj.fastaDict[readId])

            else:
                readSeq = fastaObj.fastaDict[readId]
            
            self.clippedReadDict[readId]["seq"]= readSeq

    def makeConsensusSeq(self, outDir):
        """
        multiple sequence alignment based
        """


        ## A) Single sequence
        if len(self.clippedReadDict.keys()) == 1:    
            consensusSeq = list(self.clippedReadDict.values())[0]["seq"].upper()

        ## B) Multiple sequence
        else:
            command = 'mkdir -p ' + outDir 
            os.system(command) # returns the exit status
            
            ### 1. Create fasta file containing cluster supporting reads
            fastaObj = fasta()        
            fastaDict = {}
    
            for readId in self.clippedReadDict.keys():
                fastaDict[readId] = self.clippedReadDict[readId]["seq"]

            fastaObj.fastaDict = fastaDict
            fastaPath = outDir + '/supportingReads.fa'
            fastaObj.write_fasta(fastaPath)
        
            ### 2. Make multiple sequence alignment
            msfPath = outDir + '/supportingReads.msf'
            command = 'muscle -in ' + fastaPath + ' -out ' + msfPath + ' -msf' 
            print command
            os.system(command) # returns the exit status

            ### 3. Generate consensus sequence (cons tool from EMBOSS packagge)
            consensusPath = outDir + '/consensus.fa'
    
            command = 'cons -sequence ' + msfPath + ' -outseq ' + consensusPath + ' -identity 0 -plurality 0'
            print command
            os.system(command) # returns the exit status

            ### Read consensus sequence 
            fastaObj = fasta() 
            fastaObj.fasta_reader(consensusPath)
            consensusSeq = fastaObj.fastaDict["EMBOSS_001"].upper()

            ### Do cleanup
            command = 'rm ' + fastaPath + ' ' + msfPath + ' ' + consensusPath             
            os.system(command) # returns the exit status

        ## Replace '-' by 'N' for ambiguous bases:
        consensusSeq = consensusSeq.replace('-', 'N')

        ## Convert consensus sequence into upper case:
        consensusSeq = consensusSeq.upper()

        return consensusSeq
    
              
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


def getClippedPairedClusters(chrPlus, begPlus, endPlus, chrMinus, begMinus, endMinus, rgType, bamFile, windowSize):
    """
    """
    ## 1. Extract clipped reads for positive cluster
    chrom = chrPlus

    if (rgType == "DUP"):
        beg = int(begPlus) - windowSize
        end = int(begPlus) + windowSize
    else:
        beg = int(endPlus) + 100 - windowSize
        end = int(endPlus) + 100 + windowSize

    clippedBegPlusList, clippedEndPlusList = getClippedInterval(chrom, beg, end, bamFile)

    ## 2. Extract clipped reads for negative cluster
    chrom = chrMinus

    if (rgType == "DUP"):
        beg = int(endMinus) + 100 - windowSize
        end = int(endMinus) + 100 + windowSize

    else:
        beg = int(begMinus) - windowSize
        end = int(begMinus) + windowSize
        
    print "range_-: ", chrom, beg, end
    clippedBegMinusList, clippedEndMinusList = getClippedInterval(chrom, beg, end, bamFile)

    ## 3. Merge clipped read lists:
    clippedBegList = list(set(clippedBegPlusList + clippedBegMinusList))
    clippedEndList = list(set(clippedEndPlusList + clippedEndMinusList))

    return clippedBegList, clippedEndList


def getClippedUnpairedCluster(chrPlus, begPlus, endPlus, bamFile, windowSize):
    """
    """
    ## 1. Extract clipped reads for cluster beginning
    chrom = chrPlus
    beg = int(begPlus) - windowSize
    end = int(begPlus) + windowSize
    
    print "range_beg: ", chrom, beg, end
    clippedBegClusterBegList, clippedEndClusterBegList = getClippedInterval(chrom, beg, end, bamFile)

    ## 2. Extract clipped reads for cluster ending
    chrom = chrPlus
    beg = int(endPlus) + 100 - windowSize
    end = int(endPlus) + 100 + windowSize

    print "range_end: ", chrom, beg, end
    clippedBegClusterEndList, clippedEndClusterEndList = getClippedInterval(chrom, beg, end, bamFile)

    ## 3. Merge clipped read lists:
    clippedBegList = list(set(clippedBegClusterBegList + clippedBegClusterEndList))
    clippedEndList = list(set(clippedEndClusterBegList + clippedEndClusterEndList))

    return clippedBegList, clippedEndList


def getClippedInterval(chrom, beg, end, bamFile):
    '''
    '''
    #print "** pickClipped function **"
    
    clippedBegList = []
    clippedEndList = []
    
    # Extract alignments in the interval
    iterator = bamFile.fetch(chrom, beg, end)

    # Iterate over the alignments
    for alignmentObj in iterator:

        ### Discard unmapped reads and PCR duplicates
        if (alignmentObj.is_unmapped == False) and (alignmentObj.is_duplicate == False):

            firstOperation = alignmentObj.cigartuples[0][0]
            lastOperation = alignmentObj.cigartuples[-1][0]

            #### Cleck if soft-clipped read
            # Note: soft (Operation=4) or hard clipped (Operation=5)
            # Discard reads clipped both in the beginning and ending
            ## a) Clipping at the beginning of the read while not clipping at all at the end
            # *******--------- (clipped bases: *)     
            if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation != 4) and (lastOperation != 5)):
                clippedBegList.append(alignmentObj)
                
            ## b) Clipping at the end of the read while not clipping at all at the beginning
            # ---------******* (clipped bases: *)     
            elif ((lastOperation == 4) or (lastOperation == 5)) and ((firstOperation != 4) and (firstOperation != 5)):
                clippedEndList.append(alignmentObj)      
  
    return clippedBegList, clippedEndList



def clusterCLipped(clippedList, clippedSide, minNbReads, maxNbReads):
    '''
    '''
    #print "** clusterCLipped  function **"

    ### 1. Sort the list of clipped reads in increasing coordinates order
    if (clippedSide == "beg"):    
        clippedSortedList = sorted(clippedList, key=lambda alignmentObj: alignmentObj.reference_start, reverse=False)
    else:
        clippedSortedList = sorted(clippedList, key=lambda alignmentObj: alignmentObj.reference_end, reverse=False)

    ### 2. Make clipped read clusters:
    clusterList = []

    ## For each clipped read alignment 
    for alignmentObj in clippedSortedList:

        # A) No cluster in the list -> Create first cluster
        if not clusterList:

            clusterObj = cluster(alignmentObj, clippedSide) 
            clusterObj.addClippedRead(alignmentObj)
            clusterList.append(clusterObj)

        # B) There is already at least one cluster in the list -> Check if current clipped read within the latest cluster
        else:

            ## Define bkp position:
            bkpPos = alignmentObj.reference_start if clippedSide == "beg" else alignmentObj.reference_end        

            ## Define cluster range for searching for overlap
            lastClusterObj = clusterList[-1]     
            begClusterRange = lastClusterObj.bkpPos 
            endClusterRange = lastClusterObj.bkpPos + 3
    
            #### Check if clipped read within cluster range
            overlapping = overlap(bkpPos, bkpPos, begClusterRange, endClusterRange) 
    
            ## a) Overlapping ranges, so clipped read within previous cluster interval -> add read to the cluster                                
            if overlapping:

                lastClusterObj.addClippedRead(alignmentObj)
                     
            ## b) Clipped read outside previous cluster interval -> create new cluster and add it into the list
            else:

                clusterObj = cluster(alignmentObj, clippedSide) 
                clusterObj.addClippedRead(alignmentObj)
                clusterList.append(clusterObj)


    ### 3. Filter the clusters according to the number of reads supporting them (min and max cut-offs)
    filteredClusterList = []

    for clusterObj in clusterList: 

        if (clusterObj.nbReads() >= minNbReads) and (clusterObj.nbReads() <= maxNbReads):
            filteredClusterList.append(clusterObj)            

    return filteredClusterList


def filterNbClusters(clusterBegList, clusterEndList, maxNbClusters):
    '''
    '''
    totalNbClusters = len(clusterBegList) + len(clusterEndList)
    
    ## A) Number of clipped clusters higher than the treshold -> Discard clusters as most likely are the consequence of
    # alignment artefacts. In a perfect scenario we would expect two clusters, a single one per breakpoint
    if (totalNbClusters > maxNbClusters):
        filteredClusterBegList = []
        filteredClusterEndList = []

    ## B) Pass the filter            
    else:
        filteredClusterBegList = clusterBegList
        filteredClusterEndList = clusterEndList

    return filteredClusterBegList, filteredClusterEndList


def filterDiscordantCluster(chrom, beg, end, readPairList, bamFile):
    '''
    '''
    nbDiscordant = len(readPairList)
    nbClippedBothSides = 0
    readPairFilteredList = []

    ## Extract alignments in the interval
    iterator = bamFile.fetch(chrom, beg, end)

    ## Iterate over the alignments
    for alignmentObj in iterator:

        ## Supporting discordant paired-end read and cigar available
        if (alignmentObj.query_name in readPairList) and (alignmentObj.cigartuples is not None):

            firstOperation = alignmentObj.cigartuples[0][0]
            lastOperation = alignmentObj.cigartuples[-1][0]

            ### A) Read clipped both in the beginning and ending
            if ((firstOperation == 4) or (firstOperation == 5)) and ((lastOperation == 4) or (lastOperation == 5)):
                nbClippedBothSides += 1
    
            ### B) Read not clipped in both sides
            else:
                readPairFilteredList.append(alignmentObj.query_name)

    ## Percentage of supporting paired ends that are clipped on both sides
    percClippedBothSides = float(nbClippedBothSides) / nbDiscordant * 100

    ## Recompute the number of supporting paired ends after removing problematic reads
    readPairFilteredList = list(set(readPairFilteredList))
    nbFilteredDiscordant = len(readPairFilteredList)

    ## Discard cluster if more than 50% supporting paired-ends clipped on both sides:
    if (percClippedBothSides > 50):
        print "FILTER-CLUSTER: ", nbClippedBothSides, nbDiscordant, percClippedBothSides, nbFilteredDiscordant, readPairFilteredList
        readPairFilteredList = []
        nbFilteredDiscordant = 0
        filtered = True
    else:
        filtered = False

    return filtered


#### MAIN ####

## Import modules ##
import argparse
import sys
import os
import time
from operator import itemgetter, attrgetter, methodcaller
import pysam
import itertools 
import subprocess

# Global variables:
global debugBool ## debug logging mode. Boolean.

# Environmental variables:
PICARD = os.environ['PICARD']

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('insertions', help='')
parser.add_argument('bam', help='Bam file')
parser.add_argument('--windowSize', default=50, dest='windowSize', type=int, help='Window size to search for clipped read clusters from discordant read-pair clusters ends. Default=50bp' )
parser.add_argument('--minNbReads', default=1, dest='minNbReads', type=int, help='Minimum number of clipped reads composing the cluster. Default: 1' )
parser.add_argument('--maxNbReads', default=500, dest='maxNbReads', type=int, help='Maximum number of clipped reads composing the cluster. Default: 500' )
parser.add_argument('--maxNbClusters', default=10, dest='maxNbClusters', type=int, help='Maximum number of clipped read clusters in the insertion region. Default: 10' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
insertionsPath = args.insertions
bam = args.bam
windowSize = args.windowSize
minNbReads = args.minNbReads
maxNbReads = args.maxNbReads
maxNbClusters = args.maxNbClusters
outDir = args.outDir
tmpDir = outDir + '/tmp'

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "insertionsPath: ", insertionsPath
print "bam: ", bam
print "windowSize: ", windowSize
print "minNbReads: ", minNbReads
print "maxNbReads: ", maxNbReads
print "maxNbClusters: ", maxNbClusters
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

## Open input files
insertions = open(insertionsPath, 'r')

## Open donor's BAM files for reading
bamFile = pysam.AlignmentFile(bam, "rb")

clustersDict = {}
discordantReadPairList = []

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

        ## Add discordant read pairs to the list:
        discordantReadPairList = discordantReadPairList + readPairListPlus + readPairListMinus        

        ## Define an insertion id (insertion coordinates defined by the end
        # of + cluster and beg of - cluster)
        if familyPlus == 'Other': # temporary fix
            familyPlus = 'SVA'

        insertionId = familyPlus + ":" + insertionType + ":" + chrPlus + "_" + endPlus + "_" + begMinus 

        ### 0. Refine discordant paired end clusters:
        ## A) Paired clusters
        if (begMinus != "NA") and (begMinus != "UNK"):
            filteredPlus = filterDiscordantCluster(chrPlus, int(begPlus), int(endPlus) + 100, readPairListPlus, bamFile)
            filteredMinus = filterDiscordantCluster(chrMinus, int(begMinus), int(endMinus) + 100, readPairListMinus, bamFile)

        ## B) Unpaired cluster
        else:
            filteredPlus = filterDiscordantCluster(chrPlus, int(begPlus), int(endPlus) + 100, readPairListPlus, bamFile)
            filteredMinus = False

        ## Discard those insertions with a high percentage of both-sides clipped reads supporting at least one of the clusters:
        if (filteredPlus == True) or (filteredMinus == True):
            clusterBegFilteredList = []
            clusterEndFilteredList = []

        else:
            ### 1. Search for clipped reads
            ## A) Paired clusters
            if (begMinus != "NA") and (begMinus != "UNK"):
                clippedBegList, clippedEndList = getClippedPairedClusters(chrPlus, begPlus, endPlus, chrMinus, begMinus, endMinus, rgType, bamFile, windowSize)
       
            ## B) Unpaired cluster
            else:
                clippedBegList, clippedEndList = getClippedUnpairedCluster(chrPlus, begPlus, endPlus, bamFile, windowSize)

            ### 2. Cluster clipped reads:
            ### 2.1 Tumour
            clusterBegList = clusterCLipped(clippedBegList, "beg", minNbReads, maxNbReads)
            clusterEndList = clusterCLipped(clippedEndList, "end", minNbReads, maxNbReads)       

            ### 3. Filter clusters of clipped reads:
            ## 3.1 Filter by the number of clipped-read clusters              
            clusterBegFilteredList, clusterEndFilteredList = filterNbClusters(clusterBegList, clusterEndList, maxNbClusters)

        ### 4. Add the 2 cluster lists to the dictionary:
        clustersDict[insertionId] = {}
        clustersDict[insertionId]["beg"] = clusterBegFilteredList
        clustersDict[insertionId]["end"] = clusterEndFilteredList

bamFile.close()

## 2) Make fasta containing the discordant paired-end reads + 
##############################################################
# the reads supporting the clusters of clipped reads 
####################################################

## 1. Make list containing the discordant paired-end reads
allReadPairIdList = discordantReadPairList

## 2. Add to the list the reads supporting the clusters of clipped reads 
for insertionId in clustersDict:
   
    clusterBegList = clustersDict[insertionId]["beg"] 
    clusterEndList = clustersDict[insertionId]["end"]   

    for clusterObj in clusterBegList:
        readPairIdList = [readId.split("/")[0] for readId in clusterObj.readIdList()]
        allReadPairIdList = allReadPairIdList + readPairIdList

    for clusterObj in clusterEndList:
        readPairIdList = [readId.split("/")[0] for readId in clusterObj.readIdList()]
        allReadPairIdList = allReadPairIdList + readPairIdList    

allReadPairIdList = list(set(allReadPairIdList))

## 3. Make file containing the supporting read ids
readPairsPath = outDir +'/allReadPairs.txt'
readPairsFile = open(readPairsPath, 'w')

for readPairId in allReadPairIdList:
    row = readPairId + "\n"
    readPairsFile.write(row)

## Important to close! otherwhise next step won't work properly...
readPairsFile.close()

## 4. Extract read sequences with picard and generate fasta
readPairsFasta = outDir + '/allReadPairs.fa'

command = PICARD + ' FilterSamReads I=' + bam + ' O=/dev/stdout READ_LIST_FILE=' + readPairsPath + ' FILTER=includeReadList WRITE_READS_FILES=false VALIDATION_STRINGENCY=SILENT QUIET=true | samtools fasta - > '  + readPairsFasta
print command
os.system(command)


## 3) Add to the reads supporting the clusters its complete sequence from fasta and 
####################################################################################
# generate consensus sequence
##############################
fastaObj = fasta()
fastaObj.fasta_reader(readPairsFasta)

for insertionId in clustersDict:
    print "********** ", insertionId, " *************"
    clusterBegList = clustersDict[insertionId]["beg"] 
    clusterEndList = clustersDict[insertionId]["end"]   

    print "--- clusterBeg ---"
    for clusterObj in clusterBegList:
        clusterId = clusterObj.chrom + "_" + str(clusterObj.bkpPos) + "_" + clusterObj.clippedSide + "_" + str(clusterObj.nbReads()) 
        consensusDir = tmpDir + '/' + clusterId
        clusterObj.addReadSeqs(fastaObj) 
        clusterObj.consensusSeq = clusterObj.makeConsensusSeq(consensusDir)

    #print "--- clusterEnd ---"
    for clusterObj in clusterEndList:
        clusterId = clusterObj.chrom + "_" + str(clusterObj.bkpPos) + "_" + clusterObj.clippedSide + "_" + str(clusterObj.nbReads())          
        consensusDir = tmpDir + '/' + clusterId
        clusterObj.addReadSeqs(fastaObj) 
        clusterObj.consensusSeq = clusterObj.makeConsensusSeq(consensusDir)


## 4) For each insertion generate a fasta containing the consensus sequences for each cluster
##############################################################################################
for insertionId in clustersDict:
    print "********** ", insertionId, " *************"
    
    fastaDict = {}
    clusterList = clustersDict[insertionId]["beg"] + clustersDict[insertionId]["end"]

    ## For each cluster
    for clusterObj in clusterList:

        ## Include into the header the clipped read ids..
        header = "cluster" + "_" + clusterObj.chrom + "_" + str(clusterObj.bkpPos) + "_" + clusterObj.clippedSide + "_" + str(clusterObj.nbReads()) + "\t" + ",".join(clusterObj.readIdList())
        fastaDict[header] = clusterObj.consensusSeq

    fastaObj = fasta()
    fastaObj.fastaDict = fastaDict

    ## Write into the output file
    fileName = insertionId + ".fa"
    outFilePath = outDir + "/" + fileName
    fastaObj.write_fasta(outFilePath)


### Make cleanup and finish
command = 'rm -r ' + readPairsPath + ' ' + tmpDir 
os.system(command) # returns the exit status


print "***** Finished! *****"
print

