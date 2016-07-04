#!/usr/bin/env python
#coding: utf-8 

## Load modules/libraries
import sys
import argparse
import os

## Get user's input 
parser = argparse.ArgumentParser(description= """Produces for each TE insertion two fasta (one for read pairs supporting + clusters and another one for reads supporting - clusters).""")
parser.add_argument('insertionsPath', help='TraFic output file containing the ids of the reads supporting + and - clusters.')
parser.add_argument('fastaPath', help='Fasta containing the read pairs supporting TE insertions for a given sample.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
insertionsPath = args.insertionsPath
fastaPath = args.fastaPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "insertions: ", insertionsPath
print "fasta: ", fastaPath
print "outDir: ", outDir
print 

print "***** Executing ", scriptName, " *****"
print 
print "..."
print 

## Open input files
insertions = open(insertionsPath, 'r')
fasta = open(fastaPath, 'r')

### 1) Read fasta file and store all the information in a dictionary with the 
# read pair id as key and a two elements list (element1: mate 1 sequence; 
# element2: mate 2 sequence) as value

fastaDict = {}
 
## Read fasta file line by line
for line in fasta:
    line = line.rstrip('\n')

    ## Select fasta sequence id lines
    if line.startswith(">"):
        
        ##  Extract pair id and mate info from fasta ids
        readId = line.split(">")[1] 
        pairId, mateId = readId.split("/")
        
        ## Extract the sequence corresponding to the id
        seq = next(fasta)
        seq = seq.rstrip('\n')
        
        ## Inizialize dictionary key if it does not exist
        if pairId not in fastaDict:
            fastaDict[pairId] = [None] * 2
            
	# a) Mate1
        if mateId == "1":
            fastaDict[pairId][0] = seq
        
  	# b) Mate2
        elif mateId == "2":
            fastaDict[pairId][1] = seq
      
    
### 2) Read TraFiC somatic output file and generates a nested dictionary with the following info:
# * Key 1: insertion id (A different insertion id for + and - clusters)
# * Value 1: dictionary
#         - Key 1: Read pair id
#         - Value 2: two elements list (element1: mate 1 sequence; 
#           element2: mate 2 sequence)
# NOTE: Sequences in Value 2 list extracted from the fastaDict

supportingReadsDict = {}
 
## Read insertions file line by line
for line in insertions:
    line = line.rstrip('\n')
    line = line.split("\t")

    ## A) Line with expected number of columns
    if (int(len(line)) == 12):
        chrPlus = line[0]
        begPlus = line[1]
        endPlus = int(line[2]) + 100 # Done because this coordinate is the beginning of the read. So, I need to sum the readlength. I need to add an input parameter to specify the read length. 
        nbReadsPlus = line[3]
        familyPlus = line[4]
        readPairListPlus = line[5].split(",")
        chrMinus = line[6]
        begMinus = line[7]
        endMinus = line[8]
        nbReadsMinus = line[9]
        familyMinus = line[10]
        readPairListMinus = line[11].split(",")
    
        ## Rename "Other" family insertions as SVA 
        # Plus
        if (familyPlus == "Other"):
	    familyPlus = "SVA"
	
        # Minus
        if (familyMinus == "Other"):
	    familyMinus = "SVA"
	
        ## Generate an insertion id for + and - clusters (insertion coordinates defined by the end 
        # of + cluster and beg of - cluster)    
        insertionIdPlus = familyPlus + ":" + chrPlus + "_" + str(endPlus) + "_" + begMinus + ":" + "+"
        insertionIdMinus = familyMinus + ":" + chrMinus + "_" + str(endPlus) + "_" + begMinus + ":" + "-"
    
        ## Inizialize dictionary keys for + and - clusters if they do not exist
        # a) + Cluster
        if insertionIdPlus not in supportingReadsDict:
            supportingReadsDict[insertionIdPlus] = {}
        
        # b) - Cluster    
        if insertionIdMinus not in supportingReadsDict:
            supportingReadsDict[insertionIdMinus] = {}    
    
        ## Add the list with mate 1 and mate 2 sequences as value for + and - clusters:
        # a) + Cluster
        for pairId in readPairListPlus:
            supportingReadsDict[insertionIdPlus][pairId] = fastaDict[pairId]
    
        # b) - Cluster
        for pairId in readPairListMinus:
            supportingReadsDict[insertionIdMinus][pairId] = fastaDict[pairId]

    else:
	print "[ERROR] Input line with unexpected number of columns."
        

### 3) Generate a fasta per insertion and cluster containing the 
# read pairs supporting it

## Iterate over the insertions 
for insertionId in supportingReadsDict:
    
    ## Open output file
    fileName = insertionId + ".fa"
    outFilePath = outDir + "/" + fileName
    outFile = open( outFilePath, "w" )
    
    ##  Iterate over the read pairs supporting the cluster
    for readPair in supportingReadsDict[insertionId]:
        
        ## Write supporting mate 1 and 2 into a fasta
        # a) Mate 1
        mate1Id = ">" + readPair + "/1"
        mate1Seq = supportingReadsDict[insertionId][readPair][0]
        outFile.write("%s\n" % mate1Id)
        outFile.write("%s\n" % mate1Seq)
    
        # b) Mate 2
        mate2Id = ">" + readPair + "/2"    
        mate2Seq = supportingReadsDict[insertionId][readPair][1]
        outFile.write("%s\n" % mate2Id)
        outFile.write("%s\n" % mate2Seq)
            
    # Close output fasta file
    outFile.close()


print "***** Finished! *****"
print 


