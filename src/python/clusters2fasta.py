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
parser.add_argument('genome', help='Reference genome in fasta format')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
insertionsPath = args.insertionsPath
fastaPath = args.fastaPath
genome = args.genome
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "insertions: ", insertionsPath
print "fasta: ", fastaPath
print "genome: ", genome
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print


### 0) Make list with chromosome ids.
chrIdsList = []

genome = open(genome, 'r')
for line in genome:
    line = line.rstrip('\n')

    ## Select fasta sequence id lines
    if line.startswith(">"):
        # drop the ">"
        header = line[1:]

        # drop the info
        chrom = header.split(" ")[0]

        # Add chromosome to the list
        chrIdsList.append(chrom)

### 1) Read fasta file and store all the information in a dictionary with the
# read pair id as key and a two elements list (element1: mate 1 sequence;
# element2: mate 2 sequence) as value

## Open input files
insertions = open(insertionsPath, 'r')
fasta = open(fastaPath, 'r')

## Read fasta file line by line
fastaDict = {}
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


### 2) Read TraFiC output file and generates a nested dictionary with the following info:
# * Key 1: insertion id (A different insertion id for + and - clusters)
# * Value 1: dictionary
#         - Key 1: Read pair id
#         - Value 2: two elements list (element1: mate 1 sequence;
#           element2: mate 2 sequence)
# NOTE: Sequences in Value 2 list extracted from the fastaDict

supportingReadsDict = {}

## Read insertions file line by line
for line in insertions:

    # ignore comment lines (e.g. header)
    if line.startswith('#'):
        continue
        
    line = line.rstrip('\n')
    fieldsList = line.split("\t")

    ## A) Line with expected number of columns
    if (int(len(fieldsList)) == 31):
        chrPlus = fieldsList[0]
        begPlus = fieldsList[1]
        endPlus = int(fieldsList[2]) + 100 # Done because this coordinate is the beginning of the read. So, I need to sum the readlength. I need to add an input parameter to specify the read length.
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

        ### Do more input sanity checks...
        ## A) Insertion in a chromosome not included in the provided reference genome
        if (chrPlus not in chrIdsList) or (chrMinus not in chrIdsList):
            print >>sys.stderr, "[ERROR] Filtering out an insertion in a chromosome not included in the provided reference genome: ", line

        ## B) Insertion with inconsistent number of supporting reads. Number does not match with readId list length
        elif ((nbReadsPlus != "NA") and (int(nbReadsPlus) != len(readPairListPlus))) or ((nbReadsMinus != "NA") and (int(nbReadsMinus) != len(readPairListMinus))):
            print >>sys.stderr, "[ERROR] Filtering out an insertion with inconsistent number of supporting reads. Number does not match with read identifier list length: ", line

        ## C) Insertion with everything ok
        else:
            ## Rename "Other" family insertions as SVA
            # Plus
            if (familyPlus == "Other"):
                familyPlus = "SVA"

            # Minus
            if (familyMinus == "Other"):
                familyMinus = "SVA"

            ## Generate an insertion id for + and - clusters (insertion coordinates defined by the end
            # of + cluster and beg of - cluster)
            insertionIdPlus = familyPlus + ":" + insertionType + ":" + chrPlus + "_" + str(endPlus) + "_" + begMinus + ":" + "+"
            insertionIdMinus = familyMinus + ":" + insertionType + ":" + chrMinus + "_" + str(endPlus) + "_" + begMinus + ":" + "-"
            
            ## a) Duplicated insertion. Already processed insertion with identical begin and end coordinates for + and - clusters. 
            if (insertionIdPlus in supportingReadsDict) and (insertionIdMinus in supportingReadsDict):

                print >>sys.stderr, "[WARNING] Filtering out a duplicated insertion: ", line
        
            ## b) Not duplicated insertion
            else:

                ### a) + Cluster 
                if (readPairListPlus[0] != "NA"):

                    ## Initialize dictionary key
                    supportingReadsDict[insertionIdPlus] = {}

                    ## Add the list with mate 1 and mate 2 sequences as value 
                    for pairId in readPairListPlus:
                        supportingReadsDict[insertionIdPlus][pairId] = fastaDict[pairId]



                ### b) - Cluster
                if (readPairListMinus[0] != "NA"):
                    
                    ## Initialize dictionary key
                    supportingReadsDict[insertionIdMinus] = {}

                    ## Add the list with mate 1 and mate 2 sequences as value 
                    for pairId in readPairListMinus:
                        supportingReadsDict[insertionIdMinus][pairId] = fastaDict[pairId]  

    ## Line without the expected number of columns            
    else:
        print >>sys.stderr, "[ERROR] Filtering out an insertion with unexpected number of columns: ", line
    
    
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
