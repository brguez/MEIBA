#!/usr/bin/env python
#coding: utf-8

#### FUNCTIONS ####
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
    print timeInfo, string

def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string


####### CLASSES #######
class cluster():
    """        
    """
    def __init__(self, chrom, beg, end, nbReads, rtClass, readList):
        """
        """
        self.chrom = chrom
        self.beg = int(beg)
        self.end = int(end)
        self.nbReads = int(nbReads)
        self.rtClass = rtClass
        self.readList = readList

class pairedCluster():
    """       
    """
    def __init__(self, plusClusterObj, minusClusterObj):
        """
        """
        ### Generic fields (TD0, TD1, TD2 and PSD)
        # Positive cluster
        self.chromPlus = plusClusterObj.chrom
        self.begPlus = int(plusClusterObj.beg)
        self.endPlus = int(plusClusterObj.end)
        self.nbReadsPlus = int(plusClusterObj.nbReads)
        self.classPlus = plusClusterObj.rtClass
        self.readListPlus = plusClusterObj.readList

        # Negative cluster
        self.chromMinus = minusClusterObj.chrom
        self.begMinus = int(minusClusterObj.beg)
        self.endMinus = int(minusClusterObj.end)
        self.nbReadsMinus = int(minusClusterObj.nbReads)
        self.classMinus = minusClusterObj.rtClass
        self.readListMinus = minusClusterObj.readList

        self.insertionType = "TD0"

        ### L1 transductions specific fields (TD1 and TD2). "NA" for TD0 and PSD
        self.cytobandId = "NA"        
        self.sourceType = "NA"           
        self.chromSource = "NA"           
        self.begSource = "NA"           
        self.endSource = "NA"            
        self.strandSource = "NA"           
        self.tdBeg = "NA"             
        self.tdEnd = "NA"            
        self.tdRnaLen = "NA"           
        self.tdLen = "NA"            

        ### Processed pseudogene specific fields (PSD). "NA" for TD0, TD1 and TD2
        self.psdGene = "NA"          
        self.chromExonA = "NA"           
        self.begExonA = "NA"           
        self.endExonA = "NA"           
        self.chromExonB = "NA"           
        self.begExonB = "NA"            
        self.endExonB = "NA"          

        ### L-mediated genomic rearrangement specific fields. "NA" for standard L1, Alu, SVA, ERVK and PSD insertions
        self.grType = "NA"

    def determine_rgType(self):
        """
        """
        dist = self.endPlus - self.begMinus

        # Clusters begin and end at the same position
        if (dist == 0):
            rgType = "NA"
    
        # Positive cluster ends before negative cluster begin <- Large target site deletion
        elif (dist < 0):
            rgType = "DEL"

        # Positive cluster ends after negative cluster begin <- Large target site duplication
        else:
            rgType = "DUP"
        
        return rgType


class unpairedClusters():
    """     
    """
    def __init__(self):
        """
        """
        self.clustersDict = {}
       
    def read_clusters(self, inputPath):
        """
        """
        inputFile = open(inputPath, 'r')
        info("Read input clusters file: " + inputPath)
        
        # Per iteration, create cluster object and add it to the dictionary
        for line in inputFile:
            line = line.rstrip('\n')
            line = line.split("\t")

            chrom, beg, end, nbReads, rtClass, readList = line             
            clusterObj = cluster(chrom, beg, end, nbReads, rtClass, readList)
         
            # print "clusterObj: ", clusterObj.chrom, clusterObj.beg, clusterObj.end, clusterObj.nbReads, clusterObj.rtClass, clusterObj.readList

            # A) First cluster of a given class
            if clusterObj.rtClass not in self.clustersDict:
                self.clustersDict[clusterObj.rtClass] = {}
                self.clustersDict[clusterObj.rtClass][clusterObj.chrom] = [ clusterObj ]

            # B) There are already clusters of this class
            else:

                # a) First cluster in the chromosome
                if clusterObj.chrom not in self.clustersDict[clusterObj.rtClass]:
                    self.clustersDict[clusterObj.rtClass][clusterObj.chrom] = {}
                    self.clustersDict[clusterObj.rtClass][clusterObj.chrom] = [ clusterObj ]

                # b) There are already clusters in the chromosome
                else:
                    self.clustersDict[clusterObj.rtClass][clusterObj.chrom].append(clusterObj)
        
    def filter_clusters(self, chromList):
        """
        """

        filteredClustersDict = {}

        ## For each rtClass
        for rtClass in self.clustersDict:
            
            filteredClustersDict[rtClass] = {}

            ## For each chrom
            for chrom in self.clustersDict[rtClass]:

                ## Chromosome in the target chromosomes list 
                if chrom in chromList:
                    filteredClustersDict[rtClass][chrom] = {}
                    filteredClustersDict[rtClass][chrom] = self.clustersDict[rtClass][chrom]

        return filteredClustersDict
                
    def sort_clusters(self):
        """
        """

        sortedClustersDict = {}

        ## For each rtClass
        for rtClass in self.clustersDict:
            
            sortedClustersDict[rtClass] = {}

            ## For each chrom
            for chrom in self.clustersDict[rtClass]:

                sortedClustersDict[rtClass][chrom] = {}

                clustersObjList = self.clustersDict[rtClass][chrom]
                clustersObjList.sort(key=lambda cluster: cluster.beg, reverse=False)

                sortedClustersDict[rtClass][chrom] = clustersObjList
                

        return sortedClustersDict


####### FUNCTIONS #######
class pairedClusters():
    """        
    """
    def __init__(self):
        """
        """
        self.pairedClustersDict = {}

    def pair_clusters(self, plusClustersObj, minusClustersObj, minDist, maxDist):
        """
        """
        pairedClustersDict = {}
        
        #### Iterate over positive unpaired clusters
        ## For each rtClass
        for rtClass in plusClustersObj.clustersDict:
          
            ## For each chrom
            for chrom in plusClustersObj.clustersDict[rtClass]:

                ## ---- BEGIN-POSITIVE ----
                ## for each positive unpaired cluster:
                for plusClusterObj in plusClustersObj.clustersDict[rtClass][chrom]:
                    #print "**** Positive Cluster: ", rtClass, chrom, plusClusterObj.end

                    # Initialize 
                    nbPairs = 0
                    reciprocalMinusClusterObj = ""
    
                    ## ---- BEGIN-NEGATIVE ----
                    ## If there are negative clusters of the same class and in the same chromosome:
                    if (rtClass in minusClustersObj.clustersDict) and (chrom in minusClustersObj.clustersDict[rtClass]):

                        ## Iterate over negative unpaired clusters of the same class and chromosome as the positive cluster                
                        for minusClusterObj in minusClustersObj.clustersDict[rtClass][chrom]:

                            #print "* Negative Cluster: ", rtClass, chrom, minusClusterObj.end
                            dist =  plusClusterObj.end - minusClusterObj.beg
                            distAbs = abs(dist)

                            ## Negative and positive cluster within range
                            if (minDist < distAbs) and (distAbs < maxDist):
                    
                                nbPairs+=1
                                reciprocalMinusClusterObj = minusClusterObj
                                # print "paired-candidate: ", chrom, plusClusterObj.end - 500, minusClusterObj.beg + 500, rtClass, dist, nbPairs  
  
                    ## ---- END-NEGATIVE ----
            
                    ## Positive cluster is unambiguosly paired with a negative cluster
                    if (nbPairs == 1):
                        #print "paired-cluster: ", chrom, plusClusterObj.end, reciprocalMinusClusterObj.beg, rtClass, nbPairs   
                
                        ## Create paired cluster object
                        pairedClusterObj = pairedCluster(plusClusterObj, reciprocalMinusClusterObj) 
                        pairedClusterObj.grType = pairedClusterObj.determine_rgType()
                        
                        ## Add paired cluster object to the dictionary      
                        # A) First paired cluster of a given class
                        if rtClass not in pairedClustersDict:
                            pairedClustersDict[rtClass] = {}
                            pairedClustersDict[rtClass][chrom] = [ pairedClusterObj ]
                            
                        # B) There are already paired clusters of this class
                        else:
                                    
                            # a) First pairedcluster in the chromosome
                            if chrom not in pairedClustersDict[rtClass]:
                                pairedClustersDict[rtClass][chrom] = {}
                                pairedClustersDict[rtClass][chrom] = [ pairedClusterObj ]

                            # b) There are already paired clusters in the chromosome
                            else:
                                pairedClustersDict[rtClass][chrom].append(pairedClusterObj)

        return pairedClustersDict

    def sort_clusters(self):
        """
        """

        sortedPairedClustersDict = {}

        ## For each rtClass
        for rtClass in self.pairedClustersDict:
            
            sortedPairedClustersDict[rtClass] = {}

            ## For each chrom
            for chrom in self.pairedClustersDict[rtClass]:

                sortedPairedClustersDict[rtClass][chrom] = {}

                clustersObjList = self.pairedClustersDict[rtClass][chrom]
                clustersObjList.sort(key=lambda cluster: cluster.begPlus, reverse=False)

                sortedPairedClustersDict[rtClass][chrom] = clustersObjList
                
        return sortedPairedClustersDict

    def write_clusters(self, outPath, chromList):
        """
        """
        outFile = open(outPath, 'w')

        # Print header into the output file
        header = "#chromPlus" + "\t" + "begPlus" + "\t" + "endPlus" + "\t" + "nbReadsPlus" + "\t" + "classPlus" + "\t" + "readListPlus" + "\t" + "chromMinus" + "\t" + "begMinus" + "\t" + "endMinus" + "\t" + "nbReadsMinus" + "\t" + "classMinus" + "\t" + "readListMinus" + "\t" + "insertionType" + "\t" + "cytobandId" + "\t" + "sourceType" + "\t" + "chromSource" + "\t" + "begSource" + "\t" + "endSource" + "\t" + "strandSource" + "\t" + "tdBeg" + "\t" + "tdEnd" + "\t" + "tdRnaLen" + "\t" + "tdLen" + "\t" + "psdGene" + "\t" + "chromExonA" + "\t" + "begExonA" + "\t" + "endExonA" + "\t" + "chromExonB" + "\t" + "begExonB" + "\t" + "endExonB" + "\t" + "grType" + "\n"

        outFile.write(header) 

        ## For each rtClass
        for rtClass in self.pairedClustersDict:
            
            ## For each chrom
            for chrom in chromList:
        
                ## There are pairedClusters in this chrom
                if chrom in self.pairedClustersDict[rtClass]:
                    
                    ## For each paired cluster in the chrom
                    for clusterObj in self.pairedClustersDict[rtClass][chrom]:

                        ## Write paired cluster in the output file
                        row = clusterObj.chromPlus + "\t" + str(clusterObj.begPlus) + "\t" + str(clusterObj.endPlus) + "\t" + str(clusterObj.nbReadsPlus) + "\t" + clusterObj.classPlus + "\t" + clusterObj.readListPlus + "\t" + clusterObj.chromMinus + "\t" + str(clusterObj.begMinus) + "\t" + str(clusterObj.endMinus) + "\t" + str(clusterObj.nbReadsMinus) + "\t" + clusterObj.classMinus + "\t" + clusterObj.readListMinus + "\t" + clusterObj.insertionType + "\t" + clusterObj.cytobandId + "\t" + clusterObj.sourceType + "\t" + clusterObj.chromSource + "\t" + clusterObj.begSource + "\t" + clusterObj.endSource + "\t" + clusterObj.strandSource + "\t" + clusterObj.tdBeg + "\t" + clusterObj.tdEnd + "\t" + clusterObj.tdRnaLen + "\t" + clusterObj.tdLen + "\t" + clusterObj.psdGene + "\t" + clusterObj.chromExonA + "\t" + clusterObj.begExonA + "\t" + clusterObj.endExonA + "\t" + clusterObj.chromExonB + "\t" + clusterObj.begExonB + "\t" + clusterObj.endExonB + "\t" + clusterObj.grType + "\n"        

                        outFile.write(row) 

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import re
from operator import itemgetter, attrgetter, methodcaller


## Get user's input ##
parser = argparse.ArgumentParser(description="")
parser.add_argument('plusPath', help='Not paired positive clusters tsv file')
parser.add_argument('minusPath', help='Not paired negative clusters tsv file')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('--min-dist', default=0, dest='minDist', type=int, help='Minimum distance between positive and negative clusters. Default 200bp.' )
parser.add_argument('--max-dist', default=200, dest='maxDist', type=int, help='Maximum distance between positive and negative clusters. Default 200bp.' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
plusPath = args.plusPath
minusPath = args.minusPath
fileName = args.fileName
minDist = args.minDist
maxDist = args.maxDist
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "plusPath: ", plusPath
print "minusPath: ", minusPath
print "fileName: ", fileName
print "minDist: ", minDist
print "maxDist: ", maxDist
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ##

## 1. Initialize unpaired clusters objects
###########################################
header("1. Initialize unpaired clusters objects")
unpairedPlusClustersObj = unpairedClusters()
unpairedMinusClustersObj = unpairedClusters()

#print unpairedPlusClustersObj, unpairedMinusClustersObj


## 2. Read unpaired clusters files
###################################
header("2. Read unpaired clusters files")
unpairedPlusClustersObj.read_clusters(plusPath)
unpairedMinusClustersObj.read_clusters(minusPath)


## 3. Filter unpaired clusters to select only those clusters in standard chromomosomes
#######################################################################################
header("3. Filter unpaired clusters to select only those clusters in standard chromomosomes")

targetChromList = list(range(1, 23))
targetChromList = [str(i) for i in targetChromList]
targetChromList.append("X")

unpairedPlusClustersObj.clustersDict = unpairedPlusClustersObj.filter_clusters(targetChromList)
unpairedMinusClustersObj.clustersDict = unpairedMinusClustersObj.filter_clusters(targetChromList)


## 4. Sort unpaired clusters list in the chromosomes
#####################################################
# Increasing coordinates order
header("4. Sort unpaired clusters list in the chromosomes")
unpairedPlusClustersObj.clustersDict = unpairedPlusClustersObj.sort_clusters()
unpairedMinusClustersObj.clustersDict = unpairedMinusClustersObj.sort_clusters()

#for clusterObj in unpairedPlusClustersObj.clustersDict["L1"]["1"]:
#    print "cluster-beg: ", clusterObj.beg


## 5. Pair clusters
####################
header("5. Pair clusters")

## Create paired clusters object
pairedClustersObj = pairedClusters()

pairedClustersObj.pairedClustersDict = pairedClustersObj.pair_clusters(unpairedPlusClustersObj, unpairedMinusClustersObj, minDist, maxDist)

## Sort paired clusters
pairedClustersObj.pairedClustersDict = pairedClustersObj.sort_clusters() 

## 6. Write paired clusters into the output file
#################################################
header("6. Write paired clusters into the output file")
outPath = outDir + "/" + fileName + ".tsv"

pairedClustersObj.write_clusters(outPath, targetChromList)


## End ##
print
print "***** Finished! *****"
print

