#!/usr/bin/env python
#coding: utf-8


#### CLASSES ####
class MEIcluster():
    """
    """
    def __init__(self, MEIObj, overhang):
        """
        Initialize MEI cluster:

                       cluster_range
            <------------------------------------>
            beg [bkpA - CIPOS - overhang]        end [(bkpB or bkpA) - CIPOS - overhang] 
        """

        self.chrom = MEIObj.chrom
        self.beg, self.end = insertionRange(MEIObj, overhang)        
        self.MEIlist = [ MEIObj ]

    def addMEI(self, MEIObj, overhang):
        """
        Add MEI to the cluster and extend cluster end coordinates till the new MEI
        
                      cluster_old_range                 
            <------------------------------------>--------------->
            oldBeg                            oldEnd          newEnd 

            <---------------------------------------------------->
                              cluster_new_range
        """

        beg, end = insertionRange(MEIObj, overhang)         
        self.end = end # Extend cluster end
        self.MEIlist.append(MEIObj)

#### FUNCTIONS ####
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

def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string


def insertionRange(MEIObj, overhang):
    """
    Define MEI range as follows:
                
                      range
       <------------------------------------>
       beg [bkpA - CIPOS - overhang]        end [(bkpB or bkpA) - CIPOS - overhang]  
    """
    bkpA = MEIObj.pos
    bkpB = int(MEIObj.infoDict["BKPB"]) if "BKPB" in MEIObj.infoDict else MEIObj.pos 
    CIPOS = int(MEIObj.infoDict["CIPOS"])
    begRange = bkpA - CIPOS - overhang
    endRange = bkpB + CIPOS + overhang

    return (begRange, endRange)

def overlap(begA, endA, begB, endB):

    """
    Check if both ranges overlap. 2 criteria for defining overlap: 

    ## A) Begin of the range A within the range B         
    #       *beg* <---------range_A---------->                         
    # <---------range_B----------> 
                
    #    *beg* <-------range_A----->
    # <-------------range_B------------------>

    ## B) Begin of the range B within the range A     
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->
            
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>
    """    
       
    # a) Begin of the range A within the range B   
    if ((begA >= begB) and (begA <= endB)):
        overlap = True
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True

    # c) Ranges do not overlapping
    else:
        overlap = False

    return overlap

def organizeMEI(MEIList):
        """
        Organize MEI into a nested dictionary containing for each insertion class and chromosome 
        a dictionary of MEI objects organized by coordinates. Structure: 
        
        key1(event_type) -> value(dict2) -> key2(chrId) -> value(MEIObjList)

        event_type can be L1, Alu, SVA, ERVK or PSD

        MEIObjList sorted in increasing coordinates ordering
        """
        MEIDict = {}
 
        ### Build dictionary
        # Per MEI object 
        for MEIObj in MEIList:

            eventType = MEIObj.infoDict["CLASS"] if "CLASS" in MEIObj.infoDict else MEIObj.infoDict["TYPE"]
        
            print "tiooo: ", eventType            
    
            # a) First MEI of a given class    
            if eventType not in MEIDict:
            
                MEIDict[eventType] = {}
                MEIDict[eventType][MEIObj.chrom] = [MEIObj]
            
            # b) First MEI of a given class in the chromosome
            elif MEIObj.chrom not in MEIDict[eventType]:
                MEIDict[eventType][MEIObj.chrom] = [MEIObj]

            # c) There are already MEI of this class in the chromosome
            else:
                MEIDict[eventType][MEIObj.chrom].append(MEIObj)
                
        # Sort MEI list in increasing coordinates ordering
        for eventType in MEIDict:
            for chrom in MEIDict[eventType]:

                MEIlist = MEIDict[eventType][chrom]
                MEIlistSorted = sorted(MEIlist, key=lambda MEIObj: MEIObj.pos)
                MEIDict[eventType][chrom] = MEIlistSorted
            
        return MEIDict


def clusterMEI(MEIList):
    """
    """
    
    msg = "** CLUSTER MEI **"  
    log("CLUSTER", msg)

    # 1. Organize MEI
    MEIDict = organizeMEI(MEIList)
    
    # 2. Cluster MEI objects
    #  chr1 --------*---------------*-------------*------------
    #           cluster1        cluster2      cluster3     
    #         (MEI1,MEI2,MEI3)   (MEI4)       (MEI5,MEI6)       
    clusterList = []

    # Iterate over the dictionary picking a MEI list in each iteration:
    for eventType in MEIDict:
        for chrom in MEIDict[eventType]:
        
            MEIlist = MEIDict[eventType][chrom]

            # Per MEI object:        
            for MEIObj in MEIlist:

                bkpB = int(MEIObj.infoDict["BKPB"]) if "BKPB" in MEIObj.infoDict else "UNK"
        
                msg = "MEI: " + chrom + " " + str(MEIObj.pos) + " " + str(bkpB) + " " + MEIObj.infoDict["TYPE"] + " " + MEIObj.infoDict["SCORE"] + " " + MEIObj.infoDict["CIPOS"] 
                log("CLUSTER", msg)    
        
                # A) No cluster in the list -> Create first cluster
                if not clusterList:

                    msg = "Initialize first cluster"
                    log("CLUSTER", msg) 
                    clusterObj = MEIcluster(MEIObj, 15) # Allow up to 15 nucleotides of margin
                    clusterList.append(clusterObj)
        
                # B) There is already at least one cluster in the list -> Check if current MEI within the latest cluster
                else:
                    
                    msg = "Check if MEI within latest cluster"
                    log("CLUSTER", msg) 

                    ## Define cluster range for searching for overlap
                    lastClusterObj = clusterList[-1]     
                    begClusterRange = lastClusterObj.beg
                    endClusterRange = lastClusterObj.end

                    ## Define insertion range for searching for overlap:    
                    begMEIrange, endMEIrange = insertionRange(MEIObj, 15)     
              
                    #### Check if both ranges overlap.
                    overlapping = overlap(begMEIrange, endMEIrange, begClusterRange, endClusterRange) 

                    msg = "cluster_range,MEI_range: " + str(begClusterRange) + " " + str(endClusterRange) + " " + str(begMEIrange) + " " + str(endMEIrange)
                    log("CLUSTER", msg) 

    
                    ## A) Overlapping ranges, so current MEI within previous cluster interval -> add MEI to the cluster                                
                    if overlapping:
                        msg = "MEI within cluster -> add MEI to the cluster"
                        log("CLUSTER", msg) 
                        lastClusterObj.addMEI(MEIObj, 15) # Allow up to 15 nucleotides of margin 
                     
                    ## B) Current MEI outside previous cluster interval -> create new cluster and add it into the list
                    else:
                        msg = "MEI outside the cluster -> create new cluster "
                        log("CLUSTER", msg) 
                        clusterObj = MEIcluster(MEIObj, 15) # Allow up to 15 nucleotides of margin
                        clusterList.append(clusterObj)
                    
                    msg = "Number of MEI within cluster: ", len(clusterObj.MEIlist)
                    log("CLUSTER", msg) 
                    msg = "----------------------"
                    log("CLUSTER", msg) 

    return clusterList     


def selectRepresentativeMEI(MEIList):
    """
    Select representative MEI from a list of MEI. The representative MEI will be the one with the highest score. 
    If several possible MEI with highest score. Select the one with higher number of supporting reads. 
    If still several possibilities. Raise a warning and select the first one on the list.  
    """
    
    msg = "** SELECT REPRESENTATIVE MEI **"  
    log("REPRESENTATIVE", msg)

    ## A) Single MEI in the list
    if (len(MEIList) == 1):

        msg = "Cluster composed by a single MEI"  
        log("REPRESENTATIVE", msg)
        reprMEIObj = cluster.MEIlist[0]

    # B) Multiple MEI in the list
    else:
        msg = "Cluster composed by multiple MEI: " +  str(len(cluster.MEIlist))
        log("REPRESENTATIVE", msg)
                

        #### 1. Select MEI with the highest score
        maxScore = 0
        maxScoreList = []
    
        msg = "1. Select MEI with the highest score"  
        log("REPRESENTATIVE", msg)

        ## For each MEI:
        for MEI in cluster.MEIlist:

            score = int(MEI.infoDict["SCORE"])
    
            msg = "MEI: " + MEI.chrom + " " + str(MEI.pos) + " " + MEI.genotype + " " + str(score) 
            log("REPRESENTATIVE", msg)
                    
            # a) Insertion with higher score than the maximum
            if (score > maxScore):
                msg = "MEI with higher score than the maximum" 
                log("REPRESENTATIVE", msg)
                maxScore = score
                maxScoreList = [ MEI ]
                        
            # b) Insertion with an score equal to the maximum                    
            elif (score == maxScore):
                msg = "MEI with an score equal to the maximum" 
                log("REPRESENTATIVE", msg)
                maxScoreList.append(MEI)
    
            # c) Insertion with lower score than the maximum
            else:
                msg = "MEI with lower score than the maximum"
                log("REPRESENTATIVE", msg)
    
        msg = "Number of MEI with the highest score: " +  str(len(maxScoreList)) 
        log("REPRESENTATIVE", msg)
           
        ## A) Single MEI with the highest score -> select MEI as representative 
        if (len(maxScoreList) == 1):
            reprMEIObj = maxScoreList[0]         
 
        ## B) Several possible MEI with the highest score -> 2. select the one with the highest number of supporting reads 
        else: 
            msg = "2. Select MEI with the highest number of supporting reads (discordant paired-ends + clipped reads)"  
            log("REPRESENTATIVE", msg)
            maxRC = 0
            maxRCList = []

            ## For each MEI:
            for MEI in maxScoreList:
                NDP, NDN, NCA, NCB, sampleIds = MEI.genotype.split(":")

                NCA = 0 if NCA == "." else NCA
                NCB = 0 if NCB == "." else NCB

                RC = int(NDP) + int(NDN) + int(NCA) + int(NCB)
    
                msg = "MEI: " + MEI.chrom + " " + str(MEI.pos) + " " + MEI.genotype + " " + str(RC) 
                log("REPRESENTATIVE", msg)
                
                # a) Insertion with higher read count (RC) than the maximum
                if (RC > maxRC):
                    
                    msg = "MEI with higher read count (RC) than the maximum" 
                    log("REPRESENTATIVE", msg)        
                    maxRC = RC
                        
                    # Re-initialize the maximum RC list
                    maxRCList = [ MEI ]
                            
                # b) Insertion with a RC equal to the maximum                    
                elif (RC == maxRC):
                    msg = "MEI with a RC equal to the maximum" 
                    log("REPRESENTATIVE", msg)
                    maxRCList.append(MEI)
    
                # c) Insertion with a lower RC than the maximum
                else:
                    msg = "with a lower RC than the maximum"
                    log("REPRESENTATIVE", msg)
                
            msg = "Number of MEI with the highest RC: " +  str(len(maxRCList))
            log("REPRESENTATIVE", msg)

            ## A) Single MEI with the highest number of supporting reads -> select MEI as representative 
            if (len(maxRCList) == 1):
                reprMEIObj = maxRCList[0]
        
            ## B) Still several possible MEI (really unlikely) raise a warning and pick the first one
            else:
                msg = "Warning - Several possible representative MEI (highest score and total read count)"
                log("REPRESENTATIVE", msg)  
                reprMEIObj = maxRCList[0]

    msg = "Representative MEI: " + " " + reprMEIObj.chrom + " " + str(reprMEIObj.pos) + " " + reprMEIObj.infoDict["SCORE"] + " " + reprMEIObj.genotype
    log("REPRESENTATIVE", msg)  
    return reprMEIObj

    


#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import os.path
import formats
import os.path
from operator import attrgetter

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCFs', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCFs = args.VCFs
donorId = args.donorId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "input: ", VCFs
print "donorId: ", donorId
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Create database with all MEI events
############################################

## First, make list with all the identified MEI across all the provided samples
allMEIlist = []

with open(VCFs) as VCFs:

    ### Process a VCF in each iteration 
    for VCF in VCFs:
        VCF = VCF.rstrip('\n')        
        

        ## 1. Generate VCF object
        VCFObj = formats.VCF()
        VCFObj.read_VCF(VCF)  
        VCFheader = VCFObj.header

        ## Exclude duplicated insertions
        for MEIObj in VCFObj.lineList:
        
            filterList = MEIObj.filter.split(";")     

            if 'DUP' not in filterList:
            
                allMEIlist.append(MEIObj)

## Then organize them into a dictionary
MEIDict = organizeMEI(allMEIlist)


#### 2. Generate a consensus VCF with a non-redundant list of MEI events
##########################################################################

outVCFObj = formats.VCF()              
outVCFObj.header = VCFheader

## For each possible MEI class
for eventType in MEIDict:

    ## For each possible chromosome
    for chrom in MEIDict[eventType]:

        ## 1) Select list of MEI of a given class in a given chromosome 
        msg = "1) Select list of MEI of a given class in a given chromosome "
        header(msg)
        MEIlist = MEIDict[eventType][chrom]

        ## 2) Cluster MEI
        msg = "2) Cluster MEI"
        header(msg)   
        clusterList = clusterMEI(MEIlist)
        
        ## 3) For each MEI cluster
        msg = "3) For each MEI cluster:"
        header(msg)        

        for cluster in clusterList:

            ## 3.1) Select representative MEI 
            msg="3.1) Select representative MEI"            
            subHeader(msg)
            reprMEIObj = selectRepresentativeMEI(cluster.MEIlist)

            ## 3.2) Make list of samples where the insertion is found  
            msg="3.2) Make list of samples where the insertion is found"            
            subHeader(msg)            
            sampleIdList = [MEIObj.genotype.split(":")[4] for MEIObj in cluster.MEIlist]
            newSampleIds = ";".join(sampleIdList)

            ## 3.3) Substitute the sample list field by the new one in the representative MEI
            msg="3.3) Substitute the sample list field by the new one in the representative MEI"            
            subHeader(msg)
            NDP, NDN, NCA, NCB, sampleIds = reprMEIObj.genotype.split(":")
            reprMEIObj.genotype = NDP + ":" + NDN + ":" + NCA + ":" + NCB + ":" + newSampleIds

            ## 3.4) Add the representative MEI to the output VCF
            msg="3.4) Add the representative MEI to the output VCF"            
            subHeader(msg)
            outVCFObj.addLine(reprMEIObj)
              
#### Sort VCF:
outVCFObj.lineList = outVCFObj.sort()

#### Make output VCF
outFilePath = outDir + '/' + donorId + ".vcf"

## Write header
outVCFObj.write_header(outFilePath)

## Write variants
outVCFObj.write_variants(outFilePath)

## End ##
print
print "***** Finished! *****"
print
