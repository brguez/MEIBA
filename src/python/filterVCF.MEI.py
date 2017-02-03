#!/usr/bin/env python
#coding: utf-8


#### CLASSES ####
class MEIcluster():
    """
    """
    def __init__(self, MEIObj):
        """
        """

        self.chrom = MEIObj.chrom
        self.beg = MEIObj.pos - int(MEIObj.infoDict["CIPOS"])  # Extend cluster beg-end interval with the CIPOS
        self.end = MEIObj.pos + int(MEIObj.infoDict["CIPOS"]) 
        self.MEIlist = [ MEIObj ]

    def addMEI(self, MEIObj):
        """
        """
        self.chrom = MEIObj.chrom
        self.beg = MEIObj.pos - int(MEIObj.infoDict["CIPOS"])  # Extend cluster beg-end interval with the CIPOS        
        self.end = MEIObj.pos + int(MEIObj.infoDict["CIPOS"]) 
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


def organizeMEI(MEIList):
        """
        Organize MEI into a nested dictionary containing for each chromosome a list of MEI objects. Structure: 
        
        key1(chrId) -> value(MEIObj)

        MEIObj in the list sorted in increasing coordinates ordering
        """

        MEIDict = {}
 
        ### Build dictionary
        # Per MEI object 
        for MEIObj in MEIList:

            # a) First MEI in a given chromosome
            if MEIObj.chrom not in MEIDict:
                
                MEIDict[MEIObj.chrom] = [ MEIObj ]
     
            # b) There are already MEI in the chromosome
            else:
                MEIDict[MEIObj.chrom].append(MEIObj)

        # Sort MEI in increasing coordinates ordering
        for chrom in MEIDict:
    
            MEIlist = MEIDict[chrom]
            MEIlistSorted = sorted(MEIlist, key=lambda MEIObj: MEIObj.pos)
            
        return MEIDict


def clusterMEI(MEIList):
    """
    """
    
    # 1. Organize MEI
    MEIDict = organizeMEI(MEIList)
    
    # 2. Cluster MEI objects
    #  chr1 --------*---------------*-------------*------------
    #     beg    cluster1        cluster2      cluster3      end
    #         (MEI1,MEI2,MEI3)   (MEI4)       (MEI5,MEI6)       
    clusterList = []

    # For each chromosome:
    for chrom in MEIDict:

        MEIlist = MEIDict[chrom]

        # Per MEI object:        
        for MEIObj in MEIlist:
            print "MEI-info: ", chrom, MEIObj.pos, MEIObj.infoDict["TYPE"], MEIObj.infoDict["CLASS"], MEIObj.infoDict["SCORE"], MEIObj.infoDict["CIPOS"]   
        
            # A) No cluster in the list -> Create first clustert
            if not clusterList:

                clusterObj = MEIcluster(MEIObj)
                clusterList.append(clusterObj)
        
            # B) There is already at least one cluster in the list -> Check if current MEI within the latest cluster
            else:
                lastClusterObj = clusterList[-1]     
                totalCIPOS = int(MEIObj.infoDict["CIPOS"]) + 5 # Allow a minimum of 5 nucleotides of difference (CIPOS=0, X=overhang)
                                                                                                      
                ### Define interval to search for overlap:
                ##         <-----------------------------interval------------------------------->
                # <clusterBeg-totalCIPOS)                                            <clusterEnd+totalCIPOS>
                #         beg                                                                  end
                chrom = lastClusterObj.chrom
                beg = lastClusterObj.beg - totalCIPOS
                end = lastClusterObj.end + totalCIPOS

                print "cluster_range: ", lastClusterObj.chrom, lastClusterObj.beg, lastClusterObj.end, beg, end

                ## A) Current MEI within previous cluster interval -> add MEI to the cluster                                
                if (chrom == MEIObj.chrom) and (beg <= MEIObj.pos) and (MEIObj.pos <= end):
                    msg = "MEI within cluster"
                    log("DUPLICATES", msg) 
                    lastClusterObj.addMEI(MEIObj)   
                 
                ## B) Current MEI outside previous cluster interval -> create new cluster and add it into the list
                else:
                    msg = "MEI outside the cluster"
                    log("DUPLICATES", msg) 
                    clusterObj = MEIcluster(MEIObj)
                    clusterList.append(clusterObj)
         
            print "**************"

    return clusterList     


def findDuplicates(MEIList):
    """
    There are duplicated insertions since TraFiC sometimes calls more than 1 time the same insertion but with slightly different coordinates. 
    This happens because a single positive, a negative cluster or both is disgretated and identified as two different clusters. This happens when there is not 
    overlap between all the cluster reads. Ejem:

    ----->
       ------>
     ------>
                   ------>
                      ------>
                    ------->       
     <-------->   <---------->
      cluster1      cluster2  
     <-------------------->
            cluster
 
    There are several duplication possibilities:
      - TD0/TD1. Select td1  ** Most common duplication event according to jose 
      - TD0/TD2. Select td2
      - TD1/TD2. ""         

    There also could be triplicated events:
      - TD0/TD1/TD2. Select td1

    I don't know if this kind of duplications can happen...   
      - TD0/TD0. ¿Select the one with more score and if same score use the one with more supporting reads?
      - TD1/TD1. ""
      - TD2/TD2. ""
    
    """
    
    ### Make list with the identifier of the duplicated MEI
    # Identifier in the format
    # sourceElementId = MEIObj.chrom + ':' + str(MEIObj.pos) + '-' + str(end)    
    dupList = []
    
    ### Cluster MEI
    clusterList = clusterMEI(MEIList)

    ### For each MEI cluster
    for cluster in clusterList:
        
        nbMEI = len(cluster.MEIlist)

        ## A) Select clusters composed by more than 1 MEI
        if (nbMEI > 1):
                  
            MEI1 = cluster.MEIlist[0]
            MEI2 = cluster.MEIlist[1]
    
            pos1 = MEI1.pos
            pos2 = MEI2.pos
            type1 = MEI1.infoDict["TYPE"] 
            type2 = MEI2.infoDict["TYPE"] 

            ## 2 MEI composing the cluster 
            if (nbMEI == 2):

                msg = "Cluster composed by 2 MEI: " + str(nbMEI) + " " + MEI1.chrom + " " + str(pos1) + " " +  str(pos2) + " " + type1 + " " + type2
                log("FIND-DUP", msg)                
                typeList = [type1, type2]

                ## a) TD0/TD1. Select td1 (Most common duplication type)
                if all(i in [ 'TD0', 'TD1' ] for i in typeList):
                    
                    print "***** TD0/TD1"

                    ## First MEI is not TD1    
                    if (type1 != 'TD1'):
                        MEIid = MEI1.chrom + ':' + str(MEI1.pos) + ':' + type1
                        
                    ## Second MEI is not TD1
                    else:
                        MEIid = MEI2.chrom + ':' + str(MEI2.pos) + ':' + type2

                    msg = "MEI flagged as duplicate: " + MEIid
                    log("FIND-DUP", msg)
                    dupList.append(MEIid)
                     
                ## b) TD0/TD2. Select td2
                elif all(i in [ 'TD0', 'TD2' ] for i in typeList):
                    
                    print "******** TD0/TD2"

                    ## First MEI is not TD2        
                    if (type1 != 'TD2'):
                        MEIid = MEI1.chrom + ':' + str(MEI1.pos) + ':' + type1
                        
                    ## Second MEI is not TD2
                    else:
                        MEIid = MEI2.chrom + ':' + str(MEI2.pos) + ':' + type2

                    msg = "MEI flagged as duplicate: " + MEIid
                    log("FIND-DUP", msg)
                    dupList.append(MEIid)

                ## c) TD1/TD2. Select td1  
                elif all(i in [ 'TD1', 'TD2' ] for i in typeList):

                    print "******** TD1/TD2"
            
                    ## First MEI is not TD1   
                    if (type1 != 'TD1'):
                        MEIid = MEI1.chrom + ':' + str(MEI1.pos) + ':' + type1
                        
                    ## Second MEI is not TD1
                    else:
                        MEIid = MEI2.chrom + ':' + str(MEI2.pos) + ':' + type2

                    msg = "MEI flagged as duplicate: " + MEIid
                    log("FIND-DUP", msg)
                    dupList.append(MEIid)
                
                ## d) Unexpected combination
                else:
                    msg = "Warning - unexpected combination of types"
                    log("FIND-DUP", msg)                    
                
                print "dupList: ", dupList

            ## 3 MEI composing the cluster
            elif (nbMEI == 3):

                MEI3 = cluster.MEIlist[2]
                pos3 = MEI3.pos
                type3 = MEI3.infoDict["TYPE"] 

                msg = "Cluster composed by 3 MEI: " + str(nbMEI) + " " + MEI1.chrom + " " + str(pos1) + " " +  str(pos2) + " " +  str(pos3)  + " " + type1 + " " + type2 + " " + type3
                log("FIND-DUP", msg) 
                
                typeList = [type1, type2, type3]
                
                ## a) TD0/TD1/TD2. Select td1
                if all(i in [ 'TD0', 'TD1', 'TD2' ] for i in typeList):
                    
                    print "******** TD0/TD1/TD2"
                
                    ## First MEI is TD1   
                    if (type1 == 'TD1'):
                        MEIid2 = MEI2.chrom + ':' + str(MEI2.pos) + ':' + type2
                        MEIid3 = MEI3.chrom + ':' + str(MEI3.pos) + ':' + type3
                        MEIidList = [MEIid2, MEIid3]

                    ## Second MEI is TD1
                    elif (type2 == 'TD1'):                   
                        MEIid1 = MEI1.chrom + ':' + str(MEI1.pos) + ':' + type1
                        MEIid3 = MEI3.chrom + ':' + str(MEI3.pos) + ':' + type3
                        MEIidList = [MEIid1, MEIid3]
                    
                    ## Third MEI is TD1
                    else:
                        MEIid1 = MEI1.chrom + ':' + str(MEI1.pos) + ':' + type1
                        MEIid2 = MEI2.chrom + ':' + str(MEI2.pos) + ':' + type2
                        MEIidList = [MEIid1, MEIid2]

                    msg = "2 MEI flagged as duplicate: " + str(MEIidList)
                    log("FIND-DUP", msg)

                    ## Add the two duplicated MEI to the list of duplicates
                    dupList = dupList + MEIidList

                ## b) Unexpected combination
                else:
                    msg = "Warning - unexpected combination of types"
                    log("FIND-DUP", msg) 

            ## 3 MEI composing the cluster
            else:
                msg = "Warning - more than 3 MEI composing the cluster"
                log("FIND-DUP", msg) 
                        
        ## B) Cluster composed by a single MEI
        else:
              msg = "Cluster composed by a single MEI: " + str(nbMEI)
              log("FIND-DUP", msg)   
    
    return dupList

def scoreFilter(score, minScore):
    """
    Filter mobile element insertions based on the reconstruction score. 
    If score < minScore_threshold; then filter out the insertion
    """

    msg = "score,minScore: " + str(score) + " " + str(minScore)
    log("SCORE", msg)  

    # a) Insertion do not passing the filter
    if (score < minScore):  
        filterStatus = "SCORE"

    # b) Insertion passing the filter
    else:
        filterStatus = "."

    return filterStatus


def repeatsFilter(RTclass, repeat, divergence, maxDiv):
    """    
    Filter out those solo insertions and transductions that fulfill at least one of these two conditions:

    1. Insertion overlapping a repetitive element of the same class with a millidivergence < maxDiv_threshold. Notes:
        - Millidivergence is a measure of the degree betweene a given repetitive element and a consensus reference sequence.
        - Goes from 1 (0.1%) to 1000 (100%)
        - It is the opposite of sequence identity.
        - Elements with low divergence to consensus sequence are young very repetitive elements in our genome and thus, 
        prone to missalignments what leads to false positive MEI calls
    2. Insertion overlapping one of the possible satellite regions in the database (ALR/Alpha, BSR/Beta, HSATII)
    """

    msg = "RTclass,repeat,divergence,maxDiv: " + str(RTclass) + " " + str(repeat) + " " + str(divergence) + " " + str(maxDiv)
    log("REPEATS", msg)  

    # A) Insertion overlapping a repetitive element of the same type with a millidivergence < threshold 
    if (RTclass == repeat) and (int(divergence) < maxDiv):

        filterStatus = "REP"
    
    # B) Insertion overlapping a satellite region
    elif (repeat == "ALR/Alpha") or (repeat == "BSR/Beta") or (repeat == "HSATII"):

        filterStatus = "REP"

    # C) Insertion do not fulfilling any of these conditions
    else:
        
        filterStatus = "."    
         
    return filterStatus


def duplFilter(MEIObj, dupList):
    """    
    Filter out insertions included in the duplicate insertions list
    """
    MEIid = MEIObj.chrom + ':' + str(MEIObj.pos) + ':' + MEIObj.infoDict["TYPE"]    
    filterStatus = "DUP" if MEIid in dupList else '.'

    return filterStatus   

#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import os.path
import formats

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('--min-score', default=5, dest='minScore', type=int, help='Minimum insertion score for L1, Alu and SVA. Default: 5, both breakpoints expected to be assembled (5-prime and 3-prime).' )
parser.add_argument('--min-score-ERVK', default=3, dest='minScoreERVK', type=int, help='Minimum insertion score for ERVK. Note that mechanism of retrotransposition different, no polyA, so an ERVK insertion can not have an score > 3. Default: 3.' )
parser.add_argument('--max-divergence', default=300, dest='maxDiv', type=int, help='Maximum millidivergence. Default: 300.' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
donorId = args.donorId
minScore = args.minScore
minScoreERVK = args.minScoreERVK
maxDiv = args.maxDiv
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "donorId: ", donorId
print "minScore: ", minScore
print "minScoreERVK: ", minScoreERVK
print "maxDiv: ", maxDiv
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

outFilePath = outDir + '/' + donorId + ".filtered.vcf"

#### 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)


#### 2. Find duplicated insertions
dupList = findDuplicates(VCFObj.lineList)

print "number_duplicates: ", len(dupList), dupList


#### 3. Filter MEI
# Iterate over each MEI in the VCF
for VCFlineObj in VCFObj.lineList:

    failedFiltersList = []
    insertionType = VCFlineObj.infoDict["TYPE"]
    RTclass = VCFlineObj.infoDict["CLASS"] if "CLASS" in VCFlineObj.infoDict else 'NA'

    msg ="Filter " + insertionType + ":" + RTclass + ":" + VCFlineObj.chrom + "_" + str(VCFlineObj.pos)
    subHeader(msg) 

    ### 3.1 Apply score filter:
    score = int(VCFlineObj.infoDict["SCORE"])

    msg = "Apply score filter"
    log("SCORE", msg)

    # a) ERVK insertion -> apply specific score filtering threshold 
    if (RTclass == "ERVK"):
        
        filterStatus = scoreFilter(score, minScoreERVK)

    # b) L1 (TD0, TD1 or TD2), Alu(TD0), SVA(TD0) or pseudogene insertion (PSD)
    else:
        
        filterStatus = scoreFilter(score, minScore)
 
    if (filterStatus != '.'):
        failedFiltersList.append(filterStatus) 

    msg = "Filtering status: " + str(failedFiltersList)
    log("SCORE", msg)    
    
    ### 3.2 Repeats filter:
    msg = "Apply repeats filter"
    log("REPEATS", msg)

    repeat = VCFlineObj.infoDict["REP"] if 'REP' in VCFlineObj.infoDict else '-'
    divergence = VCFlineObj.infoDict["DIV"] if 'DIV' in VCFlineObj.infoDict else '-'    
    
    # a) MEI is a solo or a transduction
    if ((insertionType == "TD0") or (insertionType == "TD1") or (insertionType == "TD2")):
    
        msg = "MEI is a solo or a transduction -> Apply filter"
        log("REPEATS", msg)
        filterStatus = repeatsFilter(RTclass, repeat, divergence, maxDiv)
    
    else:
        filterStatus = '.'
        msg = "MEI is a pseudogene -> Do not apply filter"
        log("REPEATS", msg)
    
    if (filterStatus != '.'):
        failedFiltersList.append(filterStatus) 

    msg = "Filtering status: " + str(failedFiltersList)
    log("REPEATS", msg)  


    ### 3.3 Duplicated insertions filter
    msg = "Apply duplicated insertions filter"
    log("DUPL", msg)
    
    filterStatus = duplFilter(VCFlineObj, dupList)

    if (filterStatus != '.'):
        failedFiltersList.append(filterStatus) 

    msg = "Filtering status: " + str(failedFiltersList)
    log("DUPL", msg) 

    ### 3.4 Set filter VCF field:
    if (len(failedFiltersList) == 0):
        msg = "Insertion pass all the filters"
        log("STATUS", msg)  
        VCFlineObj.filter = "PASS"
    
    else:
        filterField = ";".join(failedFiltersList)
        msg = "Insertion failed the following filters: " + filterField
        log("STATUS", msg) 
        VCFlineObj.filter = filterField


#### 4. Make output VCF
## 4.1 Write header
VCFObj.write_header(outFilePath)

## 4.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print
print "***** Finished! *****"
print
