#!/usr/bin/env python
#coding: utf-8


#### CLASSES ####
class MEIcluster():
    """
    """
    def __init__(self, MEIObj, offset):
        """
        Initialize MEI cluster:

                       cluster_range
            <------------------------------------>
            beg [bkpA - offset]        end [(bkpB or bkpA) - offset] 
        """

        self.chrom = MEIObj.chrom
        self.beg, self.end = insertionRange(MEIObj, offset)        
        self.MEIlist = [ MEIObj ]

    def addMEI(self, MEIObj, offset):
        """
        Add MEI to the cluster and extend cluster end coordinates till the new MEI
        
                      cluster_old_range                 
            <------------------------------------>--------------->
            oldBeg                            oldEnd          newEnd 

            <---------------------------------------------------->
                              cluster_new_range
        """

        beg, end = insertionRange(MEIObj, offset)         
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


def insertionRange(MEIObj, offset):
    """
    Define MEI range as follows:
                
                      range
       <------------------------------------>
       beg [bkpA - offset]        end [(bkpB or bkpA)- offset]  
    """
    bkpA = MEIObj.pos
    bkpB = int(MEIObj.infoDict["BKPB"]) if "BKPB" in MEIObj.infoDict else MEIObj.pos 
    begRange = bkpA - offset
    endRange = bkpB + offset

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
        a dictionary of germline MEI objects organized by coordinates. Structure: 
        
        key1(MEIclass) -> value(dict2) -> key2(chrId) -> value(MEIObjList)

        MEIObjList sorted in increasing coordinates ordering
        """

        MEIDict = {}
 
        ### Build dictionary
        # Per MEI object 
        for MEIObj in MEIList:

            # a) First MEI of a given class    
            if MEIObj.infoDict["CLASS"] not in MEIDict:
            
                MEIDict[MEIObj.infoDict["CLASS"]] = {}
                MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom] = [MEIObj]
            
            # b) First MEI of a given class in the chromosome
            elif MEIObj.chrom not in MEIDict[MEIObj.infoDict["CLASS"]]:
                MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom] = [MEIObj]

            # c) There are already MEI of this class in the chromosome
            else:
                MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom].append(MEIObj)
                
        # Sort MEI list in increasing coordinates ordering
        for MEIclass in MEIDict:
            for chrom in MEIDict[MEIclass]:

                MEIlist = MEIDict[MEIclass][chrom]
                MEIlistSorted = sorted(MEIlist, key=lambda MEIObj: MEIObj.pos)
                MEIDict[MEIclass][chrom] = MEIlistSorted
            
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
    for MEIclass in MEIDict:

        for chrom in MEIDict[MEIclass]:
        
            MEIlist = MEIDict[MEIclass][chrom]

            # Per MEI object:        
            for MEIObj in MEIlist:

                bkpB = int(MEIObj.infoDict["BKPB"]) if "BKPB" in MEIObj.infoDict else "UNK"
        
                msg = "MEI: " + chrom + " " + str(MEIObj.pos) + " " + str(bkpB) + " " + MEIObj.infoDict["TYPE"] + " " + MEIObj.infoDict["CLASS"] + " " + MEIObj.infoDict["SCORE"] 
                log("CLUSTER", msg)    
        
                # A) No cluster in the list -> Create first cluster
                if not clusterList:

                    msg = "Initialize first cluster"
                    log("CLUSTER", msg) 
                    clusterObj = MEIcluster(MEIObj, 5) # Allow up to 5 nucleotides of margin
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
                    begMEIrange, endMEIrange = insertionRange(MEIObj, 5)     
              
                    #### Check if both ranges overlap.
                    overlapping = overlap(begMEIrange, endMEIrange, begClusterRange, endClusterRange) 

                    msg = "cluster_range,MEI_range: " + str(begClusterRange) + " " + str(endClusterRange) + " " + str(begMEIrange) + " " + str(endMEIrange)
                    log("CLUSTER", msg) 

                    ## A) MEI in the same chromosome than the current cluster and overlapping the cluster interval -> add MEI to the cluster                                
                    if (chrom == lastClusterObj.chrom) and (overlapping):
                        msg = "MEI within cluster -> add MEI to the cluster"
                        log("CLUSTER", msg) 
                        lastClusterObj.addMEI(MEIObj, 5) # Allow up to 5 nucleotides of margin 
                     
                    ## B) Current MEI outside previous cluster interval -> create new cluster and add it into the list
                    else:
                        msg = "MEI outside the cluster -> create new cluster "
                        log("CLUSTER", msg) 
                        clusterObj = MEIcluster(MEIObj, 5) # Allow up to 5 nucleotides of margin
                        clusterList.append(clusterObj)
                    
                    msg = "Number of MEI within cluster: ", len(clusterObj.MEIlist)
                    log("CLUSTER", msg) 
                    msg = "----------------------"
                    log("CLUSTER", msg) 

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
    # Identifier in the format:
    dupList = []
    
    ### Cluster MEI
    clusterList = clusterMEI(MEIList)

    ### For each MEI cluster
    for cluster in clusterList:
        
        nbMEI = len(cluster.MEIlist)

        ## Select clusters composed by more than 1 MEI
        if (nbMEI > 1):
                          
            typeList = [MEI.infoDict["TYPE"] for MEI in cluster.MEIlist]
            
            print "typeList: ", typeList

            ## A) All the MEI in the cluster with the same insertion type (there can be 1 to N insertions; N can be infinite)
            if len(set(typeList)) == 1:
                                
                msg = "Cluster composed by " + str(nbMEI) + " MEI of the same insertion type"  
                log("FIND-DUP", msg)   
                
                scoreList = [MEI.infoDict["SCORE"] for MEI in cluster.MEIlist]
            
                print "scoreList: ", scoreList
                
                #### Select those MEI with the highest score and
                # include the others in the duplicates list
                maxScore = 0
                maxScoreList = []

                for MEI in cluster.MEIlist:
                    
                    score = int(MEI.infoDict["SCORE"])
                    
                    # a) Insertion with higher score than the maximum
                    if (score > maxScore):
                        maxScore = score
                        
                        # As the MEI in the maximum score list do not have the maximum score anymore consider them as duplicates 
                        for maxMEI in maxScoreList:
                            MEIid = maxMEI .chrom + '_' + str(maxMEI .pos) + '_' + maxMEI .infoDict["TYPE"] + '_' + maxMEI .genotype
                            dupList.append(MEIid)

                        # Re-initialize the maximum score list
                        maxScoreList = [ MEI ]
                        
                    # b) Insertion with an score equal to the maximum                    
                    elif (score == maxScore):
                        maxScoreList.append(MEI)

                    # c) Insertion with an score lower than the maximum                    
                    else:
                        MEIid = MEI.chrom + '_' + str(MEI.pos) + '_' + MEI.infoDict["TYPE"] + '_' + MEI.genotype
                        dupList.append(MEIid)
            
                print "maxScoreList: ", maxScoreList
                print "dupList: ", dupList
                
                #### If several possible MEI with the highest score -> select the one with the highest number of supporting reads 
                if (len(maxScoreList) > 1):
                    msg = str(len(maxScoreList)) + " MEI with the highest score. Select the one with more supporting reads"  
                    log("FIND-DUP", msg) 
                                                     
                    maxRC = 0
                    maxRCList = []

                    for MEI in maxScoreList:
                    
                        RCP, RCN, sampleId = MEI.genotype.split(":")                   
                        RC = int(RCP) + int(RCN)

                        print MEI.chrom + '_' + str(MEI.pos) + '_' + MEI.infoDict["TYPE"] + '_' + MEI.genotype + ' ' + str(RC)

                        # a) Insertion with higher read count (RC) than the maximum
                        if (RC > maxRC):
                            maxRC = RC
                        
                            # As the MEI in the maximum RC list do not have the maximum RC anymore consider them as duplicates 
                            for maxMEI in maxRCList:
                                MEIid = maxMEI .chrom + '_' + str(maxMEI .pos) + '_' + maxMEI .infoDict["TYPE"] + '_' + maxMEI .genotype
                                dupList.append(MEIid)

                            # Re-initialize the maximum RC list
                            maxRCList = [ MEI ]
                        
                        # b) Insertion with a RC equal to the maximum                    
                        elif (RC == maxRC):
                            maxRCList.append(MEI)

                        # c) Insertion with a RC lower than the maximum                    
                        else:
                            MEIid = MEI.chrom + '_' + str(MEI.pos) + '_' + MEI.infoDict["TYPE"] + '_' + MEI.genotype
                            dupList.append(MEIid)
                
                    print "maxRCList: ", maxRCList
                    print "dupList: ", dupList

                    ## If still several possible MEI (really unlikely) raise a warning
                    if (len(maxRCList) > 1):
                        msg = "Warning - Several possible representative MEI (highest score and total read count)"
                        log("FIND-DUP", msg)  
          
            ## B) 2 MEI with a different insertion type composing the cluster 
            elif (nbMEI == 2):
                
                MEI1 = cluster.MEIlist[0]
                MEI2 = cluster.MEIlist[1]
                pos1 = MEI1.pos
                pos2 = MEI2.pos
                type1 = MEI1.infoDict["TYPE"] 
                type2 = MEI2.infoDict["TYPE"] 
                genotype1 = MEI1.genotype
                genotype2 = MEI2.genotype
    
                msg = "Cluster composed by 2 MEI: " + str(nbMEI) + " " + MEI1.chrom + " " + str(pos1) + " " +  str(pos2) + " " + type1 + " " + type2
                log("FIND-DUP", msg)             

                ## a) TD0/TD1. Select td1 (Most common duplication type)
                if set([ 'TD0', 'TD1' ]) == set(typeList):
                    
                    msg = "TD0/TD1 combination"
                    log("FIND-DUP", msg)

                    ## First MEI is not TD1    
                    if (type1 != 'TD1'):
                        MEIid = MEI1.chrom + '_' + str(MEI1.pos) + '_' + type1 + '_' + genotype1 
                        
                    ## Second MEI is not TD1
                    else:
                        MEIid = MEI2.chrom + '_' + str(MEI2.pos) + '_' + type2 + '_' + genotype2

                    msg = "MEI flagged as duplicate: " + MEIid
                    log("FIND-DUP", msg)

                    # Add the duplicated MEI to the list
                    dupList.append(MEIid)
                     
                ## b) TD0/TD2. Select td2
                elif set([ 'TD0', 'TD2' ]) == set(typeList):

                    msg = "TD0/TD2 combination"
                    log("FIND-DUP", msg)

                    ## First MEI is not TD2        
                    if (type1 != 'TD2'):
                        MEIid = MEI1.chrom + '_' + str(MEI1.pos) + '_' + type1 + '_' + genotype1 
                        
                    ## Second MEI is not TD2
                    else:
                        MEIid = MEI2.chrom + '_' + str(MEI2.pos) + '_' + type2 + '_' + genotype2

                    msg = "MEI flagged as duplicate: " + MEIid
                    log("FIND-DUP", msg)

                    # Add the duplicated MEI to the list
                    dupList.append(MEIid)

                ## c) TD1/TD2. Select td1  
                elif set([ 'TD1', 'TD2' ]) == set(typeList):

                    msg = "TD1/TD2 combination"
                    log("FIND-DUP", msg)

                    ## First MEI is not TD1   
                    if (type1 != 'TD1'):
                        MEIid = MEI1.chrom + '_' + str(MEI1.pos) + '_' + type1 + '_' + genotype1 
                        
                    ## Second MEI is not TD1
                    else:
                        MEIid = MEI2.chrom + '_' + str(MEI2.pos) + '_' + type2 + '_' + genotype2

                    msg = "MEI flagged as duplicate: " + MEIid
                    log("FIND-DUP", msg)

                    # Add the duplicated MEI to the list
                    dupList.append(MEIid)
                
                ## d) 2 elements with unexpected combination of types
                else:
                    msg = "Warning - 2 elements with unexpected combination of types"
                    log("FIND-DUP", msg)                    
                
                print "dupList: ", dupList

            ## C) 3 MEI with a different insertion type composing the cluster (TD0/TD1/TD2)
            elif set([ 'TD0', 'TD1', 'TD2' ]) == set(typeList):

                msg = "TD0/TD1/TD2 combination"
                log("FIND-DUP", msg)

                MEI1 = cluster.MEIlist[0]
                MEI2 = cluster.MEIlist[1]
                MEI3 = cluster.MEIlist[2]
                pos1 = MEI1.pos
                pos2 = MEI2.pos
                pos3 = MEI3.pos
                type1 = MEI1.infoDict["TYPE"] 
                type2 = MEI2.infoDict["TYPE"] 
                type3 = MEI3.infoDict["TYPE"] 
                genotype1 = MEI1.genotype
                genotype2 = MEI2.genotype
                genotype3 = MEI3.genotype       
    
                msg = "Cluster composed by 3 MEI: " + str(nbMEI) + " " + MEI1.chrom + " " + str(pos1) + " " +  str(pos2) + " " +  str(pos3)  + " " + type1 + " " + type2 + " " + type3
                log("FIND-DUP", msg) 
                
                ## a) First MEI is TD1   
                if (type1 == 'TD1'):
                    MEIid2 = MEI2.chrom + '_' + str(MEI2.pos) + '_' + type2 + '_' + genotype2
                    MEIid3 = MEI3.chrom + '_' + str(MEI3.pos) + '_' + type3 + '_' + genotype3
                    MEIidList = [MEIid2, MEIid3]

                ## b) Second MEI is TD1
                elif (type2 == 'TD1'):                   
                    MEIid1 = MEI1.chrom + '_' + str(MEI1.pos) + '_' + type1 + '_' + genotype1
                    MEIid3 = MEI3.chrom + '_' + str(MEI3.pos) + '_' + type3 + '_' + genotype3
                    MEIidList = [MEIid1, MEIid3]
                    
                ## c) Third MEI is TD1
                else:
                    MEIid1 = MEI1.chrom + '_' + str(MEI1.pos) + '_' + type1 + '_' + genotype1 
                    MEIid2 = MEI2.chrom + '_' + str(MEI2.pos) + '_' + type2 + '_' + genotype2
                    MEIidList = [MEIid1, MEIid2]

                msg = "2 MEI flagged as duplicate: " + str(MEIidList)
                log("FIND-DUP", msg)

                ## Add the two duplicated MEI to the list of duplicates
                dupList = dupList + MEIidList
            
            ## D) more than 2 elements with unexpected combination of types
            else:
                msg = "Warning - more than 2 elements with unexpected combination of types"
                log("FIND-DUP", msg) 
                                
        ## Cluster composed by a single MEI
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
        filterStatus = "PASS"

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
        
        filterStatus = "PASS"    
         
    return filterStatus


def dupFilter(MEIObj, dupList):
    """    
    Filter out insertions included in the duplicate insertions list
    """
    MEIid = MEIObj.chrom + '_' + str(MEIObj.pos) + '_' + MEIObj.infoDict["TYPE"] + '_' + MEIObj.genotype 
    print "test: ", MEIid
    filterStatus = "DUP" if MEIid in dupList else 'PASS'

    return filterStatus   


def falseSomaticSrcFilter(VCFlineObj, somaticMEIDict):
    """    
    Filter out false positive transductions. These are transductions mediated by a putative somatic MEI that has not been reported as TD0 in the sample. 
    """
    
    filterStatus = 'FPSOURCE'
    MEItype = VCFlineObj.infoDict["TYPE"]
        
    ## A) Insertion is a partnered or an orphan transduction mediated by a putative somatic source element:
    if ((MEItype == "TD1") or (MEItype == "TD2")) and (("SRCTYPE" in VCFlineObj.infoDict) and (VCFlineObj.infoDict["SRCTYPE"] == "SOMATIC")):
        tdMEIObj = VCFlineObj 
        srcChrom, srcBeg, srcEnd  = tdMEIObj.infoDict["SRC"].split("_")   

        # There are L1 insertions in the same chromosome as the putative somatic source element:
        if (tdMEIObj.infoDict["CLASS"] in somaticMEIDict) and (srcChrom in somaticMEIDict[tdMEIObj.infoDict["CLASS"]]):
        
            # Select all the insertions in the same chromosome where the putative somatic source element is located 
            MEIObjList = somaticMEIDict[tdMEIObj.infoDict["CLASS"]][srcChrom]
 
            ## Search if putative source element is in the TD0 list
            for MEIObj in MEIObjList:
                    
                ## TD0 insertion. Determine if the insertion corresponds to the transduction source element:
                if (MEIObj.infoDict["TYPE"] == "TD0"):
                    soloMEIObj = MEIObj

                    ## Define transduction putative source element insertion range for search for overlap
                    begSrcRange = int(srcBeg) - 5000
                    endSrcRange = int(srcEnd) + 5000
            
                    msg = "begSrcRange,endSrcRange: " + " " + str(begSrcRange) + " " + str(endSrcRange)
                    log("FPSOURCE", msg)  
    
                    ## Define TD0 insertion range for searching for overlap:   
                    begSoloRange, endSoloRange = insertionRange(soloMEIObj, 5000)     
                    
                    msg = "begSoloRange,endSoloRange: " + " " + str(begSoloRange) + " " + str(endSoloRange)
    
                    ## Check if both ranges overlap meaning that the TD0 would correspond to the transduction source element
                    isSourceElement = overlap(begSrcRange, endSrcRange, begSoloRange, endSoloRange) 
                
                    ## a) TD0 corresponds to the source element mediating the transduction 
                    if isSourceElement == True:
                        msg = "overlap: " + " " + str(begSrcRange) + " " + str(endSrcRange) + " " + str(begSoloRange) + " " + str(endSoloRange)
                        log("FPSOURCE", msg)  
                        filterStatus = 'PASS'

                        break

                    ## b) TD0 does not correspond to the source element mediating the transduction
                    else:
                        msg = "not-overlap: " + " " + str(begSrcRange) + " " + str(endSrcRange) + " " + str(begSoloRange) + " " + str(endSoloRange)
                        log("FPSOURCE", msg)  
            
                    
                print "**************" 
            
    ## B) Insertion is not partnered or an orphan transduction mediated by a putative somatic source element.
    # So, automatically pass the filter
    else:
        filterStatus = 'PASS'

    return filterStatus


def germlineFilter(somaticMEIObj, germlineMEIDict):
    """    
    Filter out germline MEI miscalled as somatic. Check if insertion called as somatic is in a database of germline MEI. 
    """

    filterStatus = 'PASS'

    ## A) somatic MEI is in germline database
    if "GERMDB" in somaticMEIObj.infoDict:
        filterStatus = 'GERMLINE'
    
    ## B) somatic MEI is not in germline database. Check if MEI in matched normal VCF if available 
    elif (germlineMEIDict != False):
        
        ### Matched normal VCF available. Check if somatic MEI is included
        ## There are MEI in the germline VCF of the same class as the somatic MEI
        if somaticMEIObj.infoDict["CLASS"] in germlineMEIDict.keys():

            ## There are MEI in the germline VCF in the same chromosome as the somatic MEI
            if somaticMEIObj.chrom in germlineMEIDict[somaticMEIObj.infoDict["CLASS"]].keys():
            
                ## Define somatic insertion range for searching for overlap:    
                begSomaticRange, endSomaticRange = insertionRange(somaticMEIObj, 10)     
            
                msg = "begSomaticRange,endSomaticRange: " + " " + str(begSomaticRange) + " " + str(endSomaticRange)
                log("GERMLINE", msg)  
                msg = "----------------------------" 

                ## For each germline MEI 
                for germlineMEIObj in germlineMEIDict[somaticMEIObj.infoDict["CLASS"]][somaticMEIObj.chrom]:
     
                    ## Define germline insertion range for searching for overlap:    
                    begGermlineRange, endGermlineRange = insertionRange(germlineMEIObj, 10)     
                    
                    msg = "begGermlineRange,endGermlineRange: " + " " + str(begGermlineRange) + " " + str(endGermlineRange)
    
                    ## Check if the somatic and germline ranges overlap
                    isGermlineMEI = overlap(begSomaticRange, endSomaticRange, begGermlineRange, endGermlineRange) 

                    # MEI in matched normal VCF
                    if isGermlineMEI == True:
                        msg = "overlap: " + " " + str(begSomaticRange) + " " + str(endSomaticRange) + " " + str(begGermlineRange) + " " + str(endGermlineRange)
                        log("GERMLINE", msg)  

                        filterStatus = 'GERMLINE'
                        break
    
    return filterStatus


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
parser = argparse.ArgumentParser(description= "Applies a set of filters to a VCF file and set filter field as PASS or as a list with all the filters a given MEI failed to pass")
parser.add_argument('VCF', help='VCF file to be filtered')
parser.add_argument('donorId', help='Donor Id')
parser.add_argument('filters', help='List of filters to be applied out of 7 possible filtering criteria: SCORE, REP, DUP, FPSOURCE, GERMLINE, CLIPPED and MECHANISM.')
parser.add_argument('--score-L1-TD0', default=2, dest='scoreL1_TD0', type=int, help='Minimum assembly score for solo L1 insertions. Default 2.' )
parser.add_argument('--score-L1-TD1', default=2, dest='scoreL1_TD1', type=int, help='Minimum assembly score for L1 partnered transductions. Default 2.' )
parser.add_argument('--score-L1-TD2', default=2, dest='scoreL1_TD2', type=int, help='Minimum assembly score for L1 orphan transductions. Default 2.' )
parser.add_argument('--score-Alu', default=2, dest='scoreAlu', type=int, help='Minimum assembly score for Alu insertions. Default 2.' )
parser.add_argument('--score-SVA', default=2, dest='scoreSVA', type=int, help='Minimum assembly score for SVA insertions. Default 2.' )
parser.add_argument('--score-ERVK', default=2, dest='scoreERVK', type=int, help='Minimum assembly score for ERVK insertions. Default 2.' )
parser.add_argument('--score-PSD', default=2, dest='scorePSD', type=int, help='Minimum assembly score for processed-pseudogene (PSD) insertions. Default 2.' )
parser.add_argument('--min-clipped', default=2, dest='minClipped', type=int, help='Minimum number of clipped reads supporting the insertion. Default: 2.' )
parser.add_argument('--max-divergence', default=100, dest='maxDiv', type=int, help='Maximum millidivergence. Default: 100.' )
parser.add_argument('--mechanism', default='TPRT,EI,DPA', dest='mechanism', type=str, help='List of insertion mechanisms to be taken into account. 3 possible mechanisms: TPRT, EI AND DPA. Default: TPRT,EI,DPA.' )
parser.add_argument('--germline-VCF', default=False, dest='germlineVCF', help=' VCF with germline MEI calls for the same donor. If provided, input insertions are considered to be somatic. Necesary for GERMLINE filtering.' )

parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
donorId = args.donorId
filters = args.filters
mechanism = args.mechanism
scoreL1_TD0 = args.scoreL1_TD0
scoreL1_TD1 = args.scoreL1_TD1
scoreL1_TD2 = args.scoreL1_TD2
scoreAlu = args.scoreAlu
scoreSVA = args.scoreSVA
scoreERVK = args.scoreERVK
scorePSD = args.scorePSD
minClipped = args.minClipped
germlineVCF = args.germlineVCF
maxDiv = args.maxDiv
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "donorId: ", donorId
print "filters: ", filters
print "mechanism: ", mechanism
print "score-L1-TD0: ", scoreL1_TD0
print "score-L1-TD1: ", scoreL1_TD1
print "score-L1-TD2: ", scoreL1_TD2
print "score-Alu: ", scoreAlu
print "score-SVA: ", scoreSVA
print "score-ERVK: ", scoreERVK
print "score-PSD: ", scorePSD
print "min-clipped-reads: ", minClipped
print "germlineVCF: ", germlineVCF
print "maxDiv: ", maxDiv
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

outFilePath = outDir + '/' + donorId + ".filtered.vcf"

#### 0. Create database of germline MEI events 
filterList = filters.split(',')

# Germline filtering flag provided
if 'GERMLINE' in filterList:

    ## A) Normal matched VCF not provided or file does not exist
    if (germlineVCF == False):
        msg = "Matched normal VCF not provided" 
        log("WARNING", msg)
        germlineMEIDict = False

    ## B) Normal VCF provided but does not exist
    elif not (os.path.isfile(germlineVCF)):
        msg = "Matched normal VCF does not exist" 
        log("WARNING", msg)
        germlineMEIDict = False
            
    ## C) Normal VCF provided -> Organize germline MEI into a dictionary
    else:
        germlineVCFObj = formats.VCF()
        germlineVCFObj.read_VCF(germlineVCF)        
        germlineMEIDict = organizeMEI(germlineVCFObj.lineList)
        

#### 1. Create somatic VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

#### 2. Find somatic duplicated insertions
# Duplicated filtering flag provided
if "DUP" in filterList:
    dupList = findDuplicates(VCFObj.lineList)

    print "number_duplicates: ", len(dupList), dupList

#### 3. Organize somatic MEI into a dictionary. 
# False positive somatic source element filtering flag provided
if "FPSOURCE" in filterList:
    somaticMEIDict = organizeMEI(VCFObj.lineList)

#### 4. Filter somatic MEI
# Iterate over each somatic MEI in the VCF
for VCFlineObj in VCFObj.lineList:

    failedFiltersList = []
    insertionType = VCFlineObj.infoDict["TYPE"]
    RTclass = VCFlineObj.infoDict["CLASS"] if "CLASS" in VCFlineObj.infoDict else 'NA'

    msg ="Filter " + insertionType + ":" + RTclass + ":" + VCFlineObj.chrom + "_" + str(VCFlineObj.pos)
    subHeader(msg) 
 
    ### 4.1  Insertion mechanism filter:
    if "MECHANISM" in filterList:
        mechanismList = mechanism.split(',')
        insertionMechanism = VCFlineObj.infoDict["MECHANISM"] if 'MECHANISM' in VCFlineObj.infoDict else 'NA'

        # Failed filter: Mechanism available and not included into the list:
        if (insertionMechanism != "NA") and (insertionMechanism not in mechanismList):
            failedFiltersList.append("MECHANISM") 

        msg = "Filtering status: " + str(failedFiltersList)
        log("MECHANISM", msg)    

    ### 4.2  Score filter:
    if "SCORE" in filterList:
        
        score = int(VCFlineObj.infoDict["SCORE"])

        msg = "Apply score filter"
        log("SCORE", msg)

        ### Select the appropiate score threshold according to the insertion type
        # a) solo L1 insertion
        if (RTclass == "L1") and (insertionType == "TD0"):
        
            msg = "solo L1 insertion" 
            log("SCORE", msg)
            filterStatus = scoreFilter(score, scoreL1_TD0)

        # b) L1 partnered transductions
        elif (RTclass == "L1") and (insertionType == "TD1"):
        
            msg = "L1 partnered transductions" 
            log("SCORE", msg)
            filterStatus = scoreFilter(score, scoreL1_TD1)    
     
        # c) L1 orphan transductions
        elif (RTclass == "L1") and (insertionType == "TD2"):
        
            msg = "L1 orphan transductions" 
            log("SCORE", msg)
            filterStatus = scoreFilter(score, scoreL1_TD2)    

        # d) Alu insertion
        elif (RTclass == "Alu"):
        
            msg = "Alu insertion" 
            log("SCORE", msg)
            filterStatus = scoreFilter(score, scoreAlu)    

        # e) SVA insertion
        elif (RTclass == "SVA"):
            
            msg = "SVA insertion" 
            log("SCORE", msg)
            filterStatus = scoreFilter(score, scoreSVA)

        # d) ERVK insertion
        elif (RTclass == "ERVK"):
        
            msg = "ERVK insertion" 
            log("SCORE", msg)           
            filterStatus = scoreFilter(score, scoreERVK)

        # e) Processed-pseudogene insertion
        elif (insertionType == "PSD"):

            msg = "Processed-pseudogene insertion" 
            log("SCORE", msg)           
            filterStatus = scoreFilter(score, scorePSD)
    
        # f) Unexpected insertion category 
        else:
            msg = "Unexpected insertion category" 
            log("SCORE", msg)   

        # Insertion does not pass the filter
        if (filterStatus != 'PASS'):
            failedFiltersList.append(filterStatus) 

        msg = "Filtering status: " + str(failedFiltersList)
        log("SCORE", msg)    
    

    ### 4.3 Repeats filter:
    if "REP" in filterList:
        
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
            filterStatus = 'PASS'
            msg = "MEI is a pseudogene -> Do not apply filter"
            log("REPEATS", msg)
    
        # Insertion does not pass the filter
        if (filterStatus != 'PASS'):
            failedFiltersList.append(filterStatus) 

        msg = "Filtering status: " + str(failedFiltersList)
        log("REPEATS", msg)  


    ### 4.4 Duplicated insertions filter
    if "DUP" in filterList:
    
        msg = "Apply duplicated insertions filter"
        log("DUP", msg)
    
        filterStatus = dupFilter(VCFlineObj, dupList)

        # Insertion does not pass the filter
        if (filterStatus != 'PASS'):
            failedFiltersList.append(filterStatus) 

        msg = "Filtering status: " + str(failedFiltersList)
        log("DUP", msg) 

    ### 4.5 False positive somatic source element filter   
    if "FPSOURCE" in filterList:
        msg = "False positive somatic source element filter"
        log("FPSOURCE", msg)
    
        filterStatus = falseSomaticSrcFilter(VCFlineObj, somaticMEIDict)

        # Insertion does not pass the filter
        if (filterStatus != 'PASS'):
            failedFiltersList.append(filterStatus)

        msg = "Filtering status: " + str(failedFiltersList)
        log("FPSOURCE", msg) 


    ### 4.6 Germline insertions miscalled as somatic filter
    if "GERMLINE" in filterList:

        msg = "Apply germline filter"
        log("GERMLINE", msg)

        filterStatus = germlineFilter(VCFlineObj, germlineMEIDict)

        # Insertion does not pass the filter
        if (filterStatus != 'PASS'):
            failedFiltersList.append(filterStatus)

        msg = "Filtering status: " + str(failedFiltersList)
        log("DUP", msg) 

    ### 4.7  Number of clipped-reads filter:
    if "CLIPPED" in filterList:
        
        NDP, NDN, NCA, NCB, sampleId = VCFlineObj.genotype.split(":")

        NCA = 0 if NCA == "." else int(NCA)
        NCB = 0 if NCB == "." else int(NCB)
        nbClipped = NCA + NCB
        
        ## Failed filter: Number of clipped < threshold
        if (nbClipped < minClipped):
            failedFiltersList.append("CLIPPED") 

        msg = "Filtering status: " + str(failedFiltersList)
        log("CLIPPED", msg)   


    ###  Set filter VCF field:
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
