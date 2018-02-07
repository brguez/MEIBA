#!/usr/bin/env python
#coding: utf-8

## Import modules ##
import argparse
import time
import sys
import os.path
import re
import formats
import itertools 

#### FUNCTIONS ####

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

def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string


def parse_args():
    """Define and parse command line parameters."""
    parser = argparse.ArgumentParser(description="Per MEI called by TraFiC: 1) Identifies informative contigs spanning 5' and/or 3' insertion ends if possible, 2) Use informative contigs for characterizing MEI in detail (exact breakpoints, length, strand...) and 3) Produce a VCF with the MAI plus all these information")
    parser.add_argument('inputPaths', help='Text file containing, per MEI, the needed files')
    parser.add_argument('sampleId', help='Sample identifier to be incorporated in the SL field of the output VCF. In PCAWG we use the normal_wgs_aliquot_id or tumor_wgs_aliquot_id.')
    parser.add_argument('fileName', help='Output VCF name. In PCAWG we use the submitted_donor_id.')
    parser.add_argument('genome', help='Reference genome in fasta format')
    parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

    args = parser.parse_args()
    return args


def overlap(begA, endA, begB, endB, overhang):

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
    
    begA = max(0, (begA - overhang)) # Use 0 if it becomes negative 
    endA = endA + overhang
    begB = max(0, (begB - overhang)) # Use 0 if it becomes negative 
    endB = endB + overhang
 

    # A) Range A within range B
    #    *beg* <-------range_A----->   
    # <-------------range_B------------------>
    if ((begA >= begB) and (endA <= endB)):

        nbBases = endA - begA + 1
        overlap = "within"

    # B) Range B within range A
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>
    elif ((begB >= begA) and (endB <= endA)):

        nbBases = endB - begB + 1
        overlap = "within"

    # C) Range A partially overlapping range B
    #           *beg* <---------range_A----------> 
    #   <---------range_B---------->    
    elif ((begA > begB) and (begA <= endB)):
      
        nbBases = endB - begA + 1
        overlap = "partial"

    # D) Range B partially overlapping range A 
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->            
    elif ((begB > begA) and (begB <= endA)):

        nbBases = endA - begB + 1
        overlap = "partial"

    # E) Ranges do not overlapping
    else:
        nbBases = 0
        overlap = "not"

    return overlap, nbBases


#### CLASSES ####
class breakpoint():

    def __init__(self):
        self.chrom = ""
        self.pos = ""
        self.polyA = ""
        self.alignmentObj = "" ## Contig alignment demarcating the insertion breakpoint


class chained_alignment():
    """

    """

    def __init__(self, alignmentList):
        """
        """
        self.alignmentList = alignmentList

    def perc_query_covered(self):
        """
        """
        alignmentObjA = self.alignmentList[0]
        alignmentObjB = self.alignmentList[len(self.alignmentList) - 1]
    
        percQueryCovered = float(alignmentObjB.qEnd - alignmentObjA.qBeg) / alignmentObjA.qSize * 100  
        return percQueryCovered

    def rev_complement(self):
        """
        """

        alignmentListRev = []

        # Iterate over list from the back   
        for alignmentObj in reversed(self.alignmentList):

            #print "** BLAT-BEFORE: ", alignmentObj.qName, alignmentObj.qSize, alignmentObj.qBeg, alignmentObj.qEnd, alignmentObj.tName, alignmentObj.tSize, alignmentObj.tBeg, alignmentObj.tEnd, alignmentObj.blockCount

            # Reverse complement the alignment
            alignmentObj.rev_complement()

            #print "** BLAT-AFTER: ", alignmentObj.qName, alignmentObj.qSize, alignmentObj.qBeg, alignmentObj.qEnd, alignmentObj.tName, alignmentObj.tSize, alignmentObj.tBeg, alignmentObj.tEnd, alignmentObj.blockCount

            # Incorporate alignment to the list
            alignmentListRev.append(alignmentObj)

        return alignmentListRev 
        


class blat_alignment():
    """
    Blat alignment class.

    Methods:
    - in_target_region
    - partial_alignment
    """

    def __init__(self, alignment):
        """
            Initialize blat alignment object.

            Input:
            1) alignment. blat alignment in psl format

            Output:
            - Blat aligment object variables initialized
        """
        alignment = alignment.rstrip('\n')
        alignment = alignment.split("\t")

        # Define blat alignment variables
        self.matches = int(alignment[0])
        self.misMatches = int(alignment[1])
        self.repMatches = int(alignment[2])
        self.nCount = int(alignment[3])
        self.qNumInsert = int(alignment[4])
        self.qBaseInsert = int(alignment[5])
        self.tNumInsert = int(alignment[6])
        self.tBaseInsert = int(alignment[7])
        self.strand = alignment[8]
        self.qName = alignment[9]
        self.qSize = int(alignment[10])
        self.qBeg = int(alignment[11])
        self.qEnd = int(alignment[12])
        self.tName = alignment[13]
        self.tSize = int(alignment[14])
        self.tBeg = int(alignment[15])
        self.tEnd = int(alignment[16])
        self.blockCount = int(alignment[17])
        self.blockSizes = alignment[18]
        self.qStarts = alignment[19]
        self.tStarts = alignment[20]

        # Other
        self.alignType = ""

    def perc_query_covered(self):
        """

        """
        percQueryCovered = float(self.qEnd - self.qBeg) / self.qSize * 100  
        return percQueryCovered


    def rev_complement(self):
        """
            Make the reverse complementary aligment. This would be the alignment information of the reverse complementary original sequence

            Output:
            - Update alignment information for reverse complementary alignment. Need to be improved to also update tStarts variable
        """
        ## Reverse complement strand
        switchStrandDict = {'+': '-', '-': '+'}
        self.strand = switchStrandDict[self.strand]

        ## Reverse complement query start and end positions
        updatedBeg = self.qSize - self.qEnd
        updatedEnd = self.qSize - self.qBeg

        self.qBeg = updatedBeg
        self.qEnd = updatedEnd


class contig():
    """
    Methods:
    """

    def __init__(self, contigId, contigSeq):
        """
            Initialize contig object.

            Input:
            - 

            Output:
            - Contig object variables initialized
        """
        # Contig id and sequence
        self.ID = contigId
        self.seq = contigSeq
        self.length = int(len(self.seq))

        # Contig alignment information
        self.alignmentObjList = []   # List of blat alignments for this contig
        self.representativeAlignmentObj = "NA"

        # bkp information
        self.informative = ""        
        self.bkpDict = {}

        # Number supporting reads
        self.nbReads = ""     
    
        # Clipping side   
        self.clippedSide = ""

    def rev_complement(self, seq):
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

    def best_single_alignment(self):
        """

        """
        ## Sort alignment based on percentage of query covered
        alignmentObjListSorted = sorted(self.alignmentObjList, key=lambda alignmentObj: alignmentObj.perc_query_covered(), reverse=True)        

        ## Pick alignment covering the highest percentage of query sequence
        bestAlignmentObj = alignmentObjListSorted[0]

        return bestAlignmentObj

    def pair_complementary_alignments(self):
        """

        """
        ### Sort alignments based on contig beginning alignment coordinates
        alignmentObjListSorted = sorted(self.alignmentObjList, key=lambda alignmentObj: alignmentObj.qBeg, reverse=False)
        #print "sorted-alignments: ", alignmentObjListSorted, [alignmentObj.qBeg for alignmentObj in alignmentObjListSorted]

        ### Make pairs of overlapping partial alignments
        pairedAlignmentObjList = []
       
        ## Generate tuple list with all the possible alignment pairs
        alignmentObjPairsList = list(itertools.combinations(alignmentObjListSorted, 2))
         
        ## Iterate over all the possible pairs of contig's blat alignments
        for alignmentObjA, alignmentObjB in alignmentObjPairsList:
                
            ## Check if both alignments overlap:
            status, nbBases = overlap(alignmentObjA.qBeg, alignmentObjA.qEnd, alignmentObjB.qBeg, alignmentObjB.qEnd, 15)

            print status, nbBases

            ## Partial overlap
            if (status == 'partial'):    
        

                pairedAlignmentObj = chained_alignment([alignmentObjA, alignmentObjB])     
                percQueryCoveredA = alignmentObjA.perc_query_covered()
                percQueryCoveredB = alignmentObjB.perc_query_covered()
                percQueryCoveredPaired = pairedAlignmentObj.perc_query_covered()
                print "TEST: ", percQueryCoveredA, percQueryCoveredB, percQueryCoveredPaired

                ## Paired alignments covers more % of contig sequence than individual alignments:
                if (percQueryCoveredPaired > percQueryCoveredA) and (percQueryCoveredPaired > percQueryCoveredB):
                        
                    ## Add paired alignment to the list
                    pairedAlignmentObjList.append(pairedAlignmentObj)

        ## A) No paired alignment
        if (len(pairedAlignmentObjList) == 0):
            bestPairedAlignmentObj = "NA"

        ## B) At least one paired alignment
        else:
            ## Sort alignments based on percentage of query covered
            pairedAlignmentObjSortedList = sorted(pairedAlignmentObjList, key=lambda pairedAlignmentObj: pairedAlignmentObj.perc_query_covered(), reverse=True)        
        
            ## If several posible paired alignments, pick the one covering the highest percentage of query sequence
            bestPairedAlignmentObj = pairedAlignmentObjSortedList[0]
    
        return bestPairedAlignmentObj, pairedAlignmentObjList


    def chain_alignments(self, chainedAlignmentObjList):
        """

        """
        chainedAlignmentObjNewList = []

        ## Generate tuple list with all the possible chained alignment pairs
        chainedAlignmentObjPairsList = list(itertools.combinations(chainedAlignmentObjList, 2))

        ## Iterate over all the possible pairs of chained alignments
        for chainedAlignmentObjA, chainedAlignmentObjB in chainedAlignmentObjPairsList:

            alignmentObjList = chainedAlignmentObjA.alignmentList + chainedAlignmentObjB.alignmentList
                    
            ## There is some alignment in common between both chained alignments 
            if (len(alignmentObjList) != len(set(alignmentObjList))):
          
                ## Make not redundant list of alignments
                alignmentObjFilteredList = list(set(alignmentObjList))

                ## Sort alignmnents based on query beginning position
                alignmentObjSortedList = sorted(alignmentObjFilteredList, key=lambda alignmentObj: alignmentObj.qBeg, reverse=False)

                ## Create chained alignment object
                chainedAlignmentObj = chained_alignment(alignmentObjSortedList)

                ## Maximum percentage of query sequence for each paired chained alignment
                maxPerc = max([chainedAlignmentObjA.perc_query_covered(), chainedAlignmentObjB.perc_query_covered()])

                ## Chained alignment better than the best paired chained alignment
                if (chainedAlignmentObj.perc_query_covered() > maxPerc):

                    chainedAlignmentObjNewList.append(chainedAlignmentObj)
                    #print "**Chain_alignmnets: ", chainedAlignmentObj.alignmentList, len(chainedAlignmentObj.alignmentList), chainedAlignmentObj.perc_query_covered()

        nbAlignments = len(chainedAlignmentObjNewList)

        ## A) No chained alignment
        if (nbAlignments == 0):
            bestChainedAlignmentObj = "NA"
            return bestChainedAlignmentObj 

        else:
    
            ## B) Only one chained alignment (OR best one covers 100% of contig sequence) -> to implement
            if (nbAlignments == 1):
                bestChainedAlignmentObj = chainedAlignmentObjNewList[0]
                return bestChainedAlignmentObj 

            ## C) Multiple chained alignments
            else:

                ## Provisional, select the best one (In the future it should make recursion and attemp to extend again the chain...)
                chainedAlignmentObjNewListSorted = sorted(chainedAlignmentObjNewList, key=lambda chainedAlignmentObj: chainedAlignmentObj.perc_query_covered(), reverse=True)        

                bestChainedAlignmentObj = chainedAlignmentObjNewListSorted[0]
                return bestChainedAlignmentObj

     
    def representative_alignment(self):
        """
            
            
            Output: representativeAlignmentObj
        """        
        nbAlignments = len(self.alignmentObjList)

        ## A) Contig without blat hits
        if (nbAlignments == 0):
                  
            #print "A) Contig without blat hits"  
            self.representativeAlignmentObj = chained_alignment([])

        ## Contig with at least one blat hit
        else:

            ## 1. Identify single alignment spanning the largest piece of contig sequence
            bestAlignmentObj = self.best_single_alignment()            

            #print "bestAlignmentObj: ", bestAlignmentObj, bestAlignmentObj.perc_query_covered()

            ## B) Only one alignment or multiple but one spans the full contig sequence
            if (nbAlignments == 1) or (bestAlignmentObj.perc_query_covered() == 100):
                
                #print "B) Only one alignment or multiple but one spans the full contig sequence"
                self.representativeAlignmentObj = chained_alignment([bestAlignmentObj])
                           
            ## Multiple alignments and the best one does not span the full contig sequence 
            else:
        
                ## 2. Search for pairs of complementary alignments
                bestPairedAlignmentObj, pairedAlignmentObjList = self.pair_complementary_alignments()

                ## C) No paired complementary alignments OR best single alignment better/equal than best paired alignment 
                if (bestPairedAlignmentObj == "NA"):
                
                    #print "C) No paired complementary alignments OR best single alignment better/equal than best paired alignment"
                    self.representativeAlignmentObj = chained_alignment([bestAlignmentObj])

                ## D) Only one paired alignment or multiple but one spans the full contig sequence
                elif (len(pairedAlignmentObjList) == 1) or (bestPairedAlignmentObj.perc_query_covered() == 100):
    
                    #print "D) Only one paired alignment or multiple but one spans the full contig sequence" 
                    self.representativeAlignmentObj = bestPairedAlignmentObj
                    #print "bestPairedAlignmentObj: ", bestPairedAlignmentObj, bestPairedAlignmentObj.perc_query_covered()

                ## Multiple paired alignments and the best one does not span the full contig sequence 
                else:
                    chainedAlignmentObj = self.chain_alignments(pairedAlignmentObjList)
                    
                    # E) No chained complementary alignments OR best paired alignment better/equal than best paired alignment
                    if (chainedAlignmentObj == "NA"):
                        #print "E) No chained complementary alignments OR best paired alignment better/equal than best paired alignment"
                        self.representativeAlignmentObj = bestPairedAlignmentObj
        
                    # F) Chained complementary alignment better than paired alignment
                    else:
                        #print "F) Chained complementary alignment better than paired alignment"
                        self.representativeAlignmentObj = chainedAlignmentObj
   

    def fix_contig_orientation(self, chrom):
        """
         
        """     
        alignmentList = self.representativeAlignmentObj.alignmentList
    
        # Make not-redundant list with the strand for the alignments in the target region
        strandList = list(set([alignmentObj.strand for alignmentObj in alignmentList if alignmentObj.tName == chrom]))

        # A) No alignment in the target region
        if len(strandList) == 0:
            print "A) No alignment in the target region"
    
        # B) Ambiguous strand as multiple possibilities (+ and -) 
        elif len(strandList) > 1:
            print "B) Ambiguous strand as multiple possibilities (+ and -)"

        # C) Single possible strand
        else:
            print "C) Single possible strand"
        
            ## Minus strand
            if (strandList[0] == "-"):
    
                ### Make reverse complementary
                ## Contig sequence

                #print "CONTIG-BEFORE: ", self.seq
                self.seq = self.rev_complement(self.seq)
                #print "CONTIG-AFTER: ", self.seq

                ## Blat alignments
                #print "BLAT-BEFORE: ", self.representativeAlignmentObj.alignmentList
                self.representativeAlignmentObj.alignmentList = self.representativeAlignmentObj.rev_complement()
                #print "BLAT-AFTER: ", self.representativeAlignmentObj.alignmentList


    def is_complete_polyA(self, targetSeq):
        """        
        """
        ## Remove N characters
        targetSeqTrimmed = targetSeq.replace("N", "")

        targetSeqLen = len(targetSeqTrimmed)

        # Compute percentage of A and T in the clipped sequence
        if (targetSeqLen == 0):
            percA = 0
            percT = 0
        else:
            nbA = targetSeqTrimmed.count("A")
            nbG = targetSeqTrimmed.count("G")
            nbC = targetSeqTrimmed.count("C")
            nbT = targetSeqTrimmed.count("T")

            percA = float(nbA) / (nbA + nbG + nbC + nbT) * 100
            percT = float(nbT) / (nbA + nbG + nbC + nbT) * 100

        # Classify the input sequence as poly-A if at least 4 bases long and >= 80% bases are A or T
        if (targetSeqLen >= 4) and ((percA >= 80) or (percT >= 80)):
            polyASeq = targetSeq
        else:
            polyASeq = "NA"

        return polyASeq



    def is_partial_polyA(self, targetSeq, side):
        """        
        """
        windowSize = 10
        targetSeqLen = len(targetSeq)
        
        print "is_partial_polyA: ", targetSeq, targetSeqLen, side

        ## Search for polyA tails in the target seq. taking consecutive slices of X bases (X = windowSize)
        polyASeq = ""

        for pos in range(0, targetSeqLen, windowSize):
 
            # Parse sequence forward
            if (side == "beg"):
                # Take slice
                beg = pos
                end = beg + windowSize
            
                if (end > targetSeqLen):
                    end = targetSeqLen
    
            # Parse sequence backward
            else:
                # Take slice
                beg = targetSeqLen - pos - windowSize
                end = targetSeqLen - pos
            
                if (beg < 0):
                    beg = 0

            seq = targetSeq[beg:end]
            print "PAJAROO: ", beg, end, targetSeqLen, targetSeq, seq, side, polyASeq
 
            # Compute percentage of A and T in the given slice
            nbA = seq.count("A") 
            nbG = seq.count("G") 
            nbC = seq.count("C") 
            nbT = seq.count("T") 
            nbN = seq.count("N")
       
            percA = float(nbA) / (nbA + nbG + nbC + nbT + nbN) * 100
            percT = float(nbT) / (nbA + nbG + nbC + nbT + nbN) * 100
            
            # A) Classify the slide as poly-A if >= 80% are T or A
            if (percA >= 80) or (percT >= 80):
                
                # Add sequence slice to the polyA sequence
                polyASeq = polyASeq + seq 
            
            # B) Stop parsing sequence if the slide is not poly-A 
            else:
                    break

        if (polyASeq == ""):
            polyASeq = "NA"
        
        return polyASeq
            


    def is_polyA(self, targetSeq, side):
        """        
        """

        ## Check if complete polyA
        polyASeq = self.is_complete_polyA(targetSeq)

        # PolyA not found -> Check if partial polyA
        if (polyASeq == "NA"):
        
            polyASeq = self.is_partial_polyA(targetSeq, side)
        
        return polyASeq


    def is_polyA_informative(self, alignmentObj):
        """
        """
        ## Begin clipped
        #    #########-----------
        if (self.clippedSide  == "beg"):
            targetSeq = self.seq[:alignmentObj.qBeg].upper()
            polyASeq = self.is_polyA(targetSeq, 'end')

            if (polyASeq != "NA"):
                polyA = True
                self.bkpDict["polyA"] = polyASeq
                self.bkpDict["alignmentObjRT"] = "NA"
            else:
                polyA = False
                
        ## End clipped 
        #    ---------#########
        else:
            targetSeq = self.seq[alignmentObj.qEnd:].upper()
            polyASeq = self.is_polyA(targetSeq, 'beg')

            if (polyASeq != "NA"):
                polyA = True  
                self.bkpDict["polyA"] = polyASeq
                self.bkpDict["alignmentObjRT"] = "NA"
            else:
                polyA = False

        return polyA


    def is_single_informative(self, chrom):
        """
        """                
        alignmentObj = self.representativeAlignmentObj.alignmentList[0]

        # A) Contig does not align on the target region
        if (alignmentObj.tName != chrom):        
            print "CONTIG STATUS: Not informative, single alignment ouside the target region"
            informative = False

        ## B) Not informative, single alignment on the target region but spanning 100% of its sequence
        elif (self.representativeAlignmentObj.perc_query_covered() == 100):
            print "CONTIG STATUS: Not informative, single alignment on the target region but spanning 100% of its sequence"
            informative = False

        ## C) Single partial alignment -> Poly-A contig candidate
        else:
            informative = self.is_polyA_informative(alignmentObj)

        return informative
  

    def is_chained_informative(self, chrom):
        """

        """
        alignmentList = self.representativeAlignmentObj.alignmentList
        diffSeqList = list(set([ alignmentObj.tName for alignmentObj in alignmentList ]))
        nbDiffSeq = len(diffSeqList)
        
        ## A) Contig aligning multiple times in the same sequence (target region, transposon, exon or transduced region)
        if (nbDiffSeq == 1):

            print "A) Contig aligning multiple times in the same sequence (target region, transposon, exon or transduced region)"

            # a) Contig does not align on the target region
            if (diffSeqList[0] != chrom):        
                print "CONTIG STATUS: Not informative, chained alignment ouside the target region"
                informative = False

            # b) Not informative, single alignment on the target region but spanning 100% of its sequence
            elif (self.representativeAlignmentObj.perc_query_covered() == 100):
                print "CONTIG STATUS: Not informative, chained alignment on the target region but spanning 100% of its sequence"
                informative = False

            # c) Chained partial alignment on the target region -> Poly-A contig candidate   (TILL HERE)     
            else:
                print "Chained partial alignment on the target region -> Poly-A contig candidate"

                # a) End clipped
                # ------target_region------[....TE....]
                if (self.clippedSide  == "end"):
                    alignmentObj = alignmentList[0]  

                # b) Beg clipped
                # [....TE....]------target_region------
                else:
                    alignmentObj = alignmentList[len(alignmentList) - 1]

                informative = self.is_polyA_informative(alignmentObj)
            
        ## B) Contig aligning in different sequences -> Informative with or without poly-A
        else:
            print "B) Contig aligning in different sequences"    

            ### 1. Identify the alignment pair spanning the insertion breakpoint
            alignmentPairList = [] 

            ## Iterate over the aligments list
            for index, alignmentObjA in enumerate(alignmentList):

                ## Skip last alignment
                if (index + 1) != len(alignmentList):

                    ## Select next alignment
                    alignmentObjB = alignmentList[index + 1]

                    print "alignmentObjA: ", alignmentObjA.tName, alignmentObjA.qBeg, alignmentObjA.qEnd, alignmentObjA.tBeg, alignmentObjA.tEnd
                    print "alignmentObjB: ", alignmentObjB.tName, alignmentObjB.qBeg, alignmentObjB.qEnd, alignmentObjB.tBeg, alignmentObjB.tEnd

                    ## Current and next alignment in a different template. Possible combinations: 
                    # target_region -> L1
                    # target_region -> TD2 (I think I should include an alignment step for TD1 as well)
                    # target_region -> PSD
                    # all these combinations can be in the other way around
                    if alignmentObjA.tName != alignmentObjB.tName:
                            
                        # A) End clipped
                        # ------target_region------[....TE....]
                        if (self.clippedSide  == "end") and (alignmentObjA.tName != "L1") and (alignmentObjA.tName != "PSD"):
                            #print "A) Inserted sequence on the right"
                            informative = True
                            self.bkpDict["alignmentObjRT"] = alignmentObjB

                            ## Check for poly-A breakpoint
                            nbBases = alignmentObjB.qBeg - alignmentObjA.qEnd

                            if (nbBases > 0):
                                targetSeq = self.seq[alignmentObjA.qEnd : (alignmentObjA.qEnd + nbBases)].upper()
                                self.bkpDict["polyA"] = self.is_complete_polyA(targetSeq)
                            else:
                                targetSeq = self.seq[alignmentObjA.qEnd: ].upper()
                                self.bkpDict["polyA"] = self.is_complete_polyA(targetSeq)
                   
                            # Once breakpoint found stop iterating
                            break

                        # B) Beg clipped
                        # [....TE....]------target_region------
                        elif (self.clippedSide  == "beg") and (alignmentObjB.tName != "L1") and (alignmentObjB.tName != "PSD"):
                            #print "B) Inserted sequence on the left"
                            informative = True
                            self.bkpDict["alignmentObjRT"] = alignmentObjA

                            ## Check for poly-A breakpoint
                            nbBases = alignmentObjB.qBeg - alignmentObjA.qEnd

                            if (nbBases > 0):
                                targetSeq = self.seq[(alignmentObjB.qBeg - nbBases) : alignmentObjB.qBeg].upper()
                                self.bkpDict["polyA"] = self.is_complete_polyA(targetSeq)
                            else:
                                targetSeq = self.seq[: alignmentObjB.qBeg].upper()
                                self.bkpDict["polyA"] = self.is_complete_polyA(targetSeq)
                    
                            # Once breakpoint found stop iterating
                            break
                        
                        # C)
                        else:
                            informative = False

        return informative           
                
    def is_informative(self, tdType, coordinates):

        """  

        """
        chrom, beg, end =  coordinates.split("_")

        print "*** ANALYZING CONTIG: ", tdType, coordinates, chrom, beg, end, self.ID, self.seq, len(self.alignmentObjList), "***"

        ## 1. Obtain contig representative alignment (can be a single alignment or a set of complementary alignments)
        self.representative_alignment()
        
        ## 2. Reverse contig if the reverse complementary of the contig is aligning on the insertion region
        self.fix_contig_orientation(chrom)

        ## 3. Determine if informative contig based on the representative alignment
        nbAlignments = len(self.representativeAlignmentObj.alignmentList)
        
        ## A) Not informative -> Contig without blat hits 
        if (nbAlignments == 0):
            print "CONTIG STATUS: Not informative, no blat hits"
            self.informative = False

        ## B) Contig with single alignment
        elif (nbAlignments == 1):
            self.informative = self.is_single_informative(chrom)    

        ## C) Contig with chained alignment
        else:
            self.informative = self.is_chained_informative(chrom)

 
class insertion():
    """
    Transposable element insertion class.

    Methods:

    """

    def __init__(self, family, tdType, insertionCoord, contigsPath, blatPath, readPairsPlus, readPairsMinus, srcElement, transductionInfo, pseudogeneInfo, grInfo, sampleId):
        """
            Initialize insertion object.

            Input:
            1) family. TE family (L1, Alu, SVA or ERVK)
            2) tdType. Insertion type:  td0 (solo-insertion), td1 (partnered-transduccion), td2 (orphan-transduction) or psd (pseudogene insertion).
            3) coordinates. TraFiC insertion range.
            4) contigsPath. Fasta file containing the assembled contigs.
            5) blatPath. psl file containing the blat aligments for the assembled contigs.
            8) readPairsPlus. List of + cluster supporting reads.
            9) readPairsMinus. List of - cluster supporting reads.
            10) srcElement. Source element information in the format: cytobandId"_"srcElementType"_"chromSource"_"begSource"_"endSource"_"orientationSource. 
            11) transductionInfo. Transduction information in the format: chromSource"_"transductBeg"_"transductEnd"_"transductRnaLen"_"transductLen.
            12) pseudogeneInfo. Pseudogene information in the format: srcgene":"chrExonA"_"begExonA"-"endExonA":"chrExonB"_"begExonB"-"endExonB
            13) grInfo. Insertion associated genomic rearrangement (DEL, DUP, TRANS or NA) 
            14) sampleId. Sample identifier to be incorporated in the SL field of the output VCF. In PCAWG we use the normal_wgs_aliquot_id or tumor_wgs_aliquot_id.
            
            Output:
            - Insertion object variables initialized
        """
        self.family = family
        self.tdType = tdType
        self.coordinates = insertionCoord
        self.grInfo = grInfo
        self.sampleId = sampleId
        self.readPairIdsPlus = readPairsPlus
        self.nbReadPairsPlus =  len(readPairsPlus.split(','))
        self.readPairIdsMinus = readPairsMinus
        self.nbReadPairsMinus =  len(readPairsMinus.split(','))

        # Organize contigs into a dictionary
        self.contigsDict = self.read_contigs(contigsPath)

        # Add blat alignments to cluster object
        self.add_contig_alignments(blatPath)

        # Initialize insertion properties
        self.cipos = "NA"
        self.targetSiteLen = "NA"
        self.mechanism = "NA"
        self.orientation = "NA"
        self.elementLength = "NA"
        self.elementRange = "NA"
        self.structure = "NA"
        self.score = "NA"
        self.chrom = "NA"
        self.bkpAPos = "NA"
        self.bkpBPos = "NA"
        self.bkpAContigSeq = "NA"
        self.bkpBContigSeq = "NA"

        ## Set transduction or pseudogene specific fields
        # A) Solo insertion (TD0). Not applicable
        if (self.tdType == "TD0"):
            self.cytobandId = "NA"
            self.sourceElementType = "NA"
            self.sourceElementCoord = "NA"
            self.tdCoord = "NA"
            self.tdLen = "NA"
            self.tdLenRna = "NA"
            self.srcgene = "NA"

        # B) Partnered or orphan transduction (TD1 and TD2)
        elif (self.tdType in ["TD1", "TD2"]):

            ### L1 source element information
            srcElementList = srcElement.split("_")

            self.cytobandId = srcElementList[0]
            self.sourceElementType = srcElementList[1]
            self.sourceElementCoord = srcElementList[2] + "_" + srcElementList[3] + "_" + srcElementList[4]
            
            status = srcElementList[5]

            ### Mobilized region coord
            transductionInfoList = transductionInfo.split("_")
            self.tdCoord = transductionInfoList[0] + "_" + transductionInfoList[1] + "_" + transductionInfoList[2]
            self.srcgene = "NA"

            # a) Putative. Uncharacterized germline or somatic source element
            if (status == "putative"):
                self.tdLenRna = "NA"                
                self.tdLen = "NA"
              
            # b) Characterized germline source element
            else:
                self.tdLenRna = transductionInfoList[3]                
                self.tdLen = transductionInfoList[4]
                
        # C) pseudogene insertion (PSD)
        elif (self.tdType == "PSD"):
            self.cytobandId = "NA"
            self.sourceElementType = "NA"
            self.sourceElementCoord = "NA"
            self.tdLen = "NA"
            self.tdLenRna = "NA"

            # parse pseudogene info, format:
            # gene:chrA_begA-endA:chrB_begB-begB
            regex = r'(?P<gene>\w+):(?P<chrA>\w+)_(?P<begA>\d+)-(?P<endA>\d+):(?P<chrB>\w+)_(?P<begB>\d+)-(?P<endB>\d+)'

            m = re.search(regex, pseudogeneInfo)
            srcgene = m.group("gene")
            chrExonA = m.group("chrA")
            begExonA = int(m.group("begA"))
            endExonA = int(m.group("endA"))
            chrExonB = m.group("chrB")
            begExonB = int(m.group("begB"))
            endExonB = int(m.group("endB"))

            self.srcgene = srcgene

            # construct source region from exon coordinates
            minBeg = min(begExonA, begExonB)
            maxEnd = max(endExonA, endExonB)
            self.tdCoord = "%s_%d_%d" % (chrExonA, minBeg, maxEnd)

    def read_contigs(self, contigsFasta):
        """
            Read fasta file with the cluster's assembled contigs and produce a dictionary with
            the contig ids as keys and the corresponding contig objects as values.

            Input:
            1) contigsFasta. Fasta file containing the assembled contigs for the given
            cluster.

            Output:
            1) contigsDict. Dictionary containing the contig ids as keys and the corresponding
            contig objects as values.
        """
        fastaObj = formats.fasta(contigsFasta)

        contigsDict = {}

        ### For each contig create a contig object and add it to the dictionary
        # using the contig id as key
        for contigId in fastaObj.fastaDict:
            contigSeq = fastaObj.fastaDict[contigId]

            # Create contig object
            contigObj = contig(contigId, contigSeq)

            # breakpoint position (cluster_X_9823493_end_1)
            contigIdList = contigId.split('_')
            contigObj.bkpDict["chrom"] = contigIdList[1]
            contigObj.bkpDict["pos"] = int(contigIdList[2])

            # Clipping side
            contigObj.clippedSide = contigIdList[3]

            # Add number of supporting reads 
            contigObj.nbReads = int(contigIdList[4])

            # Add contig object to the dictionary
            contigsDict[contigId] = contigObj

        return contigsDict

    def read_alignments(self, blatPath):
        """
            Read a psl file containing the contig blat aligments on the reference genome, generate a blat alignment
            object per aligment and store all of them in a dictionary.

            Input:
            1) blatPath. Psl file containing blat aligments for the assembled contigs.

            Output:
            1) alignmentsDict. Dictionary containing the contig ids as keys and the list of alignment objects corresponding
                               to each contig as value
        """

        blat = open(blatPath, 'r')
        alignmentsDict = {}

        for alignment in blat:

            ## Create blat alignment object
            alignmentObject = blat_alignment(alignment)

            ## Initialize contig alignment list if it does not exists
            if alignmentObject.qName not in alignmentsDict:
                alignmentsDict[alignmentObject.qName] = []

            ## Add alignment object to the list
            alignmentsDict[alignmentObject.qName].append(alignmentObject)

        return alignmentsDict

    def add_contig_alignments(self, blatPath):
        """
            Read a psl file containing the contig blat alignments on the reference genome and TE sequences and associate
            each alignment to the corresponding contig object.

            Input:
            1) blatPath. Psl file containing blat aligments for the assembled contigs.

            Output:
            1) For each contig object in the dictionary sets the 'alignmentObjList' variable.
               This variable holds the list of alignment objects corresponding to the given contig.

        """
        alignmentsDict = self.read_alignments(blatPath)

        ## Add the alignment lists to their respective contig objects
        # Iterate over the alignments dictionary.
        for contigId in alignmentsDict:
            alignmentList = alignmentsDict[contigId]
            self.contigsDict[contigId].alignmentObjList = alignmentList


    def find_informative_contigs(self):
        """
        """

        informativeContigsClippedEndList = []
        informativeContigsClippedBegList = []

        info(str(len(self.contigsDict)) + ' input contigs')

        ## Check for each contig if it is informative or not. Make list of informative contigs
        for contigId in self.contigsDict:
            contigObj = self.contigsDict[contigId]

            # Determine if informative contig:
            contigObj.is_informative(self.tdType, self.coordinates)

            print "INFORMATIVE-CONTIG?: ", contigObj.informative, contigObj.bkpDict, contigObj.clippedSide

            # Add informative contig to the list
            if contigObj.informative:

                if (contigObj.clippedSide == "end"):
                    informativeContigsClippedEndList.append(contigObj)
                else:
                    informativeContigsClippedBegList.append(contigObj)
     

            print "------------------"

        return informativeContigsClippedEndList, informativeContigsClippedBegList


    def select_best_informative_contig(self, informativeContigsList):
        """
        """
        print "informativeContigsList: ", informativeContigsList

        # Select clusters with the highest number of supporting reads
        maxNbReads = max([ informativeContig.nbReads for informativeContig in informativeContigsList])
        bestInformativeContigsList = [informativeContig for informativeContig in informativeContigsList if informativeContig.nbReads == maxNbReads]    

        # a) Single cluster with the highest number of supporting reads 
        if (len(bestInformativeContigsList) == 1):
            print "maxNbReads: ", maxNbReads, len(bestInformativeContigsList), bestInformativeContigsList
            bestInformativeContigObj = bestInformativeContigsList[0]
        
        # b)  Multiple possible clusters... use another criteria... NA is provisional...
        else:
            bestInformativeContigObj = "NA"
 
        return bestInformativeContigObj
 

    def target_site(self, informativeContigClippedEndObj, informativeContigClippedBegObj):
        """
        """

        # a) Unknown target site 
        if (informativeContigClippedEndObj == "NA") or (informativeContigClippedBegObj == "NA"):
            targetSiteLen = "NA"

        # b) Known target site as both breakpoints characterized
        else:

            clippedEndPos = informativeContigClippedEndObj.bkpDict["pos"]
            clippedBegPos = informativeContigClippedBegObj.bkpDict["pos"]
            clippedEndAlignObj = informativeContigClippedEndObj.representativeAlignmentObj.alignmentList[0]
            clippedBegAlignObj = informativeContigClippedBegObj.representativeAlignmentObj.alignmentList[-1]   
        
            overlapping, nbBases = overlap(clippedEndAlignObj.tBeg, clippedEndAlignObj.tEnd, clippedBegAlignObj.tBeg, clippedBegAlignObj.tEnd, 0)

            ## a) Target Site Duplication
            # ----------contig----------bkp
            #             bkp-------contig-----------
            if overlapping != "not":
                targetSiteLen = abs(clippedEndPos - clippedBegPos) 
                print 'TS DUP', self.coordinates, targetSiteLen


            ## b) Target Site Deletion 
            # ----------contig----------bkp
            #                                       bkp-------contig-----------
            else:
                targetSiteLen = abs(clippedEndPos - clippedBegPos) * -1
                print 'TS DEL', self.coordinates, targetSiteLen

        return targetSiteLen


    def insertion_orientation(self, clippedEndPolyA, clippedBegPolyA):
        """
        """
        print "*** insertion_orientation ***"

        ## Clipped end RT and clipped beg poly-A
        # + orientation:
        #    Clipped-end informative contig      --------------####TE####
        #    Clipped-beg informative contig                              AAAAAAAAAA--------------
        if (clippedEndPolyA == "NA") and (clippedBegPolyA != "NA"):
            orientation = "+"

        ## Clipped end poly-A and clipped beg RT
        # - orientation 
        #    Clipped-end informative contig      --------------AAAAAAAAAA
        #    Clipped-beg informative contig                              ####TE####--------------
        else:
            orientation = "-"

        return orientation


    def element_length(self, alignmentObj):
        """
        """
        print "*** element_length ***"
            
        ## A) Solo integration or partnered transduction
        if (self.tdType == "TD0") or (self.tdType == "TD1"):
            elementLength = alignmentObj.tSize - alignmentObj.tBeg
            elementRange = str(alignmentObj.tBeg) + '-' + str(alignmentObj.tSize)

        ## B) Orphan transduction:
        else:
            elementLength = 0
            elementRange = "NA"

        return elementLength, elementRange


    def insertion_structure(self, alignmentObjElement, orientation, elementLength):
        """
        """
        print "*** insertion_structure ***"

        percLength = float(elementLength) / alignmentObjElement.tSize * 100

        ## A) Full length (>=98% consensus sequence):
        if (percLength >= 98):
            structure = "FULL"
            #print "FULL: ", self.tdType, percLength, elementLength, alignmentObjElement.tSize, alignmentObjElement.tName, orientation, alignmentObjElement.strand

        ## B) 5' inverted
        elif (orientation != alignmentObjElement.strand):
            structure = "INV"
            #print "INV: ", self.tdType, percLength, elementLength, alignmentObjElement.tSize, alignmentObjElement.tName, orientation, alignmentObjElement.strand

        ## C) 3' deleted
        else:
            structure = "DEL"
            #print "DEL: ", self.tdType, percLength, elementLength, alignmentObjElement.tSize, alignmentObjElement.tName, orientation, alignmentObjElement.strand
    
        return structure


    def insertion_orientation_noPolyA(self, alignmentObjClippedEnd, alignmentObjClippedBeg):
        """
        """
        print "*** insertion_orientation_noPolyA ***"

        ## a) + Orientation
        if (alignmentObjClippedEnd.strand == "+") and (alignmentObjClippedBeg.strand == "+"):
            orientation = "+"

        ## b) - Orientation
        elif (alignmentObjClippedEnd.strand == "-") and (alignmentObjClippedBeg.strand == "-"):
            orientation = "-"
        
        ## c) Inconsistent orientation
        else:
            orientation = "inconsistent"

        return orientation


    def element_length_noPolyA(self, alignmentObjClippedEnd, alignmentObjClippedBeg, orientation):
        """
        """
        print "*** element_length_noPolyA ***"

        ## a) <--Clipped-End--> <--Clipped-Beg-->
        #     ...#################TE#################...   
        # Ejem:
        # 10_58867421_58867517 TD0 L1 6019 5573 5649 + +
        # 10_58867421_58867517 TD0 L1 6019 5605 5682 + +  
        if (alignmentObjClippedEnd.tBeg < alignmentObjClippedBeg.tBeg):
           
            elementLength = alignmentObjClippedBeg.tEnd - alignmentObjClippedEnd.tBeg
            elementRange = str(alignmentObjClippedEnd.tBeg) + '-' + str(alignmentObjClippedBeg.tEnd)
           
        
        ## b) <--Clipped-Beg--> <--Clipped-End-->
        #     ...#################TE#################...     
        # Ejem:        
        # 1_57286341_57286435 TD0 L1 6019 5599 5678 - -
        # 1_57286341_57286435 TD0 L1 6019 5538 5614 - -
        else:
            elementLength = alignmentObjClippedEnd.tEnd - alignmentObjClippedBeg.tBeg
            elementRange = str(alignmentObjClippedBeg.tBeg) + '-' + str(alignmentObjClippedEnd.tEnd)

        return elementLength, elementRange


    def features(self, informativeContigClippedEndObj, informativeContigClippedBegObj):
        """
        """
        print "*** features_function ***"

        # A) At least one breakpoint not characterized
        if (informativeContigClippedEndObj == "NA") or (informativeContigClippedBegObj == "NA"):
            
            mechanism = "NA"
            orientation = "NA"
            elementLength = "NA"
            elementRange = "NA"
            structure = "NA"

        # B) Both breakpoints characterized
        else:
            
            # a) Two poly-A breakpoints
            if (informativeContigClippedEndObj.bkpDict["polyA"] != "NA") and (informativeContigClippedBegObj.bkpDict["polyA"] != "NA"):

                mechanism = "DPA"
                orientation = "NA"
                elementLength = "NA"
                elementRange = "NA"
                structure = "NA"

            # b) Two RT breakpoints
            elif (informativeContigClippedEndObj.bkpDict["polyA"] == "NA") and (informativeContigClippedBegObj.bkpDict["polyA"] == "NA"):

                mechanism = "EI"
                structure = "NA"

                ## 1. Orientation (+ or -)
                orientation = self.insertion_orientation_noPolyA(informativeContigClippedEndObj.bkpDict["alignmentObjRT"], informativeContigClippedBegObj.bkpDict["alignmentObjRT"])

                ## 2. Integrated element length
                elementLength, elementRange = self.element_length_noPolyA(informativeContigClippedEndObj.bkpDict["alignmentObjRT"], informativeContigClippedBegObj.bkpDict["alignmentObjRT"], orientation)

            # c) One RT and one poly-A breakpoint 
            else:
                mechanism = "TPRT"
                elementRange = "NA"

                ## 1. Orientation (+ or -)
                orientation = self.insertion_orientation(informativeContigClippedEndObj.bkpDict["polyA"], informativeContigClippedBegObj.bkpDict["polyA"])
                alignmentObjElement = informativeContigClippedEndObj.bkpDict["alignmentObjRT"] if (orientation == "+") else informativeContigClippedBegObj.bkpDict["alignmentObjRT"]

                ## 2. Integrated element length
                elementLength, elementRange = self.element_length(alignmentObjElement)

                ## 3. Structure (FULL, DEL or INV)
                structure = self.insertion_structure(alignmentObjElement, orientation, elementLength)
        
        return mechanism, orientation, elementLength, elementRange, structure
 

    def insertion_score(self, informativeContigClippedEndObj, informativeContigClippedBegObj):
        """
        """
        print "*** insertion_score_function ***"

        ## a) Inconsistent features (score 1)
        if ('inconsistent' in [self.targetSiteLen, self.orientation, self.elementLength, self.structure]):
            score = '1'

        ## b) No breakpoint identified (score 2)
        elif (informativeContigClippedEndObj == "NA") and (informativeContigClippedBegObj == "NA"):
            score = '2'

        ## c) Both breakpoints identified (score 5)
        elif (informativeContigClippedEndObj != "NA") and (informativeContigClippedBegObj != "NA"):
            score = '5'        

        ## d) PolyA breakpoint identified (score 4)
        elif ((informativeContigClippedEndObj != "NA") and (informativeContigClippedEndObj.bkpDict["polyA"] != "NA")) or ((informativeContigClippedBegObj != "NA") and (informativeContigClippedBegObj.bkpDict["polyA"] != "NA")):
            score = '4'

        ## e) Element breakpoint identified (score 3)
        else:
            score = '3'
 
        return score

    def imprecise_bkp(self):
        """
            Compute confidence interval for imprecise breakpoints from + and - TraFiC cluster coordinates.

            Imprecise breakpoints are those that do not have any
            any informative contig associated

              + cluster                      - cluster
            -------------->              <---------------
                          <-----Mean----->
                          d/2    d      d/2

            Input:
            1) insertionCoord. TraFiC insertion coordinates. Format: ${chrom}_${beg}_${end}.
                               Example: 10_108820680_108820678.

            Output:
            1) breakpoint. Three elements list (chromosome, mean_position and confidence interval)
        """

        insertionCoordList = self.coordinates.split("_")

        chrom = insertionCoordList[0]
        beg = insertionCoordList[1]
        end = insertionCoordList[2]

        ## a) Paired clusters
        if (beg != "NA") and (end != "NA"):

            beg = int(beg)
            end = int(end)
            bkpAPos = int((beg + end)/2)
            dist = abs(end - beg)
            cipos = int(dist/2)

        ## b) Independent cluster:
        else:
            bkpAPos = int(beg)
            cipos = "NA"


        return chrom, bkpAPos, cipos

    def breakpoints(self, informativeContigClippedEndObj, informativeContigClippedBegObj):
        """
        """
        print "*** breakpoints_function ***"

        # a) Both breakpoints identified (inconsistent as well)
        if (self.score == '5') or (self.score == '1'):
            cipos = 0
            chrom = informativeContigClippedEndObj.bkpDict["chrom"]

            if (informativeContigClippedEndObj.bkpDict["pos"] < informativeContigClippedBegObj.bkpDict["pos"]):

                bkpAPos = informativeContigClippedEndObj.bkpDict["pos"] 
                bkpAContigSeq = informativeContigClippedEndObj.seq

                bkpBPos = informativeContigClippedBegObj.bkpDict["pos"]
                bkpBContigSeq = informativeContigClippedBegObj.seq
                
            else:

                bkpAPos = informativeContigClippedBegObj.bkpDict["pos"] 
                bkpAContigSeq = informativeContigClippedBegObj.seq

                bkpBPos = informativeContigClippedEndObj.bkpDict["pos"]
                bkpBContigSeq = informativeContigClippedEndObj.seq


        # b) One breakpoint identified
        elif (self.score == '4') or (self.score == '3'):
            cipos = 0

            if (informativeContigClippedEndObj != "NA"):
                chrom = informativeContigClippedEndObj.bkpDict["chrom"]
                bkpAPos = informativeContigClippedEndObj.bkpDict["pos"] 
                bkpAContigSeq = informativeContigClippedEndObj.seq

                bkpBPos = "NA"
                bkpBContigSeq = "NA"
    
            else:
                chrom = informativeContigClippedBegObj.bkpDict["chrom"]
                bkpAPos = informativeContigClippedBegObj.bkpDict["pos"] 
                bkpAContigSeq = informativeContigClippedBegObj.seq

                bkpBPos = "NA"
                bkpBContigSeq = "NA"               

        # c) No breakpoint identified (Imprecise breakpoint)      
        else:
            chrom, bkpAPos, cipos = self.imprecise_bkp()
            bkpAContigSeq = "NA"
            bkpBPos = "NA"
            bkpBContigSeq = "NA"

        return chrom, bkpAPos, bkpBPos, bkpAContigSeq, bkpBContigSeq, cipos    
                 
    def find_insertionBkp(self):
        """
           
        """
        #### 1. Search for informative contigs
        informativeContigsClippedEndList, informativeContigsClippedBegList = self.find_informative_contigs()
        print "NUMBER-INFORMATIVE-CONTIGS: ", len(informativeContigsClippedEndList), len(informativeContigsClippedBegList)

        nbInformativeEnd = len(informativeContigsClippedEndList)
        nbInformativeBeg = len(informativeContigsClippedBegList)

        ## Check if any informative cluster with reads clipped at the end
        # a) Any informative
        if (nbInformativeEnd == 0):
            informativeContigClippedEndObj = "NA"

        # b) One informative
        elif (nbInformativeEnd == 1):
            informativeContigClippedEndObj = informativeContigsClippedEndList[0]

        # c) Multiple informative -> Select the best one
        else:   
            informativeContigClippedEndObj = self.select_best_informative_contig(informativeContigsClippedEndList)

        ## Check if any informative cluster with reads clipped at the begin
        # a) Any informative
        if (nbInformativeBeg == 0): 
            informativeContigClippedBegObj = "NA"

        # b) One informative
        elif (nbInformativeBeg == 1):
            informativeContigClippedBegObj = informativeContigsClippedBegList[0]
        
        # c) Multiple informative -> Select the best one
        else:
            informativeContigClippedBegObj = self.select_best_informative_contig(informativeContigsClippedBegList)

        #### 2. Define insertion properties
        ## 2.1 Target site duplication or deletion
        self.targetSiteLen = self.target_site(informativeContigClippedEndObj, informativeContigClippedBegObj)

        ## 2.2 Insertion features: mechanism (TPRT, EI, 2POLYA), orientation, structure and length   
        self.mechanism, self.orientation, self.elementLength, self.elementRange, self.structure = self.features(informativeContigClippedEndObj, informativeContigClippedBegObj)

        ## 2.3 Insertion score
        self.score = self.insertion_score(informativeContigClippedEndObj, informativeContigClippedBegObj)

        ## 2.4 Breakpoints
        self.chrom, self.bkpAPos, self.bkpBPos, self.bkpAContigSeq, self.bkpBContigSeq, self.cipos = self.breakpoints(informativeContigClippedEndObj, informativeContigClippedBegObj)

        print "FEATURES: ", self.tdType, self.coordinates, self.targetSiteLen, self.mechanism, self.orientation, self.elementLength, self.elementRange, self.structure, self.score, self.bkpAPos, self.bkpBPos, self.bkpAContigSeq, self.bkpBContigSeq, self.cipos


    def convert2VCFline(self, genomeObj):
        """
        """  
        ## Define variables with VCF line fields
        CHROM = self.chrom
        POS = self.bkpAPos
        ID = "."
        REF = genomeObj.fastaDict[CHROM][POS - 1]  # Substract 1 since python string coordinates start in 0 while bkp position in 1.
        ALT = "<MEI>"
        QUAL = "."
        FILTER = "."
        INFO = "SVTYPE=<MEI>;CLASS=" + self.family + ";TYPE=" + self.tdType + ";MECHANISM=" + self.mechanism + ";SCORE=" + str(self.score) + ";BKPB=" + str(self.bkpBPos) + ";CIPOS=" + str(self.cipos) + ";STRAND=" + self.orientation + ";STRUCT=" + self.structure + ";LEN=" + str(self.elementLength) + ";RANGE=" + str(self.elementRange) + ";TSLEN=" + str(self.targetSiteLen) + ";SRCID=" + self.cytobandId + ";SRCTYPE=" + self.sourceElementType + ";SRC=" + self.sourceElementCoord + ";TDC=" + self.tdCoord + ";TDLEN=" + self.tdLen + ";TDLENR=" + self.tdLenRna + ";SRCGENE=" + self.srcgene + ";GR=" + self.grInfo + ";CONTIGA=" + self.bkpAContigSeq + ";CONTIGB=" + self.bkpBContigSeq + ";RP=" + self.readPairIdsPlus + ";RN=" + self.readPairIdsMinus    
        FORMAT = "RCP:RCN:SL"
        RP = str(self.nbReadPairsPlus)
        RN = str(self.nbReadPairsMinus)
        GENOTYPE = RP + ":" + RN + ":" + self.sampleId

        ## Generate VCF line object
        VCFlineList = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, GENOTYPE]
        VCFlineObj = formats.VCFline(VCFlineList)
        
        ## Generate info string again without NAs
        VCFlineObj.info = VCFlineObj.make_info()

        return VCFlineObj


#### MAIN ####

if __name__ == "__main__":
    ## Get user's input ##
    args = parse_args()
    inputPaths = args.inputPaths
    sampleId = args.sampleId
    fileName = args.fileName
    genome = args.genome
    outDir = args.outDir

    scriptName = os.path.basename(sys.argv[0])

    ## Display configuration to standard output ##
    print
    print "***** ", scriptName, " configuration *****"
    print "paths2bkpAnalysis: ", inputPaths
    print "sampleId: ", sampleId
    print "fileName: ", fileName
    print "genome: ", genome
    print "outDir: ", outDir
    print
    print "***** Executing ", scriptName, " *****"
    print

    ## Start ## 

    ## 0. Create reference genome fasta object
    header("Creating reference genome fasta object")
    genomeObj = formats.fasta(genome)

    ## 1. Create VCF object and print VCF header
    header("Creating VCF object and printing VCF header into the output file")
    VCFObj = formats.VCF()
    VCFObj.create_header(genome, genomeObj)

    ## 2. Per each insertion perform breakpoint analysis
    inputFile = open(inputPaths, 'r')

    # Analyze one insertion per iteration
    for line in inputFile:
        line = line.rstrip('\n')
        line = line.split("\t")

        # Get MEI info and files
        insertionInfo = line[0]
        family, tdType, insertionCoord = insertionInfo.split(":")
        contigsPath = line[1]
        blatPath = line[2]
        readPairsPlus = line[3]
        readPairsMinus = line[4]
        srcElement = line[5]
        transductionInfo = line[6]
        pseudogeneInfo = line[7]
        grInfo = line[8]            # genomic rearrangement info

        # Perform breakpoint analysis for the TE insertion
        # A) All the input files exist. 
        if (contigsPath != "NA" and os.path.isfile(contigsPath)) and (blatPath != "NA" or os.path.isfile(blatPath)):

            header("Tranposable Element Insertion Breakpoint Analysis (TEIBA) for: " + insertionInfo)

            ## Create insertion object and identify breakpoints from assembled contigs
            insertionObj = insertion(family, tdType, insertionCoord, contigsPath, blatPath, readPairsPlus, readPairsMinus, srcElement, transductionInfo, pseudogeneInfo, grInfo, sampleId)

            ## Characterize insertion breakpoints
            insertionObj.find_insertionBkp()

            ## Create VCFline object
            VCFlineObj = insertionObj.convert2VCFline(genomeObj)

            ## Add VCFline to the list in VCF object
            VCFObj.addLine(VCFlineObj)

        # B) Some input file do not exist
        else:
            message = "Input files for " + insertionCoord + " insertion do not exist"
            log("ERROR", message)

    ## 3. Sort MEI
    VCFObj.sort()

    ## 4. Write output VCF
    outFilePath = outDir + '/' + fileName + '.vcf'
    VCFObj.write_header(outFilePath)
    VCFObj.write_variants(outFilePath)

    ## Finish ##
    print
    print "***** Finished! *****"
    print


