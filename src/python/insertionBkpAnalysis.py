#!/usr/bin/env python
#coding: utf-8 


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


#### CLASSES ####

class insertion():
    """ 
    Transposable element insertion class. 
    
    A cluster can be + or - and it has associated the contigs resulting from the assembly of TE insertion supporting
    reads identified by TraFiC.  
    
    Methods:
    - create_cluster
    - find_insertionBkp
    - insertion_orientation
    - polyA
    - find_TSD
    - imprecise_bkp
    """

    def __init__(self, family, coordinates, contigsPlusPath, blatPlusPath, contigsMinusPath, blatMinusPath):
        """ 
            Initialize insertion object.
            
            Input:
            1) family. TE family (L1, Alu or SVA)
            2) coordinates. 
            3) contigsPlusPath. Fasta file containing the assembled contigs for the positive cluster. 
            4) blatPlusPath. psl file containing the blat aligments for the positive cluster's assembled contigs.
            5) contigsMinusPath. Fasta file containing the assembled contigs for the negative cluster. 
            6) blatMinusPath. psl file containing the blat aligments for the negative cluster's assembled contigs.
            
            Output:
            - Insertion object variables initialized
        """
        self.category = category
        self.coordinates = coordinates
        self.clusterPlusObj = self.create_cluster("+", contigsPlusPath, blatPlusPath)
        self.clusterMinusObj = self.create_cluster("-", contigsMinusPath, blatMinusPath)
        
    #### FUNCTIONS ####
    def create_cluster(self, ID, contigsPath, blatPath):
        """ 
            Create cluster object.
            
            Input:
            1) ID. Cluster id (+ or -)
            2) contigsPath. Fasta file containing the assembled contigs for the given 
                            cluster.
            3) blatPath. psl file containing the blat aligments for the assembled contigs.
            
            Output:
            1) clusterObj. Cluster object 
        """
        
        # Create cluster object
        clusterObj = cluster(ID, contigsPath)
        
        # Add blat alignments to cluster object
        clusterObj.add_alignments(blatPath)

        return clusterObj


    def find_insertionBkp(self, outDir):
        """ 
            Identify TE insertion breakpoints, TSD, orientation and poly-A sequence from contigs belonging to +
            and - clusters. 
            
            Input:
            1) outDir. Output directory (provisional).
            
            Output:
            1) score. Integer from 1 to 4:
                1. Informative contig (5' and 3') for both clusters (+ and -)
                2. Informative contig (5' or 3') for only one cluster 1 (+ or -)
                3. Any informative contig identified.
                4. Inconsistent insertion. Informative contig for both clusters of the same type (5',5' or 3',3'). 
            2) breakpoint. Breakpoint coordinates (tuple: chrom and pos)
            3) TSD. Target site duplication (tuple: TSD size and sequence). 
            4) orientation. TE insertion DNA strand/orientation (+ or -) 
            5) polyA. Poly-A sequence. 
        """
        
        ## Find informative contig for + cluster 
        subHeader("Searching informative contigs for + cluster")
        
        informativeContigPlusObj = self.clusterPlusObj.find_informative_contig(self.coordinates)

        ## Find informative contig for - cluster
        subHeader("Searching informative contigs for - cluster")
        
        informativeContigMinusObj = insertionObj.clusterMinusObj.find_informative_contig(self.coordinates)
        
        ### Determine insertion breakpoints, TSD and TE orientation from informative contigs 
        
        subHeader("Determining insertion breakpoint, TSD and TE orientation from informative contigs")
        
        ## Set variables 
        
        traficId = self.category + ":" + self.coordinates
        
        # A) There is an informative contig for both + and - clusters
        if (informativeContigPlusObj != "") and (informativeContigMinusObj != ""):
            
            info('+ and - informative clusters:') 
                
            score = 1
            bkpPlus = informativeContigPlusObj.informative[1]
            bkpMinus = informativeContigMinusObj.informative[1]
            typePlus = informativeContigPlusObj.informative[0]
            infoPlus = informativeContigPlusObj.informative[2]
            typeMinus = informativeContigMinusObj.informative[0]
            infoMinus = informativeContigMinusObj.informative[2]
            contigPlus = informativeContigPlusObj.seq
            contigMinus = informativeContigMinusObj.seq
            
            # Find Target Size Duplication (TSD)
            TSDlength, TSDseq = self.find_TSD(informativeContigPlusObj, informativeContigMinusObj)
                        
        # B) There is an informative contig for + cluster
        elif (informativeContigPlusObj != ""):
           
            info('+ informative cluster:')
            
            score = 2
            bkpPlus = informativeContigPlusObj.informative[1]
            bkpMinus = [ "na", "na" ]
            typePlus = informativeContigPlusObj.informative[0]
            infoPlus = informativeContigPlusObj.informative[2]
            typeMinus = "none"
            infoMinus = "none"
            TSDlength = "na"
            TSDseq = "na"
            contigPlus = informativeContigPlusObj.seq
            contigMinus = "na"
            
        # C) There is an informative contig for - cluster
        elif (informativeContigMinusObj != ""):
            
            info('- informative cluster:')
            
            score = 2
            bkpPlus = [ "na", "na" ]
            bkpMinus = informativeContigMinusObj.informative[1]
            typePlus = "none"
            infoPlus = "none"
            typeMinus = informativeContigMinusObj.informative[0]
            infoMinus = informativeContigMinusObj.informative[2]
            TSDlength = "na"
            TSDseq = "na"
            contigPlus = "na"
            contigMinus = informativeContigMinusObj.seq
        
        # D) There are not informative contigs for any of the clusters
        else:
            
            info('no informative clusters:')
            
            score = 3
            bkpPlus = [ "na", "na" ]
            bkpMinus = [ "na", "na" ]
            typePlus = "none"
            infoPlus = "none"
            typeMinus = "none"
            infoMinus = "none"
            TSDlength = "na"
            TSDseq = "na"
            contigPlus = "na"
            contigMinus = "na"
        
        # TE insertion orientation
        orientation = self.insertion_orientation(typePlus, typeMinus)
        
        # TE insertion structure
        structure, length, percLength = self.insertion_structure(typePlus, infoPlus, typeMinus, infoMinus, orientation)
        
        # Poly-A sequence
        polyA = self.polyA(typePlus, infoPlus, typeMinus, infoMinus)     
        
        # Modify score from 1 to 4 if inconsistent orientation
        if (orientation == 'inconsistent'):    
            score = 4

	    ## Inconsistent -> no TSD info available:
	    TSDlength = "na"
            TSDseq = "na"

        ## ------ Provisional -------
        ## Print results into the standard output
        print "TraFiC-id: ", traficId
        print "Score: ", score    
        print "Bkp-plus: ", bkpPlus
        print "Bkp-minus", bkpMinus
        print "TSD-length: ", TSDlength
        print "TSD-seq: ", TSDseq 
        print "Orientation: ", orientation
        print "Structure: ", structure
        print "TE-length: ", length
        print "perc-Length: ", percLength
        print "Poly-A: ", polyA
        print "contig-Plus: ", contigPlus
        print "contig-Minus: ", contigMinus
        
        ## Print results into an output file
        fileName = "TEIBA.results.txt"
        outFilePath = outDir + "/" + fileName
        outFile = open( outFilePath, "a" )

        row = traficId + "\t" + str(score) + "\t" + str(bkpPlus) + "\t" + str(bkpMinus) + "\t" + str(TSDlength) + "\t" + TSDseq + "\t" + orientation + "\t" + structure + "\t" + str(length) + "\t" + str(percLength) + "\t" + polyA + "\t" + contigPlus + "\t" + contigMinus + "\n"
        outFile.write(row)
        
        # Close output and end
        outFile.close()
        
        
    def insertion_orientation(self, typePlus, typeMinus):
        """ 
            Determine TE insertion strand/orientation.
            
            + orientation:
                5' informative (+ cluster)      ####TE####---------------
                3' informative (- cluster)      ---------------AAAAAAAAAA
            
            - orientation (the opposite)
                3' informative (+ cluster)      ---------------AAAAAAAAAA 
                5' informative (- cluster)      ####TE####---------------

            Input:
            1) typePlus. One of 'none', '5-prime' or '3-prime'
            2) typeMinus. One of 'none', '5-prime' or '3-prime'
            
            Output:
            1) orientation. +, -, 'na' or 'inconsistent' (informative 5' or 3' for both + and - clusters)
        """
        
        # A) No informative contig in + nor - clusters
        if (typePlus == "none") and (typeMinus == "none"):
            orientation = "na"
        
        # B) Informative contig of the same type for + and - clusters 
        elif (typePlus == typeMinus):
            orientation = "inconsistent"
        
        # C) 5-prime informative for + cluster and 3-prime informative 
        # or not informative for - cluster 
        elif (typePlus == "5-prime"):
            orientation = "+"
        
        # D) 5-prime informative for - cluster and 3-prime informative 
        # or not informative for + cluster
        elif (typeMinus == "5-prime"):
            orientation = "-"
        
        # E) 3-prime informative for + cluster and not informative for - cluster
        elif (typePlus == "3-prime"):
            orientation = "-"
        
        # F) 3-prime informative for - cluster and not informative for + cluster
        elif (typeMinus == "3-prime"):
            orientation = "+"
        
        return orientation
    
    
    def insertion_structure(self, typePlus, infoPlus, typeMinus, infoMinus, orientation):
        """
            ..
        
            Input:
            1) typePlus. One of 'none', '5-prime' or '3-prime'.
            2) infoPlus. Poly-A sequence for 3' and TE alignment object for 5' 
            3) typeMinus. One of 'none', '5-prime' or '3-prime'
            4) infoMinus. Poly-A sequence for 3' and TE alignment object for 5' 
            5) orientation. +, -, 'na' or 'inconsistent' 
            
            Output:
            1) structure. One of 'na', '5'inverted', '5'truncated' or 'full-length'
            2) length. Inserted TE length, 'na' if not available.
            3) percLength. Percentage of TE consensus sequence inserted, 'na' if not available.  
        """
        
        ## 1. Gather information about contig alignment on the TE sequence (L1, Alu or SVA)
        
        ## 1.A) Informative 5' contig for + or - cluster
        if ((typePlus == "5-prime") or (typeMinus == "5-prime")) and (orientation != "inconsistent"):
        
            # 1.A.a) + cluster informative 5'
            if (typePlus == "5-prime"):
                tSize = infoPlus.tSize
                tBeg = infoPlus.tBeg 
                tEnd = infoPlus.tEnd  
                strand = infoPlus.strand
                structure = ""
        
            # 1.A.b) - cluster informative 5'
            elif (typeMinus == "5-prime"):
                tSize = infoMinus.tSize
                tBeg = infoMinus.tBeg 
                tEnd = infoMinus.tEnd 
                strand = infoMinus.strand
                structure = ""
            
        ## 1.B) No 5' informative contigs for + nor - clusters or inconsistent orientation (both + and - are 5' informative)
        else:
            tSize = "na"
            tBeg = "na" 
            tEnd = "na"
            strand = "na"
            structure = "na"
            length = "na"
            percLength = "na"
            
        ## 2. Determine TE insertion structure
        ## 2.A) L1 inverted in its 5'
        # TE in + orientation with 5'inversion    ----#######TE######AAAAA----
        #                                             <<<<<<>>>>>>>>>>
        #                                             <---->
        #                                            inversion 
        # 5-prime informative contig              --------- (5'inversion signature: the piece of contig corresponding to L1 
        #                                                    aligns in the opposite DNA strand than the TE insertion orientation)
        if (strand != orientation) and (structure != "na"):
            structure = "5'inverted"
            length = "na"
            percLength = "na"
            
        ## 2.B) L1 full length or 5' truncated
        # Full length L1                         ----#######TE######AAAAA----
        # 5' truncated L1                        --------###TE######AAAAA----
        #                                            <-->
        #                                          deletion
        # 5-prime informative contig              ---____--- (5'truncation signature: the piece of contig corresponding to L1 
        #                                                    aligns in the body of the L1 and not in the 5' extreme) 
        elif (structure != "na"):
            length = tSize - tBeg
            percLength = float(length) / tSize * 100
            
            # 2.B.a) full length TE insertion                  
            if (length > 6000):
                # Threshold set for L1 (6021 bp length, first ~300bp correspond to promoter.
                # For Alu and SVA we need to put different values
                structure = "full-length"
            
            # 2.B.b) 5' truncated 
            else:
                structure = "5'truncated"
        
        return (structure, length, percLength)
            
    def polyA(self, typePlus, infoPlus, typeMinus, infoMinus):
        """
            Report poly-A sequence
        
            Input:
            1) typePlus. One of 'none', '5-prime' or '3-prime'
            2) infoPlus. Poly-A sequence for 3' and TE alignment object for 5' 
            3) typeMinus. One of 'none', '5-prime' or '3-prime'
            4) infoMinus. Poly-A sequence for 3' and TE alignment object for 5' 
            
            Output:
            1) polyA. Poly-A sequence. 
        """
        
        # A) Informative 3-prime + and - clusters 
        if (typePlus == "3-prime") and (typeMinus == "3-prime"):
            polyA = "inconsistent"
            
        # B) Informative 3-prime + cluster
        elif (typePlus == "3-prime"):
            polyA = infoPlus
            
        # C)  Informative 3-prime - cluster
        elif (typeMinus == "3-prime"):
            polyA = infoMinus
        
        # D) No informative 3-prime clusters
        else:
            polyA = "na"
            
        return polyA
                  
    def find_TSD(self, informativeContigPlusObj, informativeContigMinusObj):
        """ 
            Determine Target Site Duplication (TSD) size:
                     
            --------------------------- bkpPlus
                    bkpMinus --------------------
                             <-------->
                                 TSD (7bp)
        
            Input:
            1) bkpPlus. Breakpoint coordinates according to + cluster informative contig (tuple: chrom and pos).
            2) bkpMinus. Breakpoint coordinates according to - cluster informative contig (tuple: chrom and pos).
                     
            Output:
            1) TSDlength. Target site duplication length
            2) TSDseq. Target site duplication sequence
        """
        
        bkpPosPlus = informativeContigPlusObj.informative[1][1]
        bkpPosMinus = informativeContigMinusObj.informative[1][1]
        alignObjPlus = informativeContigPlusObj.informative[3]
        
        # A) Target site duplication (TSD)
        if (bkpPosPlus > bkpPosMinus):
            
            ## Compute TSD length
            TSDlength = bkpPosPlus - bkpPosMinus
        
            ## Extract TSD sequence
            # A) Begin of the contig sequence aligned in the TE insertion genomic region
            #   -------------**TSD**######TE#####
            #   --------------------
            # qBeg               *qEnd*
            if (alignObjPlus.alignType == "beg"):
                beg = alignObjPlus.qEnd - TSDlength
                end = alignObjPlus.qEnd
                TSDseq = informativeContigPlusObj.seq[beg:end]
        
            # B) End of the contig sequence aligned in the TE insertion genomic region
            #   ######TE#####AAAAAAA**TSD**-------------
            #                       --------------------
            #                    *qBeg*               qEnd
            elif (alignObjPlus.alignType == "end"):
                beg = alignObjPlus.qBeg 
                end = alignObjPlus.qBeg + TSDlength
                TSDseq = informativeContigPlusObj.seq[beg:end]           
                
        # B) No target site duplication
        else:
            TSDlength = 0
            TSDseq = "na"
            
        return (TSDlength, TSDseq)
    
    def imprecise_bkp(self, insertionCoord):
        """ 
            Compute confidence interval for imprecise breakpoints. 
            
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
        
        insertionCoordList = insertionCoord.split("_")
                
        chrom = str(insertionCoordList[0])
        beg = int(insertionCoordList[1])
        end = int(insertionCoordList[2])
        
        meanPos = (beg + end)/2
        dist = abs(beg - end)
        CI = dist/2
        
        breakpoint = [chrom, meanPos, CI]
    
        return breakpoint
    
    
class cluster():
    """ 
    Transposable element insertion cluster class. 
    
    A cluster can be + or - and it has associated the contigs resulting from the assembly of TE insertion supporting
    reads identified by TraFiC.  
    
    Methods:
    - fasta_reader
    - blat_alignment_reader
    - create_contigs_dict
    - add_alignments
    - find_informative_contig
    """
    
    def __init__(self, ID, contigsFasta):
        """ 
            Initialize cluster object.
            
            Input:
            1) contigsFasta. Fasta file containing the assembled contigs for the given cluster
            
            Output:
            - Cluster object variables initialized
        """
        self.ID = ID
        self.contigsDict = self.create_contigs_dict(contigsFasta)
        

    #### FUNCTIONS ####
    def fasta_reader(self, contigsFasta):
        """ 
            Read fasta file and produce an object generator of tuples (sequenceId,sequence).
            
            Input:
            1) contigsFasta. Fasta file containing the assembled contigs for the given cluster
            
            Output:
            1) object generator of tuples (sequenceId, sequence)
        """
        fh = open(contigsFasta)
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            yield header, seq
    
    def blat_alignment_reader(self, blatPath):
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
    
    def create_contigs_dict(self, contigsFasta):
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
        contigs = self.fasta_reader(contigsFasta)
 
        contigsDict = {}
    
        ### For each contig create a contig object and add it to the dictionary
        # using the contig id as key
        for contigTuple in contigs:
            
            # Create contig object
            contigObj = contig(contigTuple)
            
            # Add contig object to the dictionary 
            contigsDict[contigObj.ID] = contigObj
        
        return contigsDict
    
    def add_alignments(self, blatPath):
        """ 
            Read a psl file containing the contig blat aligments on the reference genome and TE sequences and associate
            each alignment to the corresponding contig object. 
            
            Input:
            1) blatPath. Psl file containing blat aligments for the assembled contigs.
            
            Output:
            1) For each contig object in the dictionary sets the 'alignList' variable. 
               This variable holds the list of alignment objects corresponding to the given contig. 
            
        """
        alignmentsDict = self.blat_alignment_reader(blatPath)
        
        ## Add the alignment lists to their respective contig objects
        # Iterate over the alignments dictionary. 
        for contigId in alignmentsDict:
            alignmentList = alignmentsDict[contigId]
            self.contigsDict[contigId].alignList = alignmentList  

    def find_informative_contig(self, insertionCoord):
        """
            Identify cluster informative contig if any. An informative cluster spans 5' or 3' insertion breakpoints.
            
            Input:
            insertionCoord. Region of interest. Format: ${chrom}_${beg}_${end}. 
                               Example: 10_108820680_108820678.
            
            Output:
            1) informativeContigObj. Informative contig object
        """
        
        informativeContigObj = ""
           
        info(str(len(self.contigsDict)) + ' input contigs') 
        
        ## Iterate over each contig object checking it is informative or not
        for contigId in self.contigsDict:    
            contigObj = self.contigsDict[contigId]
    
            # Check if it is an informative contig
            informative = contigObj.is_informative(insertionCoord)
    
            # A) Informative contig
            if (informative ==  1):
                message = contigObj.ID + " " + str(contigObj.informative[0]) + " " + contigObj.seq + " " + str(contigObj.informative[1]) + " " + str(contigObj.informative[2])
                log("INFORMATIVE", message)               
                informativeContigObj = contigObj
                
                # Stop iteration when informative contig found
                break
            
            # B) Not informative contig
            else:    
                message = contigObj.ID + " " + str(contigObj.informative[0]) + " " + contigObj.seq + " " + str(contigObj.informative[1]) + " " + str(contigObj.informative[2])
                log("NOT-INFORMATIVE", message)
                  
        return informativeContigObj
        
class contig():
    """ 
    Transposable element insertion contig class. 
    
    Contig sequence results from the assembly of TraFiC + or - cluster supporting reads 
    
    Methods:
    - is_candidate
    - is_informative
    - is_3prime_bkp
    - find_polyA
    - is_5prime_bkp
    """
    
    def __init__(self, contigTuple):
        """ 
            Initialize contig object.
            
            Input:
            1) contigTuple. First element (contig Id) and second element (contig Sequence)
            
            Output:
            - Contig object variables initialized
        """
        # Contig id and sequence
        self.ID = contigTuple[0]
        self.seq = contigTuple[1]
        self.length = int(len(self.seq))
        
        # Contig alignment information
        self.alignList = []   # List of blat alignments for this contig
        self.informative = [] # Four elements list: 
                              # 1) type (5-prime, 3-prime or none) 
                              # 2) breakpoint coordinates (tuple: chrom and pos, 'na' for both if type == none)
                              # 3) info dependent on the type. 5-prime: alignment object with contig's alignment in TE sequence info; 
                              #    3-prime: PolyA sequence; none: 'na' )
                              # 4) alignment object with contig's alignment in the target region info. 
                            
    #### FUNCTIONS ####
    def is_candidate(self, insertionCoord, windowSize, maxAlignPerc):
        """
        Check if contig is candidate to be informative about the TE insertion breakpoint. Informative contigs 
        span the insertion breakpoint. 

        Candidate contigs are defined as contigs partially aligning in the TE insertion region
         
        Input:
            1) insertionCoord. Region of interest. Format: ${chrom}_${beg}_${end}. 
                               Example: 10_108820680_108820678.
            2) windowSize (integer). Window size to extend from input region coordinates to define the
                                     region of interest
            3) maxAlignPerc (Float). Threshold to consider an alignment partial or not. 
                                     Partial defined as % of aligned contig sequence < maxAlignPerc.
            
        Output: 
            1) candidate. Boolean, 1 (informative candidate) and 0 (not informative candidate) 
            2) supportingAlignList. List of alignment objects supporting the contig as informative candidate.            
        """
        
        candidate = 0
        supportingAlignList = []
        
        # Iterate over all the contig blat alignments
        for alignment in self.alignList:
            
            # 1. Check if alignment within the target region
            insertionRegion = alignment.in_target_region(insertionCoord, windowSize)
            
            # Within target region
            if (insertionRegion == 1):
              
                # 2. Check if it is a partial alignment
                partial = alignment.partial_alignment(maxAlignPerc)
                
                # Partial
                if (partial == 1):            
                    supportingAlignList.append(alignment)
                    candidate = 1
                    
        return (candidate, supportingAlignList)

    
    def is_informative(self, insertionCoord):
        """
        Check if candidate contig is 5' or 3' informative. Defined as contigs spanning 
        5' or 3' insertion breakpoints:
        
        5' informative         ####TE####---------------
        3' informative         ---------------AAAAAAAAAA
        
        
        Input:
            1) insertionCoord. Region of interest. Format: ${chrom}_${beg}_${end}. 
                               Example: 10_108820680_108820678.
                            
        Output: 
            1) informativeBoolean. Boolean, 1 (informative) and 0 (not informative)
            2) Sets 'informative' variable. Three elements list:
                2.1) type (5-prime, 3-prime or none)
                2.2) breakpoint coordinates (tuple: chrom and pos, 'na' for both if type == none)
                2.3) info dependent on the type. 5-prime: aligment object iwth contig's alignment in TE sequence info; 
                     3-prime: PolyA sequence; none: 'na'
                2.4) alignment object with contig's alignment in the target region info. 
        """
            
        # 1) Check if it is an informative candidate contig
        candidate, supportingAlignList = self.is_candidate(insertionCoord, int(5000), float(98))
            
        # IF Informative contig candidate
        if (candidate == 1):
            
            # Iterate over the alignments supporting it as an informative contig
            for alignObj in supportingAlignList:  
                
                # 2) Check if candidate contig span 3' breakpoint (have polyA tail)
                is3prime, bkpCoord, polyASeq = self.is_3prime_bkp(alignObj)
                
                # 2.A) Candidate contig is informative 3-prime, it has a polyA tail 
                if (is3prime == 1):
                    informativeBoolean = 1
                    self.informative = ["3-prime", bkpCoord, polyASeq, alignObj]
                    
                    break
                
                # 2.B) Candidate contig is not informative 3-prime
                else:
                    
                    # 3) Check if candidate contig span 5' breakpoint (remaining sequece is maps in L1
                    is5prime, bkpCoord, TEalignmentObj = self.is_5prime_bkp(alignObj)
                    
                    # 3.A) Candidate contig is informative 5-prime, it has TE sequence 
                    if (is5prime == 1):
                        informativeBoolean = 1
                        self.informative = ["5-prime", bkpCoord, TEalignmentObj, alignObj]
                        
                        break
                        
                    # 3.B) Candidate contig is not informative 3-prime neither. 
                    # Not informative contig
                    else:
                        informativeBoolean = 0
                        self.informative = ["none", bkpCoord, "na", "na"]
        else:
            informativeBoolean = 0
            self.informative = ["none", ("na", "na"), "na"]
            
        return informativeBoolean
                
        
    def is_3prime_bkp(self, alignObj):
        """
        Check if candidate contig is 3' informative. Defined as contigs spanning 
        polyA - insertion target region breakpoint:
        
        3' informative         ---------------AAAAAAAAAA
        3' informative         AAAAAAAAAA---------------
        
        Input:
            1) alignObj. Blat alignment object.
                            
        Output: 
            1) is3prime. Boolean, 1 (3' informative) and 0 (not 3' informative)
            2) bkpCoord. Two elements breakpoint coordinates list. First (bkp chromosome) and second (breakpoint position)
            3) polyASeq. PolyA sequence.            
        """
        
        ## Select contig target piece of sequence to search for poly-A. 
        # The position of the target sequence in the contig 
        # will depend on the blat alignment type
        
        # A) Begin of the contig sequence aligned in the TE insertion genomic region
        #   -------------AAAAAAAAAAAA.....
        #   -------------
        # qBeg       *qEnd*
        if (alignObj.alignType == "beg"):
            targetSeq = self.seq[alignObj.qEnd:]
            
            bkpChrom = alignObj.tName
                
            if (alignObj.strand == "+"):
                bkpPos = alignObj.tEnd
            else:
                bkpPos = alignObj.tBeg
                
            ## Search for poly-A in the contig target piece of sequence. 
            polyASeq = self.find_polyA(targetSeq, 10)
        
        # B) End of the contig sequence aligned in the TE insertion genomic region
        #   ....AAAAAAAAAAAA-------------
        #                   -------------
        #                *qBeg*        qEnd
        elif (alignObj.alignType == "end"):
            targetSeq = self.seq[:alignObj.qBeg]            
            targetSeq = targetSeq[::-1] # Make the reverse to have the poly-A at the beginning if exists
        
            bkpChrom = alignObj.tName
            
            if (alignObj.strand == "+"):
                bkpPos = alignObj.tBeg
            else:
                bkpPos = alignObj.tEnd
                
            ## Search for poly-A in the contig target piece of sequence. 
            polyASeq = self.find_polyA(targetSeq, 10)
            
            polyASeq = polyASeq[::-1] # Make the reverse to put the poly-A in its original order
            
        # C) No align type information or 'none' align type
        else:
            log("Error", "No valid alignment object provided. Alignment type variable is 'none' or not defined") 
            sys.exit(1)
        
        ## A) Poly-A sequence
        if (polyASeq != ""):
            is3prime = 1
            bkpCoord = [bkpChrom, bkpPos]
        
        ## B) No poly-A sequence
        else:    
            is3prime = 0
            bkpCoord = ["na", "na"]
            polyASeq = "na"
            
            
        return (is3prime, bkpCoord, polyASeq)
            
        
    def find_polyA(self, targetSeq, windowSize):
        """
        Identify polyA(T) tail in a target sequence. 
        
        Parse the sequence doing consecutive sliding windows and computing the percentage of T and A for each window. 
        If it is >80 classify window as poly-A and search in the next window. If no poly-A window stop parsing. 
        
        Input:
        1) targetSeq. Target DNA sequence to search for poly A. 
        2) windowSize (integer). Sliding window size to parse target sequence searching for poly A:
                                                    
                                    AAAAAAAAAAAAAAAAAAAAAAAA
                                    *---W--->---W--->---W--->
        
        Output:
        1) polyASeq. PolyA sequence.            
        """
        # Convert sequence into upper case:
        targetSeq = targetSeq.upper()
        
        ## Search for polyA tails in the target seq. taking consecutive slices of X bases (X = windowSize)
        polyASeq = ""
        for beg in range(0, len(targetSeq), windowSize):

            # Take slice
            end = beg + windowSize
            seq = targetSeq[beg:end]
        
            # Compute percentage of A and T in the given slice
            nbA = seq.count("A") 
            nbG = seq.count("G") 
            nbC = seq.count("C") 
            nbT = seq.count("T") 
            nbN = seq.count("N")
       
            percA = float(nbA) / (nbA + nbG + nbC + nbT + nbN) * 100
            percT = float(nbT) / (nbA + nbG + nbC + nbT + nbN) * 100
            
            # A) Classify the slide as poly-A if > 80% are T or A
            if (percA >= 80) or (percT >= 80):
                
                # Add sequence slice to the polyA sequence
                polyASeq = polyASeq + seq 
            
            # B) Stop parsing sequence if the slide is not poly-A 
            else:
                    break
                
        return polyASeq
      
    def is_5prime_bkp(self, alignObj):
        """
        Check if candidate contig is 5' informative. Defined as contigs spanning 
        TE sequence - insertion target region breakpoint:
        
        5' informative:         ---------------####TE####
        5' informative:         ####TE####---------------
        
        Input:
        1) alignObj. Blat alignment object.
                            
        Output:
        1) is5prime. Boolean, 1 (5' informative) and 0 (not 5' informative).
        2) bkpCoord. Two elements breakpoint coordinates list. First (bkp chromosome) and second (breakpoint position).
        3) TEalignmentObj. Blat aligment object with the alignment information of the contig in the consensus TE sequence.
                           'na' if not 5' informative.
        """
        
        ## Select contig target sequence coordinates to search for alignment in L1. 
        # The position of the target coordinates in the contig 
        # will depend on the blat alignment type
        
        # A) Begin of the contig sequence aligned in the TE insertion genomic region
        #   -------------******TE******
        #   -------------
        # qBeg        *qEnd*
        if (alignObj.alignType == "beg"):
            targetBeg = alignObj.qEnd
            targetEnd = alignObj.qSize
            
            bkpChrom = alignObj.tName
            
            if (alignObj.strand == "+"):
                bkpPos = alignObj.tEnd
            else:
                bkpPos = alignObj.tBeg
                
        # B) End of the contig sequence aligned in the TE insertion genomic region
        #   ******TE******-------------
        #                 -------------
        #              *qBeg*        qEnd
        elif (alignObj.alignType == "end"):
            targetBeg = 0
            targetEnd = alignObj.qBeg
            
            bkpChrom = alignObj.tName
            
            if (alignObj.strand == "+"):
                bkpPos = alignObj.tBeg   
            else:
                bkpPos = alignObj.tEnd
                
        # C) No align type information or 'none' align type
        else:
            log("Error", "No valid alignment object provided. Alignment type variable is 'none' or not defined") 
            sys.exit(1)
     
        
        ## Default
        is5prime = 0
        bkpCoord = ["na", "na"]
        TEalignmentObj = "na" 
                
        for alignment in self.alignList:
            
            # Contig alignment in L1
            if ( alignment.tName == "L1" ):
                
                ## Compute percentage of overlap between: 
                # Expected alignment ---------------
                #                targetBeg      targetEnd 
                
                # TE alignment           ***************           
                #                      qBeg           qEnd 
                beginList = [ targetBeg, alignment.qBeg ]
                endList = [ targetEnd, alignment.qEnd ]
                
                length = float(max(endList) - min(beginList)) 
                nbOverlapingBases = float(min(endList) - max(beginList))                
                percOverlap = nbOverlapingBases / length * 100
                
                # If percentage of overlap > 50% -> Informative 5'
                if ( percOverlap > 50 ):   
                    is5prime = 1
                    bkpCoord = [bkpChrom, bkpPos]
                    TEalignmentObj = alignment 
                    break 
        
        return (is5prime, bkpCoord, TEalignmentObj)
        
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
    
    #### FUNCTIONS ####
    def in_target_region(self, coords, windowSize):
        """ 
            Check if blat alignment within a region of interest:
            
                region     chrom   beg ------------- end
                alignment  chrom  tBeg     ------    tEnd
            
            Input:
            1) coords. Region of interest. Format: ${chrom}_${beg}_${end}. 
                       Example: 10_108820680_108820678.
            2) windowSize (integer). Window size to extend from input region coordinates to define the
                                     region of interest:  
                                     
                                     <--W-->---Input---<--W-->
                                      -----------------------
                                         region of interest
            Output:
            1) insertionRegion. Boolean, 1(in region) and 0 (outside of the region)
        """
        coordsList = coords.split("_")
        chrom = coordsList[0] 
        beg = int(coordsList[1]) - windowSize
        end = int(coordsList[2]) + windowSize
         
        # A) Within target region
        if (chrom == self.tName) and (self.tBeg >= beg) and (self.tEnd <= end):           
            insertionRegion = 1
        
        # B) Outside target region    
        else:
            insertionRegion = 0
        
        return insertionRegion
    
    def partial_alignment(self, maxAlignPerc):
        """ 
            Check if blat alignment is partial or not and classify partial alignments in one 
            of those categories: "beg" and "end". 
            
                          firstHalf    secondHalf
            contig_seq:  ------------*------------       
            partial_beg:      ---------
                             -------------

            partial_end:                ----------
                                    ------------
             
            Input:
            1) maxAlignPerc (Float). Threshold to consider an alignment partial or not. 
                                     Partial defined as % of aligned contig sequence < maxAlignPerc.
                            
            Output:
            1) partial. Boolean, 1(partial) and 0 (not partial).
            2) Sets 'alignType' variable to "none", "beg" or "end"          
        """
        
        alignLength = self.qEnd - self.qBeg
        alignPerc = float(alignLength) / float(self.qSize) * 100  
        
        # A) Partial alignment
        if (alignPerc < maxAlignPerc):
            partial = 1
            
            ## Determine type of partial alignment (beg, end, none):
            middle = float(self.qSize)/2  
            
            # a) Begin of the alignment within the first half of contig sequence:
            if self.qBeg <= middle: 
                    
                # a.a) End of the alignment within the first half of contig sequence 
                if self.qEnd <= middle:
                    self.alignType = "beg"
                
                # a.b) End of the alignment within the second half of contig sequence
                else:     
                    dist2Beg = self.qBeg
                    dist2End = self.qSize - self.qEnd 
                    
                    if (dist2Beg <= dist2End):
                        self.alignType = "beg"
                    else:
                        self.alignType = "end"
                        
            # b) Begin of the alignment within the second half of contig sequence:
            else: 
                self.alignType = "end"
                
        # B) No partial
        else:    
            partial = 0
            self.alignType = "none"
        
        return partial
                

#### MAIN ####

## Import modules ##
import argparse
import time
import sys
from itertools import groupby
import os.path

## Get user's input ## 
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('paths2bkpAnalysis', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputPath = args.paths2bkpAnalysis
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "paths2bkpAnalysis: ", inputPath
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, " *****"
print 
print "..."
print 


## Per each insertion perform breakpoint analysis ## 

inputFile = open(inputPath, 'r')

# Analyze one insertion per iteration
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")
    
    # Get TE insertion info and files
    category, insertionCoord = line[0].split(":")
    contigsPlusPath, contigsMinusPath = line[1].split(",")
    blatPlusPath, blatMinusPath = line[2].split(",")

    # Perform breakpoint analysis for the TE insertion 
    header("Tranposable Element Insertion Breakpoint Analysis (TEIBA) for: " + insertionCoord)
    
    # A) All the input files exist 
    if os.path.isfile(contigsPlusPath) and os.path.isfile(blatPlusPath) and os.path.isfile(contigsMinusPath) and os.path.isfile(blatMinusPath):  
        insertionObj = insertion(category, insertionCoord, contigsPlusPath, blatPlusPath, contigsMinusPath, blatMinusPath)
        insertionObj.find_insertionBkp(outDir)
    else:
        message = "Input files for " + insertionCoord + " insertion do not exist"
        log("ERROR", message)

## Finish ##
print 
print "***** Finished! *****"
print 


