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

def insertionRange(MEIObj, overhang):
    """
    Define MEI range as follows:
                
                      range
       <------------------------------------>
       beg [bkpA - CIPOS - overhang]        end [(bkpB or bkpA) - CIPOS - overhang]  
    """
    bkpA = MEIObj.pos
    bkpB = MEIObj.infoDict["BKPB"] if "BKPB" in MEIObj.infoDict else MEIObj.pos 
    CIPOS = MEIObj.infoDict["CIPOS"]
    begRange = int(bkpA) - int(CIPOS) - int(overhang)
    endRange = int(bkpB) + int(CIPOS) + int(overhang)

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

def clusterMEI(MEIDict):
    """
    """
    
    #msg = "** CLUSTER MEI **"  
    #log("CLUSTER", msg)

    # Cluster MEI objects
    #  chr1 --------*---------------*-------------*------------
    #           cluster1        cluster2      cluster3     
    #         (MEI1,MEI2,MEI3)   (MEI4)       (MEI5,MEI6)       
    clusterList = []

    # Iterate over the dictionary picking a MEI list in each iteration:
    for MEIclass in MEIDict:
        for chrom in MEIDict[MEIclass]:
        
            MEIlist = MEIDict[MEIclass][chrom]
            
            # Sort list of MEI objects from lower to upper coordinates
            MEIlistSorted = sorted(MEIlist, key=lambda MEIObj: MEIObj.pos)

            # Per MEI object:        
            for MEIObj in MEIlistSorted:

                bkpB = int(MEIObj.infoDict["BKPB"]) if "BKPB" in MEIObj.infoDict else "UNK"
        
                #msg = "MEI: " + chrom + " " + str(MEIObj.pos) + " " + str(bkpB) + " " + MEIObj.infoDict["TYPE"] + " " + MEIObj.infoDict["CLASS"] + " " + MEIObj.infoDict["SCORE"] + " " + MEIObj.infoDict["CIPOS"] 
                #log("CLUSTER", msg)    
        
                # A) No cluster in the list -> Create first cluster
                if not clusterList:

                    #msg = "Initialize first cluster"
                    #log("CLUSTER", msg) 
                    clusterObj = MEIcluster(MEIObj, 5) # Allow up to 5 nucleotides of margin
                    clusterList.append(clusterObj)
        
                # B) There is already at least one cluster in the list -> Check if current MEI within the latest cluster
                else:
                    
                    #msg = "Check if MEI within latest cluster"
                    #log("CLUSTER", msg) 

                    ## Define cluster range for searching for overlap
                    lastClusterObj = clusterList[-1]     
                    begClusterRange = lastClusterObj.beg
                    endClusterRange = lastClusterObj.end

                    ## Define insertion range for searching for overlap:    
                    begMEIrange, endMEIrange = insertionRange(MEIObj, 5)     
              
                    #### Check if both ranges overlap.
                    overlapping = overlap(begMEIrange, endMEIrange, begClusterRange, endClusterRange) 

                    #msg = "cluster_range,MEI_range: " + str(begClusterRange) + " " + str(endClusterRange) + " " + str(begMEIrange) + " " + str(endMEIrange)
                    #log("CLUSTER", msg) 

    
                    ## A) Overlapping ranges, so current MEI within previous cluster interval -> add MEI to the cluster                                
                    if overlapping:
                        #msg = "MEI within cluster -> add MEI to the cluster"
                        #log("CLUSTER", msg) 
                        lastClusterObj.addMEI(MEIObj, 5) # Allow up to 5 nucleotides of margin 
                     
                    ## B) Current MEI outside previous cluster interval -> create new cluster and add it into the list
                    else:
                        #msg = "MEI outside the cluster -> create new cluster "
                        #log("CLUSTER", msg) 
                        clusterObj = MEIcluster(MEIObj, 5) # Allow up to 5 nucleotides of margin
                        clusterList.append(clusterObj)
                    
                    #msg = "Number of MEI within cluster: ", len(clusterObj.MEIlist)
                    #log("CLUSTER", msg) 
                    #msg = "----------------------"
                    #log("CLUSTER", msg) 

    return clusterList     



####### CLASSES #######
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

    def consensus(self):
        """
        Select consensus MEI. Selection criteria (ordered by preference order):
        1) Lowest CIPOS
        2) Highest total number of supporting reads (+ plus - cluster supporting reads)
                chr1 --------*---------------*-------------*------------
                   beg    cluster1        cluster2      cluster3      end
                     (MEI1,MEI2,MEI3)      (MEI4)     (MEI5,MEI6)
        consensus          ****             ****       ****
        """

        ## 1. Sort by total number of supporting paired-ends (decreasing order)
        sortedList1 = sorted(self.MEIlist, key=methodcaller('nbPE'), reverse=True)

        ## 2. Sort by CIPOS (ascending order)
        MEItuples = [(MEIobj, int(MEIobj.infoDict["CIPOS"])) for MEIobj in sortedList1]
        sortedMEItuples = sorted(MEItuples, key=itemgetter(1))
        sortedList2 = [i[0] for i in sortedMEItuples]

        ## 3. Select consensus MEI (the one with lowest CIPOS and highest total number of supporting paired-ends)
        consensusMEIobj = sortedList2[0]

        return consensusMEIobj

class cohort():
    """
        .....................

        Methods:
        -
    """
    def __init__(self):
        """
        """
        self.VCFlist = []
        self.MEIDict = {}
        self.consensusMEIDict = {}

    def read_VCFs(self, inputPath):
        """
        """
        inputFile = open(inputPath, 'r')

        info("Read input VCFs ")

        # Per iteration, read a VCF, generate a VCF object and add it to the cohort
        for line in inputFile:
            line = line.rstrip('\n')
            line = line.split("\t")

            donorId = line[0]
            VCFfile = line[1]

            # Create VCF object
            VCFObj = formats.VCF()

            info("Reading " + VCFfile + "...")

            # Input VCF available
            if os.path.isfile(VCFfile):

                # Read VCF and add information to VCF object
                VCFObj.read_VCF(VCFfile)

                # Add information to the genotype field in each MEI object
                for MEIObject in VCFObj.lineList:

                    ## remove sample list format field
                    formatList = MEIObject.format.split(':')
                    genotypeList = MEIObject.genotype.split(':')

                    ## add representative sample format field
                    MEIObject.format = formatList[0] + ':' + formatList[1] + ':REPR'
                    MEIObject.genotype = genotypeList[0] + ':' + genotypeList[1] + ':' + donorId

                # Add donor VCF to cohort
                self.addDonor(VCFObj)
            else:
                print "[ERROR] Input file does not exist"

    def addDonor(self, VCFobj):
        """
        """
        self.VCFlist.append(VCFobj)

    def organizeMEI(self):
        """
        Organize all the MEI, across the donors in the cohort, into a nested dictionary containing for each chromosome and insertion class a list of MEI objects that pass all the filters:
        key1(chrId) -> value(dict2) -> key1(MEIclass) -> value(insertionObjList)
        """

        for VCFobj in self.VCFlist:

            # Per MEI VCFline object in a given donor
            for MEIObj in VCFobj.lineList:

                ## Select only those insertions passing all the filters:
                if (MEIObj.filter == "PASS"):

                    # A) First MEI of a given class
                    if MEIObj.infoDict["CLASS"] not in self.MEIDict:
                        self.MEIDict[MEIObj.infoDict["CLASS"]] = {}
                        self.MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom] = [ MEIObj ]

                    # B) There are already MEI of this class
                    else:

                        # a) First MEI in the chromosome
                        if MEIObj.chrom not in self.MEIDict[MEIObj.infoDict["CLASS"]]:
                            self.MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom] = {}
                            self.MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom] = [ MEIObj ]

                        # b) There are already MEI in the chromosome
                        else:
                            self.MEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom].append(MEIObj)

    def consensus_MEIdict(self):
        """
        Create a nested dictionary containing for each chromosome and insertion class the list of consensus MEI objects:
        
        key1(chrId) -> value(dict2) -> key1(MEIclass) -> value(consensusMEIList)
        """

        # 1. Cluster MEI
        MEIclusterList = clusterMEI(self.MEIDict)
        
        # 2. Select a consensus MEI per cluster and add it to the list within the dictionary
        # Per MEI cluster:
        for MEIclusterObj in MEIclusterList:
            
            consensusMEIObj = MEIclusterObj.consensus()

            # A) First consensus MEI of a given class
            if consensusMEIObj.infoDict["CLASS"] not in self.consensusMEIDict:
                self.consensusMEIDict[consensusMEIObj.infoDict["CLASS"]] = {}
                self.consensusMEIDict[consensusMEIObj.infoDict["CLASS"]][consensusMEIObj.chrom] = [ consensusMEIObj ]

            # B) There are already consensus MEI of this class
            else:

                # a) First consensus MEI in the chromosome
                if consensusMEIObj.chrom not in self.consensusMEIDict[consensusMEIObj.infoDict["CLASS"]]:
                    self.consensusMEIDict[consensusMEIObj.infoDict["CLASS"]][consensusMEIObj.chrom] = {}
                    self.consensusMEIDict[consensusMEIObj.infoDict["CLASS"]][consensusMEIObj.chrom] = [ consensusMEIObj ]

                # b) There are already consensus MEI in the chromosome
                else:
                    self.consensusMEIDict[consensusMEIObj.infoDict["CLASS"]][consensusMEIObj.chrom].append(consensusMEIObj)


    def write_header(self, outFilePath):
        """
            ...

            Input:
            1) ...

            Output:
            1) ...
        """

        ## Define variables
        date = time.strftime("%Y%m%d")

        context = {
         "date": date,
         "source": "TraFiCv2.0",
         "reference": "hs37d5",
         }

        ## Header template
        template = """##fileformat=VCFv4.2
##fileDate={date}
##source={source}
##reference={reference}
##contig=<ID=1,assembly=GRCh37,length=249250621,species=human>
##contig=<ID=2,assembly=GRCh37,length=243199373,species=human>
##contig=<ID=3,assembly=GRCh37,length=198022430,species=human>
##contig=<ID=4,assembly=GRCh37,length=191154276,species=human>
##contig=<ID=5,assembly=GRCh37,length=180915260,species=human>
##contig=<ID=6,assembly=GRCh37,length=171115067,species=human>
##contig=<ID=7,assembly=GRCh37,length=159138663,species=human>
##contig=<ID=8,assembly=GRCh37,length=146364022,species=human>
##contig=<ID=9,assembly=GRCh37,length=141213431,species=human>
##contig=<ID=10,assembly=GRCh37,length=135534747,species=human>
##contig=<ID=11,assembly=GRCh37,length=135006516,species=human>
##contig=<ID=12,assembly=GRCh37,length=133851895,species=human>
##contig=<ID=13,assembly=GRCh37,length=115169878,species=human>
##contig=<ID=14,assembly=GRCh37,length=107349540,species=human>
##contig=<ID=15,assembly=GRCh37,length=102531392,species=human>
##contig=<ID=16,assembly=GRCh37,length=90354753,species=human>
##contig=<ID=17,assembly=GRCh37,length=81195210,species=human>
##contig=<ID=18,assembly=GRCh37,length=78077248,species=human>
##contig=<ID=19,assembly=GRCh37,length=59128983,species=human>
##contig=<ID=20,assembly=GRCh37,length=63025520,species=human>
##contig=<ID=21,assembly=GRCh37,length=48129895,species=human>
##contig=<ID=22,assembly=GRCh37,length=51304566,species=human>
##contig=<ID=hs37d5,assembly=GRCh37,length=35477943,species=human>
##contig=<ID=GL000191.1,assembly=GRCh37,length=106433,species=human>
##contig=<ID=GL000192.1,assembly=GRCh37,length=547496,species=human>
##contig=<ID=GL000193.1,assembly=GRCh37,length=189789,species=human>
##contig=<ID=GL000194.1,assembly=GRCh37,length=191469,species=human>
##contig=<ID=GL000195.1,assembly=GRCh37,length=182896,species=human>
##contig=<ID=GL000196.1,assembly=GRCh37,length=38914,species=human>
##contig=<ID=GL000197.1,assembly=GRCh37,length=37175,species=human>
##contig=<ID=GL000198.1,assembly=GRCh37,length=90085,species=human>
##contig=<ID=GL000199.1,assembly=GRCh37,length=169874,species=human>
##contig=<ID=GL000200.1,assembly=GRCh37,length=187035,species=human>
##contig=<ID=GL000201.1,assembly=GRCh37,length=36148,species=human>
##contig=<ID=GL000202.1,assembly=GRCh37,length=40103,species=human>
##contig=<ID=GL000203.1,assembly=GRCh37,length=37498,species=human>
##contig=<ID=GL000204.1,assembly=GRCh37,length=81310,species=human>
##contig=<ID=GL000205.1,assembly=GRCh37,length=174588,species=human>
##contig=<ID=GL000206.1,assembly=GRCh37,length=41001,species=human>
##contig=<ID=GL000207.1,assembly=GRCh37,length=4262,species=human>
##contig=<ID=GL000208.1,assembly=GRCh37,length=92689,species=human>
##contig=<ID=GL000209.1,assembly=GRCh37,length=159169,species=human>
##contig=<ID=GL000210.1,assembly=GRCh37,length=27682,species=human>
##contig=<ID=GL000211.1,assembly=GRCh37,length=166566,species=human>
##contig=<ID=GL000212.1,assembly=GRCh37,length=186858,species=human>
##contig=<ID=GL000213.1,assembly=GRCh37,length=164239,species=human>
##contig=<ID=GL000214.1,assembly=GRCh37,length=137718,species=human>
##contig=<ID=GL000215.1,assembly=GRCh37,length=172545,species=human>
##contig=<ID=GL000216.1,assembly=GRCh37,length=172294,species=human>
##contig=<ID=GL000217.1,assembly=GRCh37,length=172149,species=human>
##contig=<ID=GL000218.1,assembly=GRCh37,length=161147,species=human>
##contig=<ID=GL000219.1,assembly=GRCh37,length=179198,species=human>
##contig=<ID=GL000220.1,assembly=GRCh37,length=161802,species=human>
##contig=<ID=GL000221.1,assembly=GRCh37,length=155397,species=human>
##contig=<ID=GL000222.1,assembly=GRCh37,length=186861,species=human>
##contig=<ID=GL000223.1,assembly=GRCh37,length=180455,species=human>
##contig=<ID=GL000224.1,assembly=GRCh37,length=179693,species=human>
##contig=<ID=GL000225.1,assembly=GRCh37,length=211173,species=human>
##contig=<ID=GL000226.1,assembly=GRCh37,length=15008,species=human>
##contig=<ID=GL000227.1,assembly=GRCh37,length=128374,species=human>
##contig=<ID=GL000228.1,assembly=GRCh37,length=129120,species=human>
##contig=<ID=GL000229.1,assembly=GRCh37,length=19913,species=human>
##contig=<ID=GL000230.1,assembly=GRCh37,length=43691,species=human>
##contig=<ID=GL000231.1,assembly=GRCh37,length=27386,species=human>
##contig=<ID=GL000232.1,assembly=GRCh37,length=40652,species=human>
##contig=<ID=GL000233.1,assembly=GRCh37,length=45941,species=human>
##contig=<ID=GL000234.1,assembly=GRCh37,length=40531,species=human>
##contig=<ID=GL000235.1,assembly=GRCh37,length=34474,species=human>
##contig=<ID=GL000236.1,assembly=GRCh37,length=41934,species=human>
##contig=<ID=GL000237.1,assembly=GRCh37,length=45867,species=human>
##contig=<ID=GL000238.1,assembly=GRCh37,length=39939,species=human>
##contig=<ID=GL000239.1,assembly=GRCh37,length=33824,species=human>
##contig=<ID=GL000240.1,assembly=GRCh37,length=41933,species=human>
##contig=<ID=GL000241.1,assembly=GRCh37,length=42152,species=human>
##contig=<ID=GL000242.1,assembly=GRCh37,length=43523,species=human>
##contig=<ID=GL000243.1,assembly=GRCh37,length=43341,species=human>
##contig=<ID=GL000244.1,assembly=GRCh37,length=39929,species=human>
##contig=<ID=GL000245.1,assembly=GRCh37,length=36651,species=human>
##contig=<ID=GL000246.1,assembly=GRCh37,length=38154,species=human>
##contig=<ID=GL000247.1,assembly=GRCh37,length=36422,species=human>
##contig=<ID=GL000248.1,assembly=GRCh37,length=39786,species=human>
##contig=<ID=GL000249.1,assembly=GRCh37,length=38502,species=human>
##contig=<ID=MT,assembly=GRCh37,length=16569,species=human>
##contig=<ID=NC_007605,assembly=GRCh37,length=171823,species=human>
##contig=<ID=X,assembly=GRCh37,length=155270560,species=human>
##contig=<ID=Y,assembly=GRCh37,length=59373566,species=human>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant. (All sequence is on the plus strand and in the forward direction).">
##INFO=<ID=CLASS,Number=1,Type=String,Description="Mobile element class (L1, ALU, SVA or ERVK)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Insertion type (TD0: solo, TD1: partnered-3'transduction, TD2: orphan-3'transduction, PSD: processed-pseudogene)>
##INFO=<ID=SCORE,Number=1,Type=Integer,Description="Insertion score (5: 5' and 3' breakpoints (bkp) assembled, 4: 3'bkp assembled, 3: 5'bkp assembled, 2: no bkp assembled, 1: inconsistent (contradictory orientation, bkp or TSD))">
##INFO=<ID=MANUAL,Number=0,Type=Flag,Description="MEI manually verified and curated through BAM inspection (Only used for PSD)">
##INFO=<ID=BKPB,Number=1,Type=Integer,Description="MEI right-most breakpoint position (bkp B). Left-most breakpoint position (bkp A) represented in the POS field">
##INFO=<ID=CIPOS,Number=1,Type=Integer,Description="Confidence interval around insertion breakpoints">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion DNA strand (+ or -)">
##INFO=<ID=STRUCT,Number=1,Type=String,Description="Mobile element structure (INV: 5'inverted, DEL: 5'deleted, FULL: full-length)">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Mobile element length">
##INFO=<ID=TSLEN,Number=1,Type=Integer,Description="Target site duplication (+_value) or deletion (-_value) length">
##INFO=<ID=TSSEQ,Number=1,Type=String,Description="Target site duplication sequence">
##INFO=<ID=POLYA,Number=1,Type=String,Description="Poly-A sequence">
##INFO=<ID=SRCID,Number=1,Type=String,Description="Source element cytoband identifier. Only for gemline source elements">
##INFO=<ID=SRCTYPE,Number=1,Type=String,Description="Source element type (GERMLINE or SOMATIC)">
##INFO=<ID=SRC,Number=1,Type=String,Description="Coordinates of the source element that mediated the transduction in the format: chrom_beg_end">
##INFO=<ID=TDC,Number=1,Type=String,Description="Begin and end coordinates of the integrated transduced or pseudogene sequence in the format: chrom_beg_end">
##INFO=<ID=TDLEN,Number=1,Type=Integer,Description="Transduced region length">
##INFO=<ID=TDLENR,Number=1,Type=Integer,Description="Transduced region length at RNA level">
##INFO=<ID=SRCGENE,Number=1,Type=String,Description="Source gene of the processed pseudogene insertion">
##INFO=<ID=GERMDB,Number=1,Type=String,Description="MEI already reported as germinal in a database (1KGENOMES: 1000 genomes project (source_papers_doi: 10.1038/nature15394 and 10.1073/pnas.1602336113), TRAFIC: TraFic in-house database)">
##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region where the mobile element is inserted (exonic, splicing, ncRNA, UTR5, UTR3, intronic, upstream, downstream, intergenic)">
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO gene symbol">
##INFO=<ID=ROLE,Number=1,Type=String,Description="Role in cancer (oncogene, TSG: tumor suppressor gene, oncogene/TSG: both roles)">
##INFO=<ID=COSMIC,Number=0,Type=Flag,Description="Reported as cancer driver in COSMIC cancer gene census database">
##INFO=<ID=CPG,Number=0,Type=Flag,Description="Reported as cancer predisposition gene in 10.1038/nature12981 (DOI).">
##INFO=<ID=REP,Number=1,Type=String,Description="Repetitive element overlapping the insertion breakpoint">
##INFO=<ID=DIV,Number=1,Type=Integer,Description="Millidivergence of the overlapping repetitive element with respect a consensus sequence">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Assembled contig sequence spanning 1st bkp (lowest genomic position)">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Assembled contig sequence spanning 2nd bkp (highest genomic position)">
##INFO=<ID=RP,Number=.,Type=String,Description="Reads from the tumour sample and positive cluster that support this insertion">
##INFO=<ID=RN,Number=.,Type=String,Description="Reads from the tumour sample and negative cluster that support this insertion">
##FILTER=<ID=SCORE,Description="Insertion with an score < threshold">
##FILTER=<ID=REP,Description="Insertion overlapping a satellite region or a repetitive element of the same class">
##FILTER=<ID=DUP,Description="Duplicated MEI call">
##FILTER=<ID=GERMLINE,Description="Germline MEI miscalled as somatic">
##FILTER=<ID=TD,Description="L1 transduction incorrectly identified as a processed pseudogene insertion">
##FORMAT=<ID=RCP,Number=1,Type=Integer,Description="Count of positive cluster supporting reads">
##FORMAT=<ID=RCN,Number=1,Type=Integer,Description="Count of negative cluster supporting reads">
##FORMAT=<ID=SL,Number=1,Type=String,Description="List of samples where the variant was found (specially relevant for multi-tumor donors)">
##FORMAT=<ID=REPR,Number=1,Type=String,Description="Sample selected as representative among all the samples where the variant was found (specially relevant for multi-sample VCF).">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting the variant in this sample">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads covering variant location in this sample">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  
"""

        ## Replace variables into the template and print header into the output file
        with  open(outFilePath,'w') as outFile:
            outFile.write(template.format(**context))

    def write_consensus(self, outFilePath):
        """
        """

        outFile = open(outFilePath, 'a')

        ## 1. Convert dictionary into list of consensus MEI objects
        tmpList = []

        for MEIClass in self.consensusMEIDict:
            for chrom in self.consensusMEIDict[MEIClass]:
                tmpList.append( self.consensusMEIDict[MEIClass][chrom])

        consensusMEIList = [MEI for sublist in tmpList for MEI in sublist]

        ## 2. Sort list of consensus MEI first by chromosome and then by position:
        consensusMEIList.sort(key=lambda line: (line.chrom, line.pos))

        # Iterate and print each consensus germline MEI into the output VCF file
        for MEI in consensusMEIList:

            row = MEI.chrom + "\t" + str(MEI.pos) + "\t" + MEI.id + "\t" + MEI.ref + "\t" + MEI.alt + "\t" + MEI.qual + "\t" + MEI.filter + "\t" + MEI.info + "\t" + MEI.format + "\t" + MEI.genotype + "\n"
            outFile.write(row)

        ## Close output file
        outFile.close()

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
parser = argparse.ArgumentParser(description="Takes a set of input VCFs and merge them into a non-redundant VCF containing for those MEI shared between multiple VCFs only the one selected as representative.")
parser.add_argument('inputPath', help='Tabular text file containing one row per donor with the following consecutive fields: donorId vcf_path')
parser.add_argument('fileName', help='File name of the merged VCF generated as output')
parser.add_argument('--overhang', default=5, type=int, dest='overhang', help='Maximum overhang for MEI clustering. Default: 5 base pairs.')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
fileName = args.fileName
overhang = args.overhang
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "fileName: ", fileName
print "overhang: ", overhang
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ##

## 1. Initialize cohort object
header("1. Initialize cohort object")
cohortObj = cohort()

## 2. Read VCF files and make a list of VCF objects. Add donorId information
header("2. Read VCF files and make a list of VCF objects. Add donorId information")
cohortObj.read_VCFs(inputPath)

## 3. Organize MEI by insertion class, chromosome and in increasing cromosomal coordinates
header("3. Organize MEI by insertion class, chromosome and in increasing cromosomal coordinates")
cohortObj.organizeMEI()

## 4. Make not redundant list of consensus MEI objects and organize them into a dictionary
header("4. Make not redundant list of consensus MEI objects and organize them into a dictionary")
cohortObj.consensus_MEIdict()

## 5. Make output VCF containing consensus MEI objects
header("5. Make output VCF containing consensus MEI objects")
outFilePath = outDir + '/' + fileName +'.vcf'

# 5.1 Write header
cohortObj.write_header(outFilePath)

# 5.2 Write variants
cohortObj.write_consensus(outFilePath)

## End ##
print
print "***** Finished! *****"
print

