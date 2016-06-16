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

class VCF():
    """ 
	.....................
    
	Methods:
	- 

    """

    def __init__(self):
	""" 
	    ...
            
	    Output:
            - 
	"""
	self.lineList = []  # List of VCFline objects
        	
    #### METHODS ####
    def addLine(self, VCFlineObj):
	""" 
	    ...
            
	    Input:
	    1) ...
            
	    Output:
	    1) ...
	"""
	self.lineList.append(VCFlineObj)


    def print_header(self, outFilePath):
	""" 
	    ...
            
	    Input:
	    1) ...
            
	    Output:
	    1) ...
	"""
		
	## 1. Define variables 
	date = time.strftime("%Y%m%d")
		
	context = {
	 "date":date, 
	 } 

	## 2. Header template
  	template = """##fileformat=VCFv4.2
##fileDate={date} 
##source=TraFiCvX 
##reference=hs37d5
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
##INFO=<ID=CLASS,Number=4,Type=String,Description="Transposable element class (L1, ALU, SVA or HERVK)">
##INFO=<ID=TYPE,Number=3,Type=String,Description="Insertion type (TD0: solo, TD1: partnered-3'transduction, TD2: orphan-3'transduction)">
##INFO=<ID=SCORE,Number=4,Type=String,Description="Insertion score (1: 5' and 3' breakpoints (bkp) assembled, 2A: 5'bkp assembled, 2B: 3'bkp assembled, 3: no bkp assembled, 4: inconsistent (contradictory orientation, bkp or TSD))">
##INFO=<ID=CIPOS,Number=1,Type=Integer,Description="Confidence interval around POS (insertion breakpoint)">
##INFO=<ID=STRAND,Number=2,Type=String,Description="Insertion DNA strand (+ or -)">
##INFO=<ID=STRUCT,Number=3,Type=String,Description="Transposable element structure (INV: 5'inverted, DEL: 5'deleted, FULL: full-length)">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Transposable element length">
##INFO=<ID=PLEN,Number=1,Type=Float,Description="Percentage of consensus transposable element sequence inserted">
##INFO=<ID=TSLEN,Number=1,Type=Integer,Description="Target site length (+: target site duplication, -: target site microdeletion)">
##INFO=<ID=TSDSEQ,Number=1,Type=Integer,Description="Target site duplication sequence">
##INFO=<ID=POLYA,Number=1,Type=Integer,Description="Poly-A sequence">
##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region where the transposable element is inserted (exonic, splicing, ncRNA, UTR5, UTR3, intronic, upstream, downstream, intergenic)">
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO gene symbol">
##INFO=<ID=TID,Number=1,Type=String,Description="Transcript id for transcript associated with insertion breakpoint">
##INFO=<ID=SAT,Number=1,Type=String,Description="Satellite region overlapping insertion breakpoint">
##INFO=<ID=REP,Number=1,Type=String,Description="Repetitive element overlapping insertion breakpoint">
##INFO=<ID=CONTIG5,Number=1,Type=String,Description="Assembled contig sequence spanning 5' bkp">
##INFO=<ID=CONTIG3,Number=1,Type=String,Description="Assembled contig sequence spanning 3' bkp">
##INFO=<ID=TRDS,Number=.,Type=String,Description="Reads from the tumour sample (X) that contribute to this insertion">
##FILTER=<ID=SAT,Description="Insertion breakpoint overlapping a satellite region">
##FILTER=<ID=FAM,Description="Insertion breakpoint overlapping a repetive element of the same family">
##FILTER=<ID=SCORE,Description="Insertion with an score > 2">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Count of countributing reads">
##SAMPLE=<ID=NORMAL,Description="Normal",Accession=.,Platform=ILLUMINA,Protocol=WGS,SampleName=X,Source=.>
##SAMPLE=<ID=TUMOUR,Description="Tumour",Accession=.,Platform=ILLUMINA,Protocol=WGS,SampleName=X,Source=.>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    	
""" 
	## 3. Replace variables into the template and save into an output file
	with  open(outFilePath,'w') as outFile:
	    outFile.write(template.format(**context))


    def print_lines(self, outFilePath):
	""" 
	    ...
            
	    Input:
	    1) ...
            
	    Output:
	    1) ...
	"""

	## For each VCF line (each line contains information about an insertion)
	for VCFline in self.lineList:
	   print VCFline.chrom, VCFline.pos, VCFline.id, VCFline.ref, VCFline.alt, VCFline.qual, VCFline.filter, VCFline.info
	 

class VCFline():
    """ 
	.....................
    
	Methods:
	- 

    """

    def __init__(self, insertionObj):
	""" 
	    ...
            
	    Output:
            - 
	"""
	self.chrom = insertionObj.bkpA[0]
	self.pos = insertionObj.bkpA[1]
	self.id = "."
	self.ref = "X"
	self.alt = "<MEI>"
	self.qual = "."
	self.filter = "."
	self.info = self.make_info(insertionObj)
	
    def make_info(self, insertionObj):
	""" 
	    ...
            
	    Output:
            - 
	"""
	
	## 
	infoOrder = [ "SVTYPE", "CLASS", "TYPE", "SCORE", "CIPOS", "STRAND", "STRUCT", "LEN", "TSLEN", "TSDSEQ", "POLYA", "REGION", "GENE", "TID", "SAT", "REP", "CONTIG5", "CONTIG3", "TRDS" ] 

	##
	infoDict = {}
	infoDict["SVTYPE"] = "<MEI>"
	infoDict["CLASS"] = insertionObj.TEClass
	infoDict["TYPE"] = "TD0"
	infoDict["SCORE"] = insertionObj.score
	infoDict["CIPOS"] = insertionObj.bkpA[2]
	infoDict["STRAND"] = insertionObj.orientation
	infoDict["STRUCT"] = insertionObj.structure
	infoDict["LEN"] = insertionObj.length
	infoDict["TSLEN"] = insertionObj.targetSiteSize
	infoDict["TSDSEQ"] = insertionObj.targetSiteSeq
	infoDict["POLYA"] = insertionObj.polyA
	infoDict["REGION"] = "unkn"
	infoDict["GENE"] = "unkn"
	infoDict["TID"] = "unkn"
	infoDict["SAT"] = "unkn"
	infoDict["REP"] = "unkn"
	infoDict["CONTIG5"] = insertionObj.informativeContigBkpA
	infoDict["CONTIG3"] = insertionObj.informativeContigBkpB
	infoDict["TRDS"] = "unkn"

	##
	infoList = []

	for info in infoOrder:
		
	    if (infoDict[info] != "unkn"):
  
		infoField = info + "=" +str(infoDict[info])  
		infoList.append(infoField)
	    
	info = ';'.join(infoList)
	  
	return(info)


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
        self.TEClass = TEClass
        self.coordinates = coordinates
        self.clusterPlusObj = self.create_cluster("+", contigsPlusPath, blatPlusPath)
        self.clusterMinusObj = self.create_cluster("-", contigsMinusPath, blatMinusPath)

	## Unknown values by default:
	self.traficId = "unkn"
	self.score = "unkn"
	self.bkpA = ["unkn", "unkn", "unkn"] 
	self.bkpB = ["unkn", "unkn", "unkn"] 
	self.targetSiteSize = "unkn"
	self.targetSiteSeq = "unkn"	
	self.orientation = "unkn"
	self.structure = "unkn"
	self.length = "unkn"
	self.percLength = "unkn"
	self.informativeContigIdBkpA = "unkn"
	self.informativeContigBkpA = "unkn"
	self.informativeContigIdBkpB = "unkn"
	self.informativeContigBkpB = "unkn"
	self.polyA = "unkn"
 
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
            Identify TE insertion breakpoints, TSD, orientation and poly-A sequence from assembled contigs 
            
            Input:
            1) outDir. Output directory (provisional).
            
            Output (Provisional):
            1) score. One of these values:
                1.  5' and 3' informative contigs 
                2A. 5' informative contig
                2B. 3' informative contig 
		3. No informative contig identified.
                4. Inconsistent insertion. Contradictory orientation (5' and 3' informative contigs suggest opposite orientation) or breakpoint/s () 
            2) breakpoint. Breakpoint coordinates (3 elements list: chrom, pos and confidence interval (CI))
            3) TS. Target site duplication or microdeletion (tuple: TSD size and sequence). 
            4) orientation. TE insertion DNA strand/orientation (+ or -) 
            5) polyA. Poly-A sequence. 
        """
        
        ## Find informative contig for + cluster 
        subHeader("Searching informative contigs for + cluster")
	
	bestInformative5primeContigPlusObj, bestInformative3primeContigPlusObj = self.clusterPlusObj.find_informative_contigs(self.coordinates)

        ## Find informative contig for - cluster
        subHeader("Searching informative contigs for - cluster")

        bestInformative5primeContigMinusObj, bestInformative3primeContigMinusObj = self.clusterMinusObj.find_informative_contigs(self.coordinates)
        
	### Determine insertion breakpoints, TS and TE orientation from informative contigs 
        subHeader("Determining insertion breakpoint, TS and TE orientation from informative contigs")
        
	## Set default variables         
        self.traficId = self.TEClass + ":" + self.coordinates
        
	## A) Insertion without any contig spanning one of the insertion breakpoints
	if (bestInformative5primeContigPlusObj == "unkn") and (bestInformative3primeContigPlusObj == "unkn") and (bestInformative5primeContigMinusObj == "unkn") and (bestInformative3primeContigMinusObj == "unkn"): 

	    info("no-informative-contigs:")
	    
	    # No informative contigs -> Use TraFiC + and - cluster coordinates
	    insertionCoordList = self.coordinates.split("_")
	    chrom = str(insertionCoordList[0])
            beg = int(insertionCoordList[1])
            end = int(insertionCoordList[2])
	    
	    # Output variables
	    self.score = '3'
	    self.bkpA = [chrom, beg, "unkn"]
	    self.bkpB = [chrom, end, "unkn"]

	## B) Insertion with at least one contig spanning one of the insertion breakpoints
	else:

	    info("informative-contig:")

	    ## 1. Determine 5' informative contig and insertion breakpoint
	    # a) 5' informative contigs in + and - clusters:
	    if (bestInformative5primeContigPlusObj != "unkn") and (bestInformative5primeContigMinusObj != "unkn"):
	        
		## Evaluate breakpoint consistency between both 5' informative contigs (+ and -)
		informative5primeContigObj = bestInformative5primeContigPlusObj	
		
		bkpPos5primePlus = bestInformative5primeContigPlusObj.informativeDict["bkp"][1]
		bkpPos5primeMinus = bestInformative5primeContigMinusObj.informativeDict["bkp"][1]
		
              	# Compute consensus breakpoint position.
		# Consensus defined as mean position + confidence interval (CI):
		#     + bkp                       - bkp	
		#  -----------*               *------------
		#              <-----Mean----->
                #              <-CI-> 
		bkpPos5prime = int(bkpPos5primePlus + bkpPos5primeMinus) / 2

		CI = int(abs(bkpPos5primePlus - bkpPos5primeMinus)) / 2
		bkp5prime = [ informative5primeContigObj.informativeDict["bkp"][0], bkpPos5prime, CI]
	
	    # b) 5' informative contig in + cluster
            elif (bestInformative5primeContigPlusObj != "unkn"):
		informative5primeContigObj = bestInformative5primeContigPlusObj
		bkp5prime = informative5primeContigObj.informativeDict["bkp"] + [0]

	    # c) 5' informative contig in - cluster	    
	    elif (bestInformative5primeContigMinusObj != "unkn"):
		informative5primeContigObj = bestInformative5primeContigMinusObj 
		bkp5prime = informative5primeContigObj.informativeDict["bkp"] + [0] 

            # d) none 5' informative contig
	    else:
         	informative5primeContigObj = "unkn"
		bkp5prime = ["unkn", "unkn", "unkn"]

	    ## 2. Determine 3' informative contig and insertion breakpoint
	    # a) 3' informative contigs in + and - clusters:
	    if (bestInformative3primeContigPlusObj != "unkn") and (bestInformative3primeContigMinusObj != "unkn"):
	        
		## Evaluate breakpoint consistency between both 3' informative contigs (+ and -)
		informative3primeContigObj = bestInformative3primeContigPlusObj	
		
		bkpPos3primePlus = bestInformative3primeContigPlusObj.informativeDict["bkp"][1]
		bkpPos3primeMinus = bestInformative3primeContigMinusObj.informativeDict["bkp"][1]
		
              	# Compute consensus breakpoint position.
		# Consensus defined as mean position + confidence interval (CI):
		#     + bkp                       - bkp	
		#  -----------*               *------------
		#              <-----Mean----->
                #              <-CI-> 
		bkpPos3prime = int(bkpPos3primePlus + bkpPos3primeMinus) / 2

		CI = int(abs(bkpPos3primePlus - bkpPos3primeMinus)) / 2
		bkp3prime = [ informative3primeContigObj.informativeDict["bkp"][0], bkpPos3prime, CI]

		## Poly-A
		self.polyA = informative3primeContigObj.informativeDict["info"] 

	    # b) 3' informative contig in + cluster
            elif (bestInformative3primeContigPlusObj != "unkn"):
		informative3primeContigObj = bestInformative3primeContigPlusObj
		bkp3prime = informative3primeContigObj.informativeDict["bkp"] + [0]
		self.polyA = informative3primeContigObj.informativeDict["info"] 

	    # c) 3' informative contig in - cluster	    
	    elif (bestInformative3primeContigMinusObj != "unkn"):
		informative3primeContigObj = bestInformative3primeContigMinusObj 
		bkp3prime = informative3primeContigObj.informativeDict["bkp"] + [0] 
		self.polyA = informative3primeContigObj.informativeDict["info"] 

	    # d) none 3' informative contig
	    else:
         	informative3primeContigObj = "unkn"
		bkp3prime = ["unkn", "unkn", "unkn"]
		polyA = "unkn"


	    ## 3. Determine insertion score and Target Site Duplication (TSD) or microdeletion (TSM) if possible:
	    CI5prime = bkp5prime[2]
	    CI3prime = bkp3prime[2]

	    # a) Inconsistent bkp 5'
	    if (CI5prime > 8) and (CI5prime != "unkn"): 
		info("inconsistent bkp 5'")
		self.score = '4'
		self.targetSiteSize = "unkn" 
		self.targetSiteSeq = "unkn"
		self.orientation = "unkn"
		self.structure = "unkn"
		self.length = "unkn"
		self.percLength = "unkn"

	    # b) Inconsistent bkp 3'
	    elif (CI3prime > 8) and (CI3prime != "unkn"): 
		info("inconsistent bkp 3'")
		self.score = '4'
		self.targetSiteSize = "unkn" 
		self.targetSiteSeq = "unkn"
		self.orientation = "unkn"
		self.structure = "unkn"
		self.length = "unkn"
		self.percLength = "unkn"

	    # c) Consistent 5' and 3' bkps/informative_contigs 
	    elif (informative5primeContigObj != "unkn") and (informative3primeContigObj != "unkn"):   	    
		info("5' and 3' informative contigs:")
		self.score = '1'	

		# Find TSD or TSM
            	self.targetSiteSize, self.targetSiteSeq = self.target_site(informative5primeContigObj, informative3primeContigObj)

		# Inconsistent TSD: 
		if (self.targetSiteSize == "inconsistent"):
		    self.score = '4'
		    self.orientation = "unkn"
		    self.structure = "unkn"
		    self.length = "unkn"
        	    self.percLength = "unkn"

	    # d) 5' bkp/informative_contig
	    elif (informative5primeContigObj != "unkn"):
		info("5' informative contig:")
		self.score = '2A'

	    # e) 3' bkp/informative_contig
	    else:
		info("3' informative contig:")
		self.score = '2B'

	    ## 4. Compute TE insertion orientation and structure if not inconsistent
	    if (self.score != '4'):
		    # TE insertion orientation
		    self.orientation = self.insertion_orientation(informative5primeContigObj, informative3primeContigObj)

	    	    # TE insertion structure
	    	    self.structure, self.length, self.percLength = self.insertion_structure(informative5primeContigObj)

		    # Inconsistent orientation: 
		    if (self.orientation == "inconsistent"):
			self.score = '4'
			self.targetSiteSize = "unkn" 
			self.targetSiteSeq = "unkn"
			self.structure = "unkn"
			self.length = "unkn"
			self.percLength = "unkn"

	    ## 5. Order breakpoints by coordinates
	    bkpCoord5prime  = bkp5prime[1]
	    bkpCoord3prime  = bkp3prime[1]
	    
	    # a) 5' bkp characterized
	    if (bkpCoord3prime == "unkn"):
		
		self.informativeContigIdBkpA = informative5primeContigObj.ID
		self.bkpA = bkp5prime
		self.informativeContigBkpA = informative5primeContigObj.seq
	    		
	    # b) 3' bkp characterized
	    elif (bkpCoord5prime == "unkn"):

		self.informativeContigIdBkpA = informative3primeContigObj.ID
	    	self.bkpA = bkp3prime
		self.informativeContigBkpA = informative3primeContigObj.seq
	    	
	    # c) 5' and 3' bkp characterized
	    else:

		# c.a) 5' bkp < 3' bkp
		if (bkpCoord5prime < bkpCoord3prime):
		
		    self.informativeContigIdBkpA = informative5primeContigObj.ID
	    	    self.informativeContigIdBkpB = informative3primeContigObj.ID
   	    	    self.bkpA = bkp5prime
	    	    self.bkpB = bkp3prime
		    self.informativeContigBkpA = informative5primeContigObj.seq
	    	    self.informativeContigBkpB = informative3primeContigObj.seq
	
		# c.b) 3' bkp < 5' bkp
		else:

		    self.informativeContigIdBkpA = informative3primeContigObj.ID
	    	    self.informativeContigIdBkpB = informative5primeContigObj.ID
   	    	    self.bkpA = bkp3prime
	    	    self.bkpB = bkp5prime
		    self.informativeContigBkpA = informative3primeContigObj.seq
	    	    self.informativeContigBkpB = informative5primeContigObj.seq

    def target_site(self, informative5primeContigObj, informative3primeContigObj):
        """ 
            Determine Target Site Duplication (TSD) or Target Site Microdeletion Size (TSM).

	    1) TSD:
                     
            --------------------------- bkpPlus
                    bkpMinus --------------------
                             <-------->
                                 TSD (7bp)
	    2) TSM:
		      bkpPlus
	    -----------------          bkpMinus
                                       --------------------
                             <-------->
                                 TSM (7bp)

            Input:
            1) informative5primeContigObj. 
            2) informative3primeContigObj
                     
            Output:
            1) targetSiteSize. Target site duplication or microdeletion length. 'inconsistent' if expected TSD size different to TSD sequence length.
            2) targetSiteSeq. Target site duplication sequence or 'na' if no TSD or TSM. 'inconsistent' if expected TSD size different to TSD sequence length.
        """
        
	bkpPos5prime = informative5primeContigObj.informativeDict["bkp"][1]
        bkpPos3prime = informative3primeContigObj.informativeDict["bkp"][1]	
	alignObj5prime = informative5primeContigObj.informativeDict["targetRegionAlignObj"]
        
        # A) Target site duplication (TSD)
        if (bkpPos5prime > bkpPos3prime):
            
            ## Compute TSD length
            targetSiteSize = bkpPos5prime - bkpPos3prime
        
            ## Extract TSD sequence
            # A.a) Begin of the contig sequence aligned in the TE insertion genomic region
            #   -------------**TSD**######TE#####
            #   --------------------
            # qBeg               *qEnd*
            if (alignObj5prime.alignType == "beg"):
                beg = alignObj5prime.qEnd - targetSiteSize
                end = alignObj5prime.qEnd
                targetSiteSeq = informative5primeContigObj.seq[beg:end]
        
            # A.b) End of the contig sequence aligned in the TE insertion genomic region
            #   ######TE#####AAAAAAA**TSD**-------------
            #                       --------------------
            #                    *qBeg*               qEnd
            elif (alignObj5prime.alignType == "end"):
                beg = alignObj5prime.qBeg 
                end = alignObj5prime.qBeg + targetSiteSize
                targetSiteSeq = informative5primeContigObj.seq[beg:end]     
	    
	    ## Inconsistent TSD if sequence has not the expected length or TSD > 100 bp or TSD 
	    if (targetSiteSize != len(targetSiteSeq)) or (targetSiteSize > 100):
		targetSiteSize = "inconsistent"
		targetSiteSeq = "inconsistent"

	# B) Target site microdeletion (TSM)        
	elif (bkpPos5prime < bkpPos3prime):    
	    
	    ## Compute TSM length
            targetSiteSize =  bkpPos5prime - bkpPos3prime 
	    targetSiteSeq = "unkn"

	    ## Inconsistent TSM if TSM < -100 bp
	    if (targetSiteSize < -100):
		targetSiteSize = "inconsistent"
		targetSiteSeq = "inconsistent"

        # C) No TSD or TSM
        else:
            targetSiteSize = 0
            targetSiteSeq = "unkn"
    
        return (targetSiteSize, targetSiteSeq)
 
       
    def insertion_orientation(self, informative5primeContigObj, informative3primeContigObj):
        """ 
            Determine TE insertion strand/orientation.
            
            1) + orientation:
                5' informative contig      -------beg-------####TE####
                3' informative contig      AAAAAAAAAA--------end------
            
            2) - orientation (the opposite)
                3' informative contig      --------beg------AAAAAAAAAA 
                5' informative contig      ####TE####-------end-------

            Input:
            1) informative5primeContigObj 
            2) informative3primeContigObj
            
            Output:
            1) orientation. +, -, na or inconsistent dna strand. 

	    Note: Inconsistent when 5' and 3' informative contigs suggest contradictory orientations (should not happen). 
        """
	
	## A) 5' and 3' informative contigs
	if (informative5primeContigObj != "unkn") and (informative3primeContigObj != "unkn"):   	
            alignType5prime = informative5primeContigObj.informativeDict["targetRegionAlignObj"].alignType 
            alignType3prime = informative3primeContigObj.informativeDict["targetRegionAlignObj"].alignType 
        
	    # a) + strand
	    if (alignType5prime == "beg") and (alignType3prime == "end"):	 
		orientation = "+"

            # b) - strand
            elif (alignType5prime == "end") and (alignType3prime == "beg"):
		orientation = "-"

	    # c) Contradictory/inconsistent
            else:   	
		orientation = "inconsistent"

	## B) 5' informative contig 
 	elif (informative5primeContigObj != "unkn"):
	    alignType5prime = informative5primeContigObj.informativeDict["targetRegionAlignObj"].alignType 
            
	    # a) + strand
	    if (alignType5prime == "beg"):
		orientation = "+"
	    
            # b) - strand
	    else:
		orientation = "-"

	## C) 3' informative contig
	elif (informative3primeContigObj != "unkn"):
            alignType3prime = informative3primeContigObj.informativeDict["targetRegionAlignObj"].alignType 	    

	    # a) + strand
	    if (alignType3prime == "end"):
		orientation = "+"
	    
            # b) - strand
	    else:
		orientation = "-"
     
	## D) None informative contig
	else:    
	    orientation = "unkn"

	return orientation
    
    
    def insertion_structure(self, informative5primeContigObj):
        """
            Determine TE (L1, Alu or SVA) insertion structure.

            1) 5'inverted:

	    	TE in + orientation with 5'inversion    ----#######TE######AAAAA----
                                                            <<<<<<>>>>>>>>>>
        	                                            <---->
                                                           inversion 
        	5-prime informative contig              --------- (5'inversion signature: the piece of contig corresponding to the TE 
                                                                   aligns in the opposite DNA strand than the piece of contig aligning in the target region)

	    2) Full length insetion:

                                                     -------#######TE######AAAAA----
            5-prime informative contig                  ----------
	
	    3) 5' truncated insertion:

	                                               --------###TE######AAAAA----
         					 	   <-->
                                                  	 deletion
                5-prime informative contig              ---____--- (5'truncation signature: the piece of contig corresponding to TE 
                                                                    aligns in the body of the TE and not in the 5' extreme) 

            Input:
            1) informative5primeContigObj
            
            Output:
            1) structure. One of 'na', '5'inverted', '5'truncated' or 'full-length'
            2) length. Inserted TE length, 'na' if not available.
            3) percLength. Percentage of TE consensus sequence inserted, 'na' if not available.  
        """

	## A) 5' informative contig
	if (informative5primeContigObj != "unkn"):

	    strand = informative5primeContigObj.informativeDict["info"].strand

	    ## Determine TE insertion structure
            
	    # a) TE inverted in its 5'
	    if (strand == "-"):
		structure = "5'inverted"
		length = "unkn"
		percLength = "unkn"

	    # b) Full length or 5' truncated
	    else:
	        tBeg = informative5primeContigObj.informativeDict["info"].tBeg
                tSize = informative5primeContigObj.informativeDict["info"].tSize

                length = tSize - tBeg
                percLength = float(length) / tSize * 100
            
                # b.a) full length TE insertion                  
            	if (percLength > 95):
                    # Threshold set for L1 (6021 bp length + 30bp polyA, first ~300bp correspond to promoter)
                    # and Alu (282 bp length + 30bp polyA). For SVA we need to put different values (or not...)
                    structure = "full-length"    

		# b.b) 5' truncated 
       		else:
       		    structure = "5'truncated"

	## B) No 5' informative contig
	else:
	    structure = "unkn"
            length = "unkn"
            percLength = "unkn"

	return (structure, length, percLength)

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

    def find_informative_contigs(self, insertionCoord):
        """
            Identify 5' and 3' informative contig belonging to the cluster. Informative 5' and 3' contigs span 5' and 3' insertion breakpoints, respectively.
            
            Input:
            insertionCoord. Region of interest. Format: ${chrom}_${beg}_${end}. 
                               Example: 10_108820680_108820678.
            
            Output:
            1) bestInformative5primeContigObj. Best 5' Informative contig object. 'na' if not found     
	    2) bestInformative3primeContigObj. Best 3' Informative contig object. 'na' if not found
        """
        
	## Initial status -> none informative contig 
	bestInformative5primeContigObj = "unkn"
	bestInformative3primeContigObj = "unkn"

	informative5primeContigObjList = []
	informative3primeContigObjList = []

	## Determine TraFiC target position
	# A) Target position for positive cluster
	if (self.ID == "+"):
	    targetPos = int(insertionCoord.split("_")[1])
	# B) Target position for negative cluster	
	else:
	    targetPos = int(insertionCoord.split("_")[2])
           
        info(str(len(self.contigsDict)) + ' input contigs') 

        ## Iterate over each contig object checking if it is informative or not. 
	# Make two list of 5' and 3' informative contigs 
        for contigId in self.contigsDict:    
            contigObj = self.contigsDict[contigId]
    	    contigObj.cluster = self.ID

            # Check if it is an informative contig
            informative = contigObj.is_informative(insertionCoord)
	
	    # A) Informative contig
            if (informative ==  1):
	
		# A.a) informative 5-prime 		
		if contigObj.informativeDict['type'] == "5-prime":
	                informative5primeContigObjList.append(contigObj)

		# A.b) informative 3-prime 		
		else:   
			informative3primeContigObjList.append(contigObj)
	
		message = str(contigObj) + " " + contigObj.ID + " " + contigObj.informativeDict['type'] + " " + contigObj.seq + " " + str(contigObj.informativeDict["bkp"]) + " " + str(contigObj.informativeDict["info"])
		log("INFORMATIVE", message)               
               
            # B) Not informative contig
            else:    
                message = str(contigObj) + " " + contigObj.ID + " " + contigObj.informativeDict['type'] + " " + contigObj.seq + " " + str(contigObj.informativeDict["bkp"]) + " " + str(contigObj.informativeDict["info"])
                log("NOT-INFORMATIVE", message)

	## Select the best 5' and 3' informative contigs
	# Note: Best defined as contig spanning putative insertion  
	# breakpoint closer to the insertion target region
	
	# 1) informative 5-prime 
	
	bestDist = ""

	for contigObj in  informative5primeContigObjList:
	    
	    bkpCoord = contigObj.informativeDict["bkp"][1]
	    dist = abs(targetPos - bkpCoord)

	    if (bestInformative5primeContigObj == "unkn") or (dist < bestDist):
			
		bestInformative5primeContigObj = contigObj
		bestDist = dist
	
	# 2) informative 3-prime  
	
	bestDist = ""

	for contigObj in  informative3primeContigObjList:
		
            bkpCoord = contigObj.informativeDict["bkp"][1]
	    dist = abs(targetPos - bkpCoord)

	    if (bestInformative3primeContigObj == "unkn") or (dist < bestDist):
			
		bestInformative3primeContigObj = contigObj
		bestDist = dist

	return (bestInformative5primeContigObj, bestInformative3primeContigObj)


class contig():
    """ 
    Transposable element insertion contig class. 
    
    Contig sequence results from the assembly of TraFiC + or - cluster supporting reads 
    
    Methods:
    - is_candidate
    - is_informative
    - is_3prime_bkp
    - is_polyA
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

	# Cluster the contig belongs to
	self.cluster = ""

        # Contig alignment information
        self.alignList = []   # List of blat alignments for this contig
        self.informativeDict = {} # Dictionary key -> value pairs:
        self.informativeDict["type"] = "unkn"       	      # type -> 5-prime, 3-prime or none
        self.informativeDict["bkp"] = "unkn"          	      # bkp -> list: chrom and pos, 'na' for both if type == none)
        self.informativeDict["info"] = "unkn"        	      # info -> 5-prime: aligment object with contig's alignment in TE sequence info; 3-prime: PolyA sequence; none: 'na'
        self.informativeDict["targetRegionAlignObj"] = "unkn"   # targetRegionAlignObj -> alignment object with contig's alignment in the target region info. 
                             
    #### FUNCTIONS ####
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


    def is_candidate(self, insertionCoord, windowSize, maxAlignPerc):
        """
        Check if contig is candidate to be informative about the TE insertion breakpoint. Informative contigs 
        span the insertion breakpoint. 

        Candidate contigs are defined as contigs partially aligning in the TE insertion region. 
	Contigs completely alignment on the TE sequence (L1, Alu or SVA) are not considered as candidates.
         
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
            
	    alignPerc = float(alignment.qEnd - alignment.qBeg) / alignment.qSize * 100
    
	    ## A) Discard contig completely aligning in the TE sequence as informative candidate  
	    # TE: L1, Alu or SVA (SVA need to be included)
            if (( alignment.tName == "L1" ) or ( alignment.tName == "Alu" )) and ( alignPerc > 99 ):
		
		candidate = 0
	        supportingAlignList = []
		break 
	   
	    ## B) Contig do not alignining in the TE sequence
            else:
	
	        # 1. Check if alignment within the target region
	        insertionRegion = alignment.in_target_region(insertionCoord, windowSize)
            
                # Within target region
                if (insertionRegion == 1):
              
                    # 2. Check if it is a partial alignment
                    partial = alignment.partial_alignment(maxAlignPerc)
                
                    # Partial
                    if (partial == 1):            
                        supportingAlignList.append(alignment)

			# Informative candidate contig -> partially alignining on the target region and do not aligning completely on the TE sequence                        
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
            2) Sets 'informativeDict' variable. Dictionary key -> value pairs:
                type -> 5-prime, 3-prime or none
                bkp -> list: chrom and pos, 'na' for both if type == none)
                info -> 5-prime: aligment object with contig's alignment in TE sequence info; 
                        3-prime: PolyA sequence; none: 'na'
                targetRegionAlignObj -> alignment object with contig's alignment in the target region info. 
        """
            
	## Initial status -> no informative	
	bestInformative5primeDict = "unkn"
	bestInformative3primeDict = "unkn"
	informative5primeDict = {}
        informative3primeDict = {}  

	## Determine TraFiC target position
	# A) Target position for positive cluster
	if (self.cluster == "+"):
	    targetPos = int(insertionCoord.split("_")[1])
	# B) Target position for negative cluster	
	else:
	    targetPos = int(insertionCoord.split("_")[2])

        ## Check if it is an informative candidate contig
        candidate, supportingAlignList = self.is_candidate(insertionCoord, int(500), float(98))
           
        ## Informative candidate contig:
        if (candidate == 1):

	    ## Iterate over the alignments supporting the contig as informative. 
            # Check if the alignment support the contig as informative 5' and/or 3'. 
	    # Add alignment to the list of 5' and/or 3' alignments 
            for alignObj in supportingAlignList:  
                
		## Check if alignment support the contig as informative 5' 
                is5prime, bkpCoord, TEalignmentObj = self.is_5prime_informative(alignObj)
                    
                # Contig informative 5-prime, it has TE sequence 
                if (is5prime == 1):
		    
		    informative5primeDict[alignObj] = {}
		    informative5primeDict[alignObj]["type"] = "5-prime"
		    informative5primeDict[alignObj]["bkp"] = bkpCoord
		    informative5primeDict[alignObj]["info"] = TEalignmentObj
		    informative5primeDict[alignObj]["targetRegionAlignObj"] = alignObj

		## Check if alignment support the contig as informative 3' 
                is3prime, bkpCoord, polyASeq = self.is_3prime_informative(alignObj)
                
                # Contig informative 3-prime, it has a polyA tail 
                if (is3prime == 1):

		    informative3primeDict[alignObj] = {}
		    informative3primeDict[alignObj]["type"] = "3-prime"
		    informative3primeDict[alignObj]["bkp"] = bkpCoord
		    informative3primeDict[alignObj]["info"] = polyASeq
		    informative3primeDict[alignObj]["targetRegionAlignObj"] = alignObj
             
	## Select the best alignments supporting the contig as 5' and/or 3' informative contigs
	# Note: Best defined as alignment supporting putative insertion breakpoint closer to the
	# insertion target region
	
	# 1) informative 5-prime 
	bestDist5Prime = ""

	for alignment in  informative5primeDict:
	     bkpCoord = informative5primeDict[alignment]["bkp"][1]	
	     dist = abs(targetPos - bkpCoord)

	     if (bestInformative5primeDict == "unkn") or (dist < bestDist5Prime):
		bestInformative5primeDict = informative5primeDict[alignment]
		bestDist5Prime = dist
	
	# 2) informative 3-prime 
	bestDist3Prime = ""

	for alignment in informative3primeDict:
	     bkpCoord = informative3primeDict[alignment]["bkp"][1]	
	     dist = abs(targetPos - bkpCoord)

	     if (bestInformative3primeDict == "unkn") or (dist < bestDist3Prime):
		bestInformative3primeDict = informative3primeDict[alignment]
		bestDist3Prime = dist

	## In case there are alignments supporting the contig is 5' and 3' select the best
	# Note: Best defined as alignment supporting putative insertion breakpoint closer to the
	# insertion target region

	# A) Any alignment supporting the contig as 5' or 3' informative
	if (bestDist5Prime == "") and (bestDist3Prime == ""):
	    informativeBoolean = 0

	# B) Best 5' 
	elif (bestDist5Prime <= bestDist3Prime) or ((bestDist5Prime != "") and (bestDist3Prime == "")):
	    informativeBoolean = 1
	    self.informativeDict = bestInformative5primeDict

	# C) Best 3'
	else:
	    informativeBoolean = 1
	    self.informativeDict = bestInformative3primeDict

	## Fix contig orientation if necessary (contig aligning in -)
	if (informativeBoolean == 1):

	    ## Contig aligning in - strand
	    if (self.informativeDict["targetRegionAlignObj"].strand == "-"):

		## Substitute contig by its reverse complementary sequence		
		self.seq = self.rev_complement(self.seq)		
		
		## Fix contig alignment coordinates on the TE insertion genomic region 
		switchTypeDict = {'beg': 'end', 'end': 'beg'}
		alignType = self.informativeDict["targetRegionAlignObj"].alignType
		self.informativeDict["targetRegionAlignObj"].alignType = switchTypeDict[alignType] 
		self.informativeDict["targetRegionAlignObj"].rev_complement()

		## A) Informative 5-prime -> fix contig alignment coordinates on the transposon sequence
		if (self.informativeDict["type"] == "5-prime"):
	
		    self.informativeDict["info"].rev_complement()

		## B) Informative 3-prime -> make the complementary reverse of the poly-A sequence
		else:
		    self.informativeDict["info"] = self.rev_complement(self.informativeDict["info"])
	    
	return informativeBoolean 
	      
                
    def is_3prime_informative(self, alignObj):
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
            polyASeq = self.is_polyA(targetSeq)
        
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
                
            ## Check if the contig target piece of sequence corresponds to polyA 
            polyASeq = self.is_polyA(targetSeq)
            
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
            bkpCoord = ["unkn", "unkn"]
            polyASeq = "unkn"
            
        return (is3prime, bkpCoord, polyASeq)
            
        
    def is_polyA(self, targetSeq):
        """
        Check if input sequence is a poly-A tail. 
        
        Compute the percentage of T and A for the input sequence. 
        If it is >80 classify sequence as poly-A.
        
        Input:
        1) targetSeq. Target DNA sequence to search for poly A. 
        
        Output:
        1) polyASeq. PolyA sequence.            
        """
        
	polyASeq = ""

	# Convert sequence into upper case:
        targetSeq = targetSeq.upper()
        
	# Compute percentage of A and T in the given slice
        nbA = targetSeq.count("A") 
        nbG = targetSeq.count("G") 
        nbC = targetSeq.count("C") 
        nbT = targetSeq.count("T") 
        nbN = targetSeq.count("N")
       
        percA = float(nbA) / (nbA + nbG + nbC + nbT + nbN) * 100
        percT = float(nbT) / (nbA + nbG + nbC + nbT + nbN) * 100
            
        # Classify the input sequence as poly-A if > 80% are T or A
        if (percA >= 80) or (percT >= 80):    
            polyASeq = targetSeq
                
        return polyASeq
      
    def is_5prime_informative(self, alignObj):
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
        
        ## Select contig target sequence coordinates to search for alignment in TE sequence (L1, Alu or SVA)
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
        bkpCoord = ["unkn", "unkn"]
        TEalignmentObj = "unkn" 
                
        for alignment in self.alignList:
            
            # Contig alignment in TE sequence (L1, Alu or SVA (SVA need to be included))
            if ( alignment.tName == "L1" ) or ( alignment.tName == "Alu" ):
                
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
parser.add_argument('inputPaths', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('genome', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputPaths = args.inputPaths
donorId = args.donorId
genome = args.genome
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "paths2bkpAnalysis: ", inputPaths
print "donorId: ", donorId
print "genome: ", genome
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, " *****"
print 
print "..."
print 


## Start ## 

outFilePath = outDir + '/' + donorId + '.vcf'

## 1. Create VCF output file and print VCF header

VCFObj = VCF()

VCFObj.print_header(outFilePath)

## 2. Per each insertion perform breakpoint analysis 

inputFile = open(inputPaths, 'r')

# Analyze one insertion per iteration
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")
    
    # Get TE insertion info and files
    TEClass, insertionCoord = line[0].split(":")
    contigsPlusPath, contigsMinusPath = line[1].split(",")
    blatPlusPath, blatMinusPath = line[2].split(",")

    # Perform breakpoint analysis for the TE insertion 
    header("Tranposable Element Insertion Breakpoint Analysis (TEIBA) for: " + insertionCoord)
    
    # A) All the input files exist 
    if os.path.isfile(contigsPlusPath) and os.path.isfile(blatPlusPath) and os.path.isfile(contigsMinusPath) and os.path.isfile(blatMinusPath):  

	## Create insertion object and identify breakpoints from assembled contigs        
	insertionObj = insertion(TEClass, insertionCoord, contigsPlusPath, blatPlusPath, contigsMinusPath, blatMinusPath)
        insertionObj.find_insertionBkp(outDir)

	## Create VCF line object 
	VCFlineObj = VCFline(insertionObj)
	VCFObj.addLine(VCFlineObj)

    else:
        message = "Input files for " + insertionCoord + " insertion do not exist"
        log("ERROR", message)


print 
print 
print "****************************"

print VCFObj.print_lines(outFilePath)


## Finish ##
print 
print "***** Finished! *****"
print 






