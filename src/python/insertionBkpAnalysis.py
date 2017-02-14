#!/usr/bin/env python
#coding: utf-8

## Import modules ##
import argparse
import time
import sys
from itertools import groupby
import os.path
import re

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


    def sort(self):
        """
            ...

            Input:
            1) ...

            Output:
            1) ...
        """

        lineListSorted = sorted(self.lineList, key=lambda line: (line.chrom, line.pos))

        return lineListSorted

    def print_header(self, outFilePath):
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
         "reference": "hs37d5"
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
##INFO=<ID=CLASS,Number=1,Type=String,Description="Transposable element class (L1, ALU, SVA or ERVK)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Insertion type (TD0: solo, TD1: partnered-3'transduction, TD2: orphan-3'transduction), PSD: processed-pseudogene">
##INFO=<ID=SCORE,Number=1,Type=String,Description="Insertion score (5: 5' and 3' breakpoints (bkp) assembled, 4: 3'bkp assembled, 3: 5'bkp assembled, 2: no bkp assembled, 1: inconsistent (contradictory orientation, bkp or TSD))">
##INFO=<ID=BKPB,Number=1,Type=String,Description="MEI right-most breakpoint position (bkp B). Left-most breakpoint position (bkp A) represented in the POS field">
##INFO=<ID=CIPOS,Number=1,Type=Integer,Description="Confidence interval around insertion breakpoints">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion DNA strand (+ or -)">
##INFO=<ID=STRUCT,Number=1,Type=String,Description="Transposable element structure (INV: 5'inverted, DEL: 5'deleted, FULL: full-length)">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Transposable element length">
##INFO=<ID=TSLEN,Number=1,Type=Integer,Description="Target site duplication length">
##INFO=<ID=TSSEQ,Number=1,Type=String,Description="Target site duplication sequence">
##INFO=<ID=POLYA,Number=1,Type=String,Description="Poly-A sequence">
##INFO=<ID=SRC,Number=1,Type=String,Description="Coordinates of the source element that mediated the transduction in the format: chrom_beg_end_strand">
##INFO=<ID=TDC,Number=1,Type=String,Description="Begin and end coordinates of the transduced region in the format: chrom_beg_end">
##INFO=<ID=TDLEN,Number=1,Type=String,Description="Transduced region length">
##INFO=<ID=TDLENR,Number=1,Type=String,Description="Transduced region length at RNA level">
##INFO=<ID=PSDGENE,Number=1,Type=String,Description="Source gene of the processed pseudogene insertion">
##INFO=<ID=GERMDB,Number=1,Type=String,Description="MEI already reported as germinal in a database (1KGENOMES: 1000 genomes project (source_papers_doi: 10.1038/nature15394 and 10.1073/pnas.1602336113), TRAFIC: TraFic in-house database)">
##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region where the transposable element is inserted (exonic, splicing, ncRNA, UTR5, UTR3, intronic, upstream, downstream, intergenic)">
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO gene symbol">
##INFO=<ID=ROLE,Number=1,Type=String,Description="Role in cancer (oncogene, TSG: tumor suppressor gene, oncogene/TSG: both roles)">
##INFO=<ID=COSMIC,Number=0,Type=Flag,Description="Reported as cancer driver in COSMIC cancer gene census database">
##INFO=<ID=CPG,Number=0,Type=Flag,Description="Reported as cancer predisposition gene in 10.1038/nature12981 (DOI).">
##INFO=<ID=REP,Number=1,Type=String,Description="Repetitive element overlapping the insertion breakpoint">
##INFO=<ID=DIV,Number=1,Type=Integer,Description="Millidivergence of the overlapping repetitive element with respect a consensus sequence">
##INFO=<ID=MEISEQ,Number=1,Type=String,Description="Assembled sequence of the mobile element 5' boundary">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Assembled contig sequence spanning 1st bkp (lowest genomic position)">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Assembled contig sequence spanning 2nd bkp (highest genomic position)">
##INFO=<ID=RP,Number=.,Type=String,Description="Reads from the tumour sample and positive cluster that support this insertion">
##INFO=<ID=RN,Number=.,Type=String,Description="Reads from the tumour sample and negative cluster that support this insertion">
##FILTER=<ID=SCORE,Description="Insertion with an score < threshold">
##FILTER=<ID=REP,Description="Insertion overlapping a satellite region or a repetitive element of the same class">
##FILTER=<ID=DUP,Description="Duplicated MEI call">
##FILTER=<ID=GERMLINE,Description="Germline MEI miscalled as somatic">
##FORMAT=<ID=RCP,Number=1,Type=Integer,Description="Count of positive cluster supporting reads">
##FORMAT=<ID=RCN,Number=1,Type=Integer,Description="Count of negative cluster supporting reads">
##FORMAT=<ID=SL,Number=1,Type=Integer,Description="List of samples where the variant was found (relevant for multi-tumor donors)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting the variant in this sample">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads covering variant location in this sample">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  
"""
        ## Replace variables into the template and print header into the output file
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

        ## Open outfile
        outFile = open(outFilePath, 'a')

        ## Sort VCF lines
        lineListSorted = self.sort()

        ## Iterate and print each VCF line into the output VCF file
        for VCFline in lineListSorted:

            row = VCFline.chrom + "\t" + str(VCFline.pos) + "\t" + VCFline.id + "\t" + VCFline.ref + "\t" + VCFline.alt + "\t" + VCFline.qual + "\t" + VCFline.filter + "\t" + VCFline.info + "\t" + VCFline.format + "\t" + VCFline.genoType + "\n"
            outFile.write(row)

        ## Close output file
        outFile.close()


class VCFline():
    """
        .....................

        Methods:
        -

    """

    def __init__(self, insertionObj, genomeObj):
        """
            ...

            Output:
            -
        """
        self.chrom = insertionObj.bkpA[0]
        self.pos = insertionObj.bkpA[1]
        self.id = "."
        self.ref = genomeObj.fastaDict[self.chrom][insertionObj.bkpA[1] - 1]  # Substract 1 since python string coordinates start in 0 while bkp position in 1.
        self.alt = "<MEI>"
        self.qual = "."
        self.filter = "."
        self.info = self.make_info(insertionObj)
        self.format = "RCP:RCN:SL"
        self.genoType = str(insertionObj.clusterPlusObj.nbPairs) + ":" + str(insertionObj.clusterMinusObj.nbPairs) + ":" + insertionObj.sampleId

    def make_info(self, insertionObj):
        """
            ...

            Output:
            -
        """

        ## Create list containing the order of info fields
        infoOrder = [ "SVTYPE", "CLASS", "TYPE", "SCORE", "BKPB", "CIPOS", "STRAND", "STRUCT", "LEN", "TSLEN", "TSSEQ", "POLYA", "SRC", "TDC", "TDLEN", "TDLENR", "PSDGENE", "GERMDB", "REGION", "GENE", "ROLE", "COSMIC", "CPG", "REP", "DIV", "MEISEQ", "CONTIGA", "CONTIGB", "RP", "RN" ]

        ## Build dictionary with info tags as keys
        infoDict = {}
        infoDict["SVTYPE"] = "<MEI>"
        infoDict["CLASS"] = insertionObj.TEClass
        infoDict["TYPE"] = insertionObj.tdType
        infoDict["SCORE"] = insertionObj.score
        infoDict["BKPB"] = insertionObj.bkpB[1]
        infoDict["CIPOS"] = insertionObj.bkpA[2]
        infoDict["STRAND"] = insertionObj.orientation
        infoDict["STRUCT"] = insertionObj.structure
        infoDict["LEN"] = insertionObj.length
        infoDict["TSLEN"] = insertionObj.targetSiteSize
        infoDict["TSSEQ"] = insertionObj.targetSiteSeq
        infoDict["POLYA"] = insertionObj.polyA
        infoDict["SRC"] = insertionObj.srcElement
        infoDict["TDC"] = insertionObj.srcCoord
        infoDict["TDLEN"] = insertionObj.tdLen
        infoDict["TDLENR"] = insertionObj.tdLenRna
        infoDict["PSDGENE"] = insertionObj.psdGene
        infoDict["MEISEQ"] = insertionObj.MEISeq
        infoDict["CONTIGA"] = insertionObj.informativeContigBkpA
        infoDict["CONTIGB"] = insertionObj.informativeContigBkpB
        infoDict["RP"] = insertionObj.clusterPlusObj.readPairIds
        infoDict["RN"] = insertionObj.clusterMinusObj.readPairIds

        # handle specifics of pseudogene insertions
        if (insertionObj.tdType == "PSD"):
            infoDict["TDC"] = "UNK"
            infoDict["MEISEQ"] = "UNK"

        ## Create info string in the correct order from dictionary
        infoList = []

        for info in infoOrder:

            # Information about current info field available and applicable
            if (info in infoDict.keys()) and (infoDict[info] != "UNK") and (infoDict[info] != "NA"):

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
    - insertion_orientation_polyA
    - insertion_orientation_ERVK
    - polyA
    - target_site
    - imprecise_bkp
    """

    def __init__(self, family, tdType, coordinates, contigsPlusPath, blatPlusPath, contigsMinusPath, blatMinusPath, readPairsPlus, readPairsMinus, srcElement, transductionInfo, pseudogeneInfo, sampleId):
        """
            Initialize insertion object.

            Input:
            1) family. TE family (L1, Alu, SVA or ERVK)
            2) tdType. Insertion type:  td0 (solo-insertion), td1 (partnered-transduccion), td2 (orphan-transduction) or psd (pseudogene insertion).
            3) coordinates. TraFiC insertion range.
            4) contigsPlusPath. Fasta file containing the assembled contigs for the positive cluster.
            5) blatPlusPath. psl file containing the blat aligments for the positive cluster's assembled contigs.
            6) contigsMinusPath. Fasta file containing the assembled contigs for the negative cluster.
            7) blatMinusPath. psl file containing the blat aligments for the negative cluster's assembled contigs.
            8) readPairsPlus. List of + cluster supporting reads.
            9) readPairsMinus. List of - cluster supporting reads.
            10) srcElement.
            11) transductionInfo.
            12) pseudogeneInfo.
            13) sampleId
            
            Output:
            - Insertion object variables initialized
        """
        self.TEClass = TEClass
        self.tdType = tdType
        self.coordinates = coordinates
        self.clusterPlusObj = self.create_cluster("+", contigsPlusPath, blatPlusPath, readPairsPlus)
        self.clusterMinusObj = self.create_cluster("-", contigsMinusPath, blatMinusPath, readPairsMinus)
        self.srcElement = srcElement
        self.sampleId = sampleId

        # A) Solo insertion (TD0). Not applicable
        if (self.tdType == "TD0"):

            self.srcCoord = "NA"
            self.tdLen = "NA"
            self.tdLenRna = "NA"
            self.psdGene = "UNK"

        # B) Partnered or orphan transduction (TD1 and TD2)
        elif (self.tdType in ["TD1", "TD2"]):
            
            transductionInfoList = transductionInfo.split("_")
            self.srcCoord = transductionInfoList[0] + "_" + transductionInfoList[1] + "_" + transductionInfoList[2]
            self.psdGene = "UNK"

            status = self.srcElement.split("_")[3]

            print "test_status: ", status

            # A) Putative. Uncharacterized germline or somatic source element
            if (status == "putative"):
                self.tdLenRna = "UNK"                
                self.tdLen = "UNK"
                

            # B) Characterized germline source element
            else:
                self.tdLenRna = transductionInfoList[3]                
                self.tdLen = transductionInfoList[4]
                
        # C) pseudogene insertion (PSD)
        elif (self.tdType == "PSD"):
            self.tdLen = "UNK"
            self.tdLenRna = "UNK"

            # parse pseudogene info, format:
            # gene:chrA_begA-endA:chrB_begB-begB
            regex = r'(?P<gene>\w+):(?P<chrA>\w+)_(?P<begA>\d+)-(?P<endA>\d+):(?P<chrB>\w+)_(?P<begB>\d+)-(?P<endB>\d+)'
            m = re.search(regex, pseudogeneInfo)
            psdGene = m.group("gene")
            chrExonA = m.group("chrA")
            begExonA = int(m.group("begA"))
            endExonA = int(m.group("endA"))
            chrExonB = m.group("chrB")
            begExonB = int(m.group("begB"))
            endExonB = int(m.group("endB"))

            self.psdGene = psdGene

            # construct source region from exon coordinates
            minBeg = min(begExonA, begExonB)
            maxEnd = max(endExonA, endExonB)
            self.srcCoord = "%s_%d_%d" % (chrExonA, minBeg, maxEnd)

        ## Unkown values by default:
        self.traficId = "UNK"
        self.score = "UNK"
        self.bkpA = ["UNK", "UNK", "UNK"]
        self.bkpB = ["UNK", "UNK", "UNK"]
        self.targetSiteSize = "UNK"
        self.targetSiteSeq = "UNK"
        self.orientation = "UNK"
        self.structure = "UNK"
        self.length = "UNK"
        self.percLength = "UNK"
        self.MEISeq = "UNK"
        self.informativeContigIdBkpA = "UNK"
        self.informativeContigBkpA = "UNK"
        self.informativeContigIdBkpB = "UNK"
        self.informativeContigBkpB = "UNK"
        self.polyA = "UNK"

    #### FUNCTIONS ####
    def create_cluster(self, ID, contigsPath, blatPath, readPairsList):
        """
            Create cluster object.

            Input:
            1) ID. Cluster id (+ or -)
            2) contigsPath. Fasta file containing the assembled contigs for the given
                            cluster.
            3) blatPath. psl file containing the blat aligments for the assembled contigs.
            4) readPairsList. List of cluster supporting reads.

            Output:
            1) clusterObj. Cluster object
        """

        # Create cluster object
        clusterObj = cluster(ID, contigsPath, readPairsList)

        # Add blat alignments to cluster object
        clusterObj.add_alignments(blatPath)

        return clusterObj


    def find_insertionBkp(self, genomeObj, outDir):
        """
            Identify TE insertion breakpoints, TSD, orientation and poly-A sequence from assembled contigs

            Input:
            1) outDir. Output directory (provisional).

            Output (Provisional):
            1) score. One of these values:
                5: 5' and 3' breakpoints (bkp) assembled
                4: 3'bkp assembled
                3: 5'bkp assembled
                2: no bkp assembled
                1: inconsistent (contradictory orientation, bkp or TSD)

            2) breakpoint. Breakpoint coordinates (3 elements list: chrom, pos and confidence interval (CI))
            3) TS. Target site duplication or deletion (tuple: TSD size and sequence).
            4) orientation. TE insertion DNA strand/orientation (+ or -)
            5) polyA. Poly-A sequence.
        """

        ## Find informative contig for + cluster 
        subHeader("Searching informative contigs for + cluster")

        bestInformative5primeContigPlusObj, bestInformative3primeContigPlusObj = self.clusterPlusObj.find_informative_contigs(self.coordinates, self.tdType, self.srcCoord)

        ## Find informative contig for - cluster
        subHeader("Searching informative contigs for - cluster")

        bestInformative5primeContigMinusObj, bestInformative3primeContigMinusObj = self.clusterMinusObj.find_informative_contigs(self.coordinates, self.tdType, self.srcCoord)

        ### Determine insertion breakpoints, TS and TE orientation from informative contigs 
        subHeader("Determining insertion breakpoint, TS and TE orientation from informative contigs")

        ## Set default variables
        self.traficId = self.TEClass + ":" + self.coordinates

        ## A) Insertion without any contig spanning one of the insertion breakpoints
        if (bestInformative5primeContigPlusObj == "UNK") and (bestInformative3primeContigPlusObj == "UNK") and (bestInformative5primeContigMinusObj == "UNK") and (bestInformative3primeContigMinusObj == "UNK"):

            info("no-informative-contigs:")

            # No informative contigs, imprecise breakpoint
            self.score = '2'
            self.bkpA = self.imprecise_bkp(self.coordinates)

        ## B) Insertion with at least one contig spanning one of the insertion breakpoints
        else:

            info("informative-contig:")

            ## There are two main classes of MEI events:
            #   1. ERVK-like (insertion without poly-A tail)
            #   2. L1-like (insertion with poly-A tail)

            # 1. ERVK: two 5' informative contigs are detected,
            #          but one represents 3' end of insertion
            if self.TEClass == "ERVK":
                # sanity check: we are not expecting 3' informative contigs here
                if (bestInformative3primeContigPlusObj != "UNK"):
                    log("WARNING", "ERVK event \"%s\" has 3'-informative contig (+ cluster), which is not expected." % self.traficId)
                if (bestInformative3primeContigMinusObj != "UNK"):
                    log("WARNING", "ERVK event \"%s\" has 3'-informative contig (- cluster), which is not expected." % self.traficId)

                # treat 5' informative contig of minus cluster as 3' informative
                bestInformative3primeContigMinusObj = bestInformative5primeContigMinusObj
                bestInformative5primeContigMinusObj = "UNK"

            ## 1. Determine 5' informative contig and insertion breakpoint
            # a) 5' informative contigs in + and - clusters:
            if (bestInformative5primeContigPlusObj != "UNK") and (bestInformative5primeContigMinusObj != "UNK"):

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
            elif (bestInformative5primeContigPlusObj != "UNK"):
                informative5primeContigObj = bestInformative5primeContigPlusObj
                bkp5prime = informative5primeContigObj.informativeDict["bkp"] + [0]

            # c) 5' informative contig in - cluster
            elif (bestInformative5primeContigMinusObj != "UNK"):
                informative5primeContigObj = bestInformative5primeContigMinusObj
                bkp5prime = informative5primeContigObj.informativeDict["bkp"] + [0]

            # d) none 5' informative contig
            else:
                informative5primeContigObj = "UNK"
                bkp5prime = ["UNK", "UNK", "UNK"]

            ## 2. Determine 3' informative contig and insertion breakpoint
            # a) 3' informative contigs in + and - clusters:
            if (bestInformative3primeContigPlusObj != "UNK") and (bestInformative3primeContigMinusObj != "UNK"):

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
            elif (bestInformative3primeContigPlusObj != "UNK"):
                informative3primeContigObj = bestInformative3primeContigPlusObj
                bkp3prime = informative3primeContigObj.informativeDict["bkp"] + [0]
                self.polyA = informative3primeContigObj.informativeDict["info"]

            # c) 3' informative contig in - cluster
            elif (bestInformative3primeContigMinusObj != "UNK"):
                informative3primeContigObj = bestInformative3primeContigMinusObj
                bkp3prime = informative3primeContigObj.informativeDict["bkp"] + [0]
                self.polyA = informative3primeContigObj.informativeDict["info"]

            # d) none 3' informative contig
            else:
                informative3primeContigObj = "UNK"
                bkp3prime = ["UNK", "UNK", "UNK"]
                polyA = "UNK"

            # ERVK-like events do not have poly-A signal
            if self.TEClass == "ERVK":
                self.polyA = "UNK"


            ## 3. Determine insertion score and Target Site Duplication if possible:
            CI5prime = bkp5prime[2]
            CI3prime = bkp3prime[2]

            # a) Inconsistent bkp 5'
            if (CI5prime > 8) and (CI5prime != "UNK"):
                info("inconsistent bkp 5'")
                self.score = '1'
                self.targetSiteSize = "UNK"
                self.targetSiteSeq = "UNK"
                self.orientation = "UNK"
                self.structure = "UNK"
                self.length = "UNK"
                self.percLength = "UNK"

            # b) Inconsistent bkp 3'
            elif (CI3prime > 8) and (CI3prime != "UNK"):
                info("inconsistent bkp 3'")
                self.score = '1'
                self.targetSiteSize = "UNK"
                self.targetSiteSeq = "UNK"
                self.orientation = "UNK"
                self.structure = "UNK"
                self.length = "UNK"
                self.percLength = "UNK"

            # c) Consistent 5' and 3' bkps/informative_contigs
            elif (informative5primeContigObj != "UNK") and (informative3primeContigObj != "UNK"):
                info("5' and 3' informative contigs:")
                self.score = '5'

                # Find Target site Duplication
                self.targetSiteSize, self.targetSiteSeq = self.target_site(informative5primeContigObj, informative3primeContigObj)

                # Inconsistent TSD:
                if (self.targetSiteSize == "inconsistent"):
                    self.score = '1'
                    self.orientation = "UNK"
                    self.structure = "UNK"
                    self.length = "UNK"
                    self.percLength = "UNK"

            # d) 5' bkp/informative_contig
            elif (informative5primeContigObj != "UNK"):
                info("5' informative contig:")
                self.score = '3'

            # e) 3' bkp/informative_contig
            else:
                info("3' informative contig:")
                # ERVK-like events do not have poly-A signal, so they do not receive score 4
                if (self.TEClass != "ERVK"):
                    self.score = '4'
                else:
                    self.score = '3'

            ## 4. Compute TE insertion orientation and structure if not inconsistent
            if (self.score != '1'):

                if (self.TEClass != "ERVK"):
                    # TE insertion orientation
                    self.orientation = self.insertion_orientation_polyA(informative5primeContigObj, informative3primeContigObj)

                    if (self.tdType != "PSD"):
                        # TE insertion structure
                        self.structure, self.length, self.percLength = self.insertion_structure(informative5primeContigObj)
                else:
                    # TE insertion orientation
                    self.orientation = self.insertion_orientation_ERVK(informative5primeContigObj, informative3primeContigObj)
                    #self.structure = "UNK"

                # Inconsistent orientation:
                if (self.orientation == "inconsistent"):
                    self.score = '1'
                    self.targetSiteSize = "UNK"
                    self.targetSiteSeq = "UNK"
                    self.structure = "UNK"
                    self.length = "UNK"
                    self.percLength = "UNK"

            ## 5. Order breakpoints by coordinates
            bkpCoord5prime  = bkp5prime[1]
            bkpCoord3prime  = bkp3prime[1]

            # a) 5' bkp characterized
            if (bkpCoord3prime == "UNK"):

                self.informativeContigIdBkpA = informative5primeContigObj.ID
                self.bkpA = bkp5prime
                self.informativeContigBkpA = informative5primeContigObj.seq

                 # MEI 5' boundary assembled sequence
                self.MEISeq = informative5primeContigObj.informativeDict["MEISeq"]

            # b) 3' bkp characterized
            elif (bkpCoord5prime == "UNK"):

                self.informativeContigIdBkpA = informative3primeContigObj.ID
                self.bkpA = bkp3prime
                self.informativeContigBkpA = informative3primeContigObj.seq

            # c) 5' and 3' bkp characterized
            else:

                # MEI 5' boundary assembled sequence
                self.MEISeq = informative5primeContigObj.informativeDict["MEISeq"]

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

        ## Print results into the standard output
        print "TraFiC-id: ", self.traficId
        print "Score: ", self.score
        print "bkpA: ", self.bkpA
        print "bkpB", self.bkpB
        print "TS-length: ", self.targetSiteSize
        print "TS-seq: ", self.targetSiteSeq
        print "Orientation: ", self.orientation
        print "Structure: ", self.structure
        print "TE-length: ", self.length
        print "perc-Length: ", self.percLength
        print "bkpAContigId: ", self.informativeContigIdBkpA
        print "bkpAContig: ", self.informativeContigBkpA
        print "bkpBContigId: ", self.informativeContigIdBkpB
        print "bkpBContig: ", self.informativeContigBkpB
        print "poly-A: ", self.polyA

    def target_site(self, informative5primeContigObj, informative3primeContigObj):
        """
            Determine Target Site Duplication (TSD) (or microdeletion??).


            *** Target Site Duplication *** 

            + strand)

            --------------------------- bkp5'
                       bkp3' --------------------
                             <-------->
                             TSD (7bp)

            - strand)

            --------------------------- bkp3'
                       bkp5' --------------------
                             <-------->
                             TSD (7bp)


            *** ¿Target Site deletion? ***

            + strand)

            --------------------------- bkp5'
                                                                bkp3' --------------------
                                            <------------------>
                                              deletion (10bp)
            - strand)

            --------------------------- bkp3'
                                                                bkp5' --------------------
                                           <------------------>
                                              deletion (10bp)

            Input:
            1) informative5primeContigObj.
            2) informative3primeContigObj.

            Output:
            1) targetSiteSize. Target site duplication length. 'inconsistent' if expected TSD size different to TSD sequence length.
            2) targetSiteSeq. Target site duplication sequence or 'na' if no TSD. 'inconsistent' if expected TSD size different to TSD sequence length.
        """

        bkpPos5prime = informative5primeContigObj.informativeDict["bkp"][1]
        bkpPos3prime = informative3primeContigObj.informativeDict["bkp"][1]
        alignObj5prime = informative5primeContigObj.informativeDict["targetRegionAlignObj"]

        ## Compute TSD length
        targetSiteSize = abs(bkpPos5prime - bkpPos3prime)

        ## Extract TSD sequence
        # A) Begin of the contig sequence aligned in the TE insertion genomic region
        #   -------------**TSD**######TE#####
        #   --------------------
        # qBeg               *qEnd*
        if (alignObj5prime.alignType == "beg"):
            beg = alignObj5prime.qEnd - targetSiteSize
            end = alignObj5prime.qEnd
            targetSiteSeq = informative5primeContigObj.seq[beg:end]

        # B) End of the contig sequence aligned in the TE insertion genomic region
        #   ######TE#####AAAAAAA**TSD**-------------
        #                       --------------------
        #                    *qBeg*               qEnd
        else:
            beg = alignObj5prime.qBeg    # (no substract 1 since psl coordinates are 0-based as python strings)
            end = alignObj5prime.qBeg + targetSiteSize
            targetSiteSeq = informative5primeContigObj.seq[beg:end]

        ## Inconsistent TSD if sequence has not the expected length or longer than 100 bp
        if (targetSiteSize != len(targetSiteSeq)) or (targetSiteSize > 100):
            targetSiteSize = "inconsistent"
            targetSiteSeq = "inconsistent"

        return (targetSiteSize, targetSiteSeq)


    def insertion_orientation_polyA(self, informative5primeContigObj, informative3primeContigObj):
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
        if (informative5primeContigObj != "UNK") and (informative3primeContigObj != "UNK"):
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
        elif (informative5primeContigObj != "UNK"):
            alignType5prime = informative5primeContigObj.informativeDict["targetRegionAlignObj"].alignType

            # a) + strand
            if (alignType5prime == "beg"):
                orientation = "+"

            # b) - strand
            else:
                orientation = "-"

        ## C) 3' informative contig
        elif (informative3primeContigObj != "UNK"):
            alignType3prime = informative3primeContigObj.informativeDict["targetRegionAlignObj"].alignType

            # a) + strand
            if (alignType3prime == "end"):
                orientation = "+"

            # b) - strand
            else:
                orientation = "-"

        ## D) None informative contig
        else:
            orientation = "UNK"

        return orientation

    def insertion_orientation_ERVK(self, informative5primeContigObj, informative3primeContigObj):
        """
            Determine ERVK insertion strand/orientation.

            If informative contigs have a consistent mapping direction
            over their length, the element has been in normal sense:

            --------|#####TE#####
            >>>>>>>>|>>>>>>>>>>>>

            A switch in mapping direction indicates that the element has
            been inserted in an inverted direction:

            --------|#####TE#####
            >>>>>>>>|<<<<<<<<<<<<

            Input:
            1) informative5primeContigObj
            2) informative3primeContigObj

            Output:
            1) orientation. +, -, UNK or inconsistent dna strand.

            Note: Inconsistent when 5' and 3' informative contigs suggest
                  contradictory orientations (should not happen).
        """

        # check if required alignment members are present
        has_req_aln_5prime = (
            (informative5primeContigObj == "UNK") or
            (
             (informative5primeContigObj.informativeDict["targetRegionAlignObj"] and
              isinstance(informative5primeContigObj.informativeDict["targetRegionAlignObj"], blat_alignment))
             and
             (informative5primeContigObj.informativeDict["info"] and
              isinstance(informative5primeContigObj.informativeDict["info"], blat_alignment))
            )
        )
        has_req_aln_3prime = (
            (informative3primeContigObj == "UNK") or
            (
             (informative3primeContigObj.informativeDict["targetRegionAlignObj"] and
              isinstance(informative3primeContigObj.informativeDict["targetRegionAlignObj"], blat_alignment))
             and
             (informative3primeContigObj.informativeDict["info"] and
              isinstance(informative3primeContigObj.informativeDict["info"], blat_alignment))
            )
        )
        has_required_alignments = has_req_aln_5prime and has_req_aln_3prime

        # if any of the expected alignments is missing:
        # log error and set orientation to "inconsistent"
        if (not has_required_alignments):
            log("ERROR", "Alignment object missing for insertion '%s'" % self.traficId)
            orientation = "inconsistent"
            return orientation

        ## A) 5' and 3' informative contigs
        if (informative5primeContigObj != "UNK") and (informative3primeContigObj != "UNK"):

            # get alignment directions
            targDir5prime = informative5primeContigObj.informativeDict["targetRegionAlignObj"].strand
            consDir5prime = informative5primeContigObj.informativeDict["info"].strand
            targDir3prime = informative3primeContigObj.informativeDict["targetRegionAlignObj"].strand
            consDir3prime = informative3primeContigObj.informativeDict["info"].strand

            # a) + strand
            if (targDir5prime == consDir5prime) and (targDir3prime == consDir3prime):
                orientation = "+"

            # b) - strand
            elif (targDir5prime != consDir5prime) and (targDir3prime != consDir3prime):
                orientation = "-"

            # c) Contradictory/inconsistent
            else:
                #orientation = "inconsistent"
                orientation = "UNK"

        ## B) 5' informative contig
        elif (informative5primeContigObj != "UNK"):

            # get alignment directions
            targDir5prime = informative5primeContigObj.informativeDict["targetRegionAlignObj"].strand
            consDir5prime = informative5primeContigObj.informativeDict["info"].strand

            # a) + strand
            if (targDir5Prime == consDir5prime):
                orientation = "+"

            # b) - strand
            else:
                orientation = "-"

        ## C) 3' informative contig
        elif (informative3primeContigObj != "UNK"):

            # get alignment directions
            targDir3prime = informative3primeContigObj.informativeDict["targetRegionAlignObj"].strand
            consDir3prime = informative3primeContigObj.informativeDict["info"].strand

            # a) + strand
            if (targDir3prime == consDir3prime):
                orientation = "+"

            # b) - strand
            else:
                orientation = "-"

        ## D) None informative contig
        else:
            orientation = "UNK"

        return orientation


    def insertion_structure(self, informative5primeContigObj):
        """
            Determine TE (L1, Alu, SVA or ERVK) insertion structure.

            1) 5'inverted:

                TE in + orientation with 5'inversion    ----#######TE######AAAAA----
                                                            <<<<<<>>>>>>>>>>
                                                            <---->
                                                           inversion
                5-prime informative contig              --------- (5'inversion signature: the piece of contig corresponding to the TE
                                                                   aligns in the opposite DNA strand than the piece of contig aligning in the target region)

            2) Full length insetion:

                                                     -------#######TE######AAAAA----
            5-prime informative contig                  ---------

            3) 5' truncated insertion:

                                                       --------###TE######AAAAA----
                                                           <-->
                                                         deletion
                5-prime informative contig              ---____--- (5'truncation signature: the piece of contig corresponding to TE
                                                                    aligns in the body of the TE and not in the 5' extreme)

            Input:
            1) informative5primeContigObj

            Output:
            1) structure. One of 'na', 'INV', 'DEL' or 'FULL'
            2) length. Inserted TE length, 'na' if not available.
            3) percLength. Percentage of TE consensus sequence inserted, 'na' if not available.
        """

        ## A) Not TD2, not ERVK and 5' informative contig
        if (self.tdType != "TD2") and (self.TEClass != "ERVK") and (informative5primeContigObj != "UNK"):

            strand = informative5primeContigObj.informativeDict["info"].strand

            ## Determine TE insertion structure

            # a) TE inverted in its 5'
            if (strand == "-"):
                structure = "INV"
                length = "UNK"
                percLength = "UNK"

            # b) Full length or 5' truncated
            else:
                tBeg = informative5primeContigObj.informativeDict["info"].tBeg
                tSize = informative5primeContigObj.informativeDict["info"].tSize

                length = tSize - tBeg
                percLength = float(length) / tSize * 100

                # b.a) full length TE insertion
                if (percLength > 95):
                    # L1 (6021 bp length + 30bp polyA, first ~300bp correspond to promoter)
                    # Alu (282 bp length + 30bp polyA).
                    # SVA (X bp length + 30bp polyA).
                    # ERVK (X bp length + 30bp polyA).
                    structure = "FULL"

                # b.b) 5' truncated
                else:
                    structure = "DEL"

        ## B) TD2 or ERVK or No 5' informative contig 
        else:
            structure = "UNK"
            length = "UNK"
            percLength = "UNK"

        return (structure, length, percLength)


    def imprecise_bkp(self, insertionCoord):
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

        insertionCoordList = insertionCoord.split("_")

        chrom = str(insertionCoordList[0])
        beg = int(insertionCoordList[1])
        end = int(insertionCoordList[2])

        meanPos = (beg + end)/2
        dist = abs(end - beg)
        CI = dist/2

        breakpoint = [chrom, meanPos, CI]

        return breakpoint

class cluster():
    """
    Transposable element insertion cluster class.

    A cluster can be + or - and it has associated the contigs resulting from the assembly of TE insertion supporting
    reads identified by TraFiC.

    Methods:
    - blat_alignment_reader
    - create_contigs_dict
    - add_alignments
    - find_informative_contig
    """

    def __init__(self, ID, contigsFasta, readPairsList):
        """
            Initialize cluster object.

            Input:
            1) ID. Cluster id.
            2) contigsFasta. Fasta file containing the assembled contigs for the given cluster
            3) readPairsList. List of cluster supporting reads

            Output:
            - Cluster object variables initialized
        """
        self.ID = ID
        self.contigsDict = self.create_contigs_dict(contigsFasta)
        self.readPairIds = readPairsList
        self.nbPairs =  len(readPairsList.split(','))

    #### FUNCTIONS ####
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
        fastaObj = fasta(contigsFasta)

        contigsDict = {}

        ### For each contig create a contig object and add it to the dictionary
        # using the contig id as key
        for contigId in fastaObj.fastaDict:

            contigSeq = fastaObj.fastaDict[contigId]

            # Create contig object
            contigObj = contig(contigId, contigSeq)

            # Add contig object to the dictionary
            contigsDict[contigId] = contigObj

        return contigsDict

    def add_alignments(self, blatPath):
        """
            Read a psl file containing the contig blat alignments on the reference genome and TE sequences and associate
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

    def find_informative_contigs(self, insertionCoord, tdType, srcCoord):
        """
            Identify 5' and 3' informative contig belonging to the cluster. Informative 5' and 3' contigs span 5' and 3' insertion breakpoints, respectively.

            Input:
            1) insertionCoord. Region of interest. Format: ${chrom}_${beg}_${end}.
                               Example: 10_108820680_108820678.
            2) tdType. Type of transduction event ("TD0": solo, "TD1": partnered, "TD2": orphan)
            3) srcCoord. Mobilized region (relevant for TD2 and PSD events)

            Output:
            1) bestInformative5primeContigObj. Best 5' Informative contig object. 'na' if not found
            2) bestInformative3primeContigObj. Best 3' Informative contig object. 'na' if not found
        """

        ## Initial status -> none informative contig
        bestInformative5primeContigObj = "UNK"
        bestInformative3primeContigObj = "UNK"

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
            informative = contigObj.is_informative(insertionCoord, tdType, srcCoord)

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

            if (bestInformative5primeContigObj == "UNK") or (dist < bestDist):

                bestInformative5primeContigObj = contigObj
                bestDist = dist

        # 2) informative 3-prime

        bestDist = ""

        for contigObj in  informative3primeContigObjList:

            bkpCoord = contigObj.informativeDict["bkp"][1]
            dist = abs(targetPos - bkpCoord)

            if (bestInformative3primeContigObj == "UNK") or (dist < bestDist):

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

    def __init__(self, contigId, contigSeq):
        """
            Initialize contig object.

            Input:
            1) contigTuple. First element (contig Id) and second element (contig Sequence)

            Output:
            - Contig object variables initialized
        """
        # Contig id and sequence
        self.ID = contigId
        self.seq = contigSeq
        self.length = int(len(self.seq))

        # Cluster the contig belongs to
        self.cluster = ""

        # Contig alignment information
        self.alignList = []   # List of blat alignments for this contig
        self.informativeDict = {} # Dictionary key -> value pairs:
        self.informativeDict["type"] = "UNK"                 # type -> 5-prime, 3-prime or none
        self.informativeDict["bkp"] = "UNK"                  # bkp -> list: chrom and pos, 'na' for both if type == none)
        self.informativeDict["info"] = "UNK"                 # info -> 5-prime: aligment object with contig's alignment in TE sequence info; 3-prime: PolyA sequence; none: 'na'
        self.informativeDict["targetRegionAlignObj"] = "UNK" # targetRegionAlignObj -> alignment object with contig's alignment in the target region info.
        self.informativeDict["MEISeq"] = "UNK"               # MEISeq -> Assembled sequence of the mobile element 5' boundary

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
        Contigs completely aligning to either the insertion site or the source element (TE or transduced seq) are not considered as candidates.

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

            # 1. Check if contig completely aligns on the consensus TE sequence
            inConsensusTE = alignment.in_consensus_TE()

            # A) Contig not completely aligning on the consensus TE
            if (inConsensusTE == 0):

                # 2. Check if alignment within the target region
                insertionRegion = alignment.in_target_region(insertionCoord, windowSize)

                # Within target region
                if (insertionRegion == 1):

                    # 3. Check if it is a partial alignment
                    partial = alignment.partial_alignment(maxAlignPerc)

                    # Partial
                    if (partial == 1):
                        supportingAlignList.append(alignment)

                        # Informative candidate contig -> partially alignining on the target region and do not aligning completely on the TE sequence
                        candidate = 1

            # B) Not informative contig as it completely aligns on the consensus TE
            else:
                candidate = 0
                supportingAlignList = []
                break

        return (candidate, supportingAlignList)


    def is_informative(self, insertionCoord, tdType, srcCoord):
        """
        Check if candidate contig is 5' or 3' informative. Defined as contigs spanning
        5' or 3' insertion breakpoints:

        5' informative         ####TE####---------------
        3' informative         ---------------AAAAAAAAAA


        Input:
            1) insertionCoord. Region of interest. Format: ${chrom}_${beg}_${end}.
                               Example: 10_108820680_108820678.
            2) tdType. Type of transduction event
                       ("TD0": solo,
                        "TD1": partnered,
                        "TD2": orphan,
                        "PSD": pseudogene insertion)
            3) srcCoord. Mobilized region (relevant for TD2 and PSD events)

        Output:
            1) informativeBoolean. Boolean, 1 (informative) and 0 (not informative)
            2) Sets 'informativeDict' variable. Dictionary key -> value pairs:
                type -> 5-prime, 3-prime or none
                bkp -> list: chrom and pos, 'na' for both if type == none)
                info -> 5-prime: aligment object with contig's alignment in TE sequence info;
                        3-prime: PolyA sequence; none: 'na'
                targetRegionAlignObj -> alignment object with contig's alignment in the target region info.
                MEISeq -> Assembled sequence of the mobile element 5' boundary
        """

        ## Initial status -> no informative
        bestInformative5primeDict = "UNK"
        bestInformative3primeDict = "UNK"
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
                is5prime, bkpCoord, srcAlignmentObj, MEISeq = self.is_5prime_informative(alignObj, tdType, srcCoord)

                # Contig informative 5-prime, it has TE sequence
                if (is5prime == 1):

                    informative5primeDict[alignObj] = {}
                    informative5primeDict[alignObj]["type"] = "5-prime"
                    informative5primeDict[alignObj]["bkp"] = bkpCoord
                    informative5primeDict[alignObj]["info"] = srcAlignmentObj
                    informative5primeDict[alignObj]["targetRegionAlignObj"] = alignObj
                    informative5primeDict[alignObj]["MEISeq"] = MEISeq

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

        for alignment in informative5primeDict:
            bkpCoord = informative5primeDict[alignment]["bkp"][1]
            dist = abs(targetPos - bkpCoord)

            if (bestInformative5primeDict == "UNK") or (dist < bestDist5Prime):
                bestInformative5primeDict = informative5primeDict[alignment]
                bestDist5Prime = dist

        # 2) informative 3-prime
        bestDist3Prime = ""

        for alignment in informative3primeDict:
            bkpCoord = informative3primeDict[alignment]["bkp"][1]
            dist = abs(targetPos - bkpCoord)

            if (bestInformative3primeDict == "UNK") or (dist < bestDist3Prime):
                bestInformative3primeDict = informative3primeDict[alignment]
                bestDist3Prime = dist

        ## In case there are alignments supporting the contig is 5' and 3' select the best
        # Note: Best defined as alignment supporting putative insertion breakpoint closer to the
        # insertion target region

        # A) Any alignment supporting the contig as 5' or 3' informative
        if (bestDist5Prime == "") and (bestDist3Prime == ""):
            informativeBoolean = 0

        # B) Best 5'
        elif (bestDist5Prime < bestDist3Prime) or ((bestDist5Prime != "") and (bestDist3Prime == "")):
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
            bkpCoord = ["UNK", "UNK"]
            polyASeq = "UNK"

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

    def is_5prime_informative(self, alignObj, tdType, srcCoord):
        """
        Check if candidate contig is 5' informative. Defined as contigs spanning
        source element (TE or transduced seq.) - insertion target region breakpoint:

        5' informative:         ---------------###SRC####
        5' informative:         ###SRC####---------------

        Input:
        1) alignObj. Blat alignment object.
        2) tdType. Type of transduction event
                   ("TD0": solo,
                    "TD1": partnered,
                    "TD2": orphan,
                    "PSD": pseudogene insertion)
        3) srcCoord. mobilized region (relevant for TD2 and PSD events)

        Output:
        1) is5prime. Boolean, 1 (5' informative) and 0 (not 5' informative).
        2) bkpCoord. Two elements breakpoint coordinates list. First (bkp chromosome) and second (breakpoint position).
        3) srcAlignmentObj. Blat aligment object with the alignment information of the contig in the source sequence (TE or TD).
                           'na' if not 5' informative.
        4) MEISeq. Assembled sequence of the mobile element 5' boundary
        """

        ## Select contig target sequence coordinates to search for alignment in source sequence:
        #   TD0, TD1 events: consensus L1, Alu or SVA for RD
        #   TD2 events: transduced region downstream of TE
        # The position of the target coordinates in the contig
        # will depend on the blat alignment type

        # A) Begin of the contig sequence aligned in the TE insertion genomic region
        #   -------------******SRC******
        #   -------------
        # qBeg        *qEnd*

        if (alignObj.alignType == "beg"):
            targetBeg = alignObj.qEnd
            targetEnd = alignObj.qSize
            targetSeq = self.seq[alignObj.qEnd:]
            bkpChrom = alignObj.tName

            # a) Positive dna strand
            if (alignObj.strand == "+"):
                bkpPos = alignObj.tEnd
                MEISeq = targetSeq

            # b) Negative dna strand
            else:
                bkpPos = alignObj.tBeg
                MEISeq = self.rev_complement(targetSeq)

        # B) End of the contig sequence aligned in the TE insertion genomic region
        #   *****SRC******-------------
        #                 -------------
        #              *qBeg*        qEnd
        elif (alignObj.alignType == "end"):
            targetBeg = 0
            targetEnd = alignObj.qBeg
            targetSeq = self.seq[:alignObj.qBeg]

            bkpChrom = alignObj.tName

            # a) Positive dna strand
            if (alignObj.strand == "+"):
                bkpPos = alignObj.tBeg
                MEISeq = targetSeq

            # b) Negative dna strand
            else:
                bkpPos = alignObj.tEnd
                MEISeq = self.rev_complement(targetSeq)

        # C) No align type information or 'none' align type
        else:
            log("Error", "No valid alignment object provided. Alignment type variable is 'none' or not defined")
            sys.exit(1)


        ## Default
        is5prime = 0
        bkpCoord = ["UNK", "UNK"]
        srcAlignmentObj = "UNK"
        targetMEIs = ["L1", "Alu", "SVA", "ERVK"]

        for alignment in self.alignList:

            # a breakpoint-spanning contig is expected to align partly to
            #   TD0 (solo insertion) : consensus MEI sequence
            #   TD1 (partnered TD)   : consensus MEI sequence
            #   TD2 (orphan TD)      : transduced region downstream of TE
            has_expected_target = (
              ( (tdType not in ["TD2", "PSD"]) and (alignment.tName in targetMEIs) ) or
              ( (tdType in ["TD2", "PSD"]) and (alignment.in_target_region(srcCoord, 1000)) )
            )

            # Contig alignment in TE sequence (L1, Alu, SVA or ERVK) or transduced region
            if has_expected_target:

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
                    srcAlignmentObj = alignment
                    # NOTE: break from loop here? could following alignments be better?


        return (is5prime, bkpCoord, srcAlignmentObj, MEISeq)

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

    def in_consensus_TE(self):
        """
            Check if blat alignment is completely on the mobile element consensus sequence. No gap is allowed
        """
        alignLength = self.qEnd - self.qBeg
        alignPerc = float(alignLength) / float(self.qSize) * 100      

        MEIlist = ["L1", "Alu", "SVA", "ERVK"]

        # A) More than 99% Contig aligning on the mobile element consensus sequence. No gap allowed       
        if (self.tName in MEIlist) and (alignPerc > 99) and (self.blockCount == 1):
            inConsensusTE = 1

        # B) Contig do not aligning on the consensus seq.
        else:
            inConsensusTE = 0            

        return inConsensusTE
        
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

        # A) One block partial alignment 
        #if (alignPerc < maxAlignPerc) and (self.blockCount == 1):
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



class fasta():
    """
    .... class.

    .....

    Methods:
    - fasta_reader

    """

    def __init__(self, fastaFile):
        """
            Initialize fasta object.

            Input:
            1)

            Output:
            -
        """
        self.fastaDict = self.fasta_reader(fastaFile)

    #### FUNCTIONS ####
    def fasta_reader(self, fastaFile):
        """


            Input:
            1)

            Output:
            1)
        """
        fastaDict = {}

        subHeader("Fasta reader")

        fh = open(fastaFile)
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()

            # drop the info
            header = header.split(" ")[0]

            info("Reading " + header + "...")
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            fastaDict[header] = seq

        return fastaDict

def parse_args():
    """Define and parse command line parameters."""
    parser = argparse.ArgumentParser(description="Per MEI called by TraFiC: 1) Identifies informative contigs spanning 5' and/or 3' insertion ends if possible, 2) Use informative contigs for characterizing MEI in detail (exact breakpoints, length, strand...) and 3) Produce a VCF with the MAI plus all these information")
    parser.add_argument('inputPaths', help='Text file containing, per MEI, the needed files')
    parser.add_argument('sampleId', help='Sample identifier. The output vcf will be named accordingly')
    parser.add_argument('genome', help='Reference genome in fasta format')
    parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

    args = parser.parse_args()
    return args

#### MAIN ####

if __name__ == "__main__":
    ## Get user's input ##
    args = parse_args()
    inputPaths = args.inputPaths
    sampleId = args.sampleId
    genome = args.genome
    outDir = args.outDir

    scriptName = os.path.basename(sys.argv[0])

    ## Display configuration to standard output ##
    print
    print "***** ", scriptName, " configuration *****"
    print "paths2bkpAnalysis: ", inputPaths
    print "sampleId: ", sampleId
    print "genome: ", genome
    print "outDir: ", outDir
    print
    print "***** Executing ", scriptName, " *****"
    print


    ## Start ## 

    outFilePath = outDir + '/' + sampleId + '.vcf'

    ## 0. Create reference genome fasta object

    header("Creating reference genome fasta object")
    genomeObj = fasta(genome)

    ## 1. Create VCF object and print VCF header

    header("Creating VCF object and printing VCF header into the output file")
    VCFObj = VCF()
    VCFObj.print_header(outFilePath)

    ## 2. Per each insertion perform breakpoint analysis

    inputFile = open(inputPaths, 'r')

    # Analyze one insertion per iteration
    for line in inputFile:
        line = line.rstrip('\n')
        line = line.split("\t")

        # Get TE insertion info and files
        insertionInfo = line[0]
        TEClass, tdType, insertionCoord = insertionInfo.split(":")
        contigsPlusPath, contigsMinusPath = line[1].split(",")
        blatPlusPath, blatMinusPath = line[2].split(",")
        readPairsPlus = line[3]
        readPairsMinus = line[4]
        srcElement = line[5]
        transductionInfo = line[6]
        pseudogeneInfo = line[7]

        # Perform breakpoint analysis for the TE insertion
        header("Tranposable Element Insertion Breakpoint Analysis (TEIBA) for: " + insertionCoord)

        # A) All the input files exist
        if os.path.isfile(contigsPlusPath) and os.path.isfile(blatPlusPath) and os.path.isfile(contigsMinusPath) and os.path.isfile(blatMinusPath):

            ## Create insertion object and identify breakpoints from assembled contigs
            insertionObj = insertion(TEClass, tdType, insertionCoord, contigsPlusPath, blatPlusPath, contigsMinusPath, blatMinusPath, readPairsPlus, readPairsMinus, srcElement, transductionInfo, pseudogeneInfo, sampleId)
            insertionObj.find_insertionBkp(genomeObj, outDir)

            ## Create VCFline object
            VCFlineObj = VCFline(insertionObj, genomeObj)

            ## Add VCFline to the list in VCF object
            VCFObj.addLine(VCFlineObj)

        else:
            message = "Input files for " + insertionCoord + " insertion do not exist"
            log("ERROR", message)

    ## 3. Write lines describing the TE insertions into the VCF file
    VCFObj.print_lines(outFilePath)


    ## Finish ##
    print
    print "***** Finished! *****"
    print
