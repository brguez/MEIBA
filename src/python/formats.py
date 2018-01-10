#!/usr/bin/env python
#coding: utf-8

import re
import time
import itertools 

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


class VCF():
    """
        .....................

        Methods:
        -

    """

    def __init__(self):
        """

        """
        self.header = ""
        self.lineList = []  # List of VCFline objects

    #### METHODS ####
    def addLine(self, VCFlineObj):
        """

        """
        self.lineList.append(VCFlineObj)


    def sort(self):
        """

        """

        lineListSorted = sorted(self.lineList, key=lambda line: (line.chrom, line.pos))

        return lineListSorted


    def read_VCF(self, VCFpath):
        """

        """
        VCFfile = open(VCFpath, 'r')
        headerList = []

        for line in VCFfile:

            # A) Header
            if line.startswith("##"):
                headerList.append(line)

            # B) Variants
            elif not line.startswith("#"):
                line = line.rstrip('\n')
                line = line.split('\t')
                VCFlineObj = VCFline(line)
                self.addLine(VCFlineObj)

        self.header = "".join(headerList)

    def read_VCF_multiSample(self, VCFpath):
        """
        """

        VCFfile = open(VCFpath, 'r')
        headerList = []

        for line in VCFfile:

            # A) Meta-information
            if line.startswith("##"):
                headerList.append(line)

            # B) Header
            elif line.startswith("#CHROM"):
                line = line.rstrip('\n')
                line = line.split('\t')
                donorIdList = [line[i] for i in range(9, len(line))]

            # C) Variants
            elif not line.startswith("#"):
                line = line.rstrip('\n')
                line = line.split('\t')
                VCFlineObj = VCFline(line)
                gnTypeList = [line[i] for i in range(9, len(line))]

                for i in range(0, len(gnTypeList)):

                    donorId = donorIdList[i]
                    gnType = gnTypeList[i]
                    VCFlineObj.genotypesDict[donorId] = gnType

                self.addLine(VCFlineObj)

        self.header = "".join(headerList)

        return(donorIdList)

    def write_header(self, outFilePath):
        """

        """

        outFile = open(outFilePath, 'w')
        outFile.write(self.header)
        outFile.close()

    def write_variants(self, outFilePath):
        """

        """

        outFile = open(outFilePath, 'a')

        ## 1. Write VCF header
        header = "#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT"  + "\n"

        outFile.write(header)

        ## 2. Write variants
        # Iterate and write each VCF data line into the output VCF file
        for VCFline in self.lineList:

            row = VCFline.chrom + "\t" + str(VCFline.pos) + "\t" + VCFline.id + "\t" + VCFline.ref + "\t" + VCFline.alt + "\t" + VCFline.qual + "\t" + VCFline.filter + "\t" + VCFline.info + "\t" + VCFline.format + "\t" + VCFline.genotype + "\n"
            outFile.write(row)

        ## Close output file
        outFile.close()


    def write_variants_multiSample(self, donorIdList, outFilePath):
        """
        """

        outFile = open(outFilePath, 'a')

        ## 1. Write VCF header
        # Make fixed fields
        header = "#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT"

        # Add donor ids to the header
        for donorId in donorIdList:
            header =  header + "\t" + donorId

        # Add trailing new line and write header
        header = header + "\n"
        outFile.write(header)

        ## 2. Write genotyped variants
        # Iterate and write each VCF data line into the output VCF file
        for VCFline in self.lineList:

            row = VCFline.chrom + "\t" + str(VCFline.pos) + "\t" + VCFline.id + "\t" + VCFline.ref + "\t" + VCFline.alt + "\t" + VCFline.qual + "\t" + VCFline.filter + "\t" + VCFline.info + "\t" + "GT:NV:NR"

            # Add donor's genotypes to the VCF data line. One per iteration
            for donorId in donorIdList:
                genotype = VCFline.genotypesDict[donorId]
                row =  row + "\t" + genotype

            # Add trailing new line and write VCF data line
            row = row + "\n"
            outFile.write(row)


    def create_header(self):
        """
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
##INFO=<ID=CLASS,Number=1,Type=String,Description="Mobile element class (L1, ALU, SVA or ERVK)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Insertion type (TD0: solo, TD1: partnered-3'transduction, TD2: orphan-3'transduction, PSD: processed-pseudogene)>
##INFO=<ID=MECHANISM,Number=1,Type=String,Description="Insertion mechanism (TPRT: target primed reverse transcription, EI: endonuclease independent, DPA: double poly-A)>
##INFO=<ID=SCORE,Number=1,Type=Integer,Description="Insertion score (5: 5' and 3' breakpoints (bkp) identified, 4: 3'bkp identified, 3: 5'bkp identified, 2: no bkp identified, 1: inconsistent (contradictory insertion features))">
##INFO=<ID=MANUAL,Number=0,Type=Flag,Description="MEI manually verified and curated through BAM inspection (Only used for PSD)">
##INFO=<ID=BKPB,Number=1,Type=Integer,Description="MEI right-most breakpoint position (bkp B). Left-most breakpoint position (bkp A) represented in the POS field">
##INFO=<ID=CIPOS,Number=1,Type=Integer,Description="Confidence interval around insertion breakpoint">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion DNA strand (+ or -)">
##INFO=<ID=STRUCT,Number=1,Type=String,Description="Mobile element structure (INV: 5'inverted, DEL: 5'deleted, FULL: full-length)">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Mobile element length">
##INFO=<ID=TSLEN,Number=1,Type=Integer,Description="Target site duplication (+_value) or deletion (-_value) length">
##INFO=<ID=SRCID,Number=1,Type=String,Description="Source element cytoband identifier. Only for gemline source elements">
##INFO=<ID=SRCTYPE,Number=1,Type=String,Description="Source element type (GERMLINE or SOMATIC)">
##INFO=<ID=SRC,Number=1,Type=String,Description="Coordinates of the source element that mediated the transduction in the format: chrom_beg_end">
##INFO=<ID=TDC,Number=1,Type=String,Description="Begin and end coordinates of the integrated transduced or pseudogene sequence in the format: chrom_beg_end">
##INFO=<ID=TDLEN,Number=1,Type=Integer,Description="Transduced region length">
##INFO=<ID=TDLENR,Number=1,Type=Integer,Description="Transduced region length at RNA level">
##INFO=<ID=SRCGENE,Number=1,Type=String,Description="Source gene of the processed pseudogene insertion">
##INFO=<ID=GR,Number=1,Type=String,Description="L1-mediated genomic rearrangement (DEL: deletion, DUP: duplication or TRANS: translocation)">
##INFO=<ID=GERMDB,Number=1,Type=String,Description="MEI reported as germinal in a database">
##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region where the mobile element is inserted (exonic, splicing, ncRNA, UTR5, UTR3, intronic, upstream, downstream, intergenic)">
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO gene symbol">
##INFO=<ID=ROLE,Number=1,Type=String,Description="Role in cancer (oncogene, TSG: tumor suppressor gene, oncogene/TSG: both roles)">
##INFO=<ID=COSMIC,Number=0,Type=Flag,Description="Reported as cancer driver in COSMIC cancer gene census database">
##INFO=<ID=CPG,Number=0,Type=Flag,Description="Reported as cancer predisposition gene in 10.1038/nature12981 (DOI).">
##INFO=<ID=REP,Number=1,Type=String,Description="Repetitive element overlapping the insertion breakpoints">
##INFO=<ID=DIV,Number=1,Type=Integer,Description="Millidivergence of the overlapping repetitive element with respect a consensus sequence">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Assembled contig sequence spanning 1st bkp (lowest genomic position)">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Assembled contig sequence spanning 2nd bkp (highest genomic position)">
##INFO=<ID=RP,Number=.,Type=String,Description="Reads from the tumour sample and positive cluster that support this insertion">
##INFO=<ID=RN,Number=.,Type=String,Description="Reads from the tumour sample and negative cluster that support this insertion">
##FILTER=<ID=SCORE,Description="Insertion with an score < threshold">
##FILTER=<ID=FPSOURCE,Description="Transduction mediated by a false somatic source element">
##FILTER=<ID=REP,Description="Insertion overlapping a satellite region or a repetitive element of the same class">
##FILTER=<ID=DUP,Description="Duplicated MEI call">
##FILTER=<ID=GERMLINE,Description="Germline MEI miscalled as somatic">
##FILTER=<ID=TD,Description="L1 transduction incorrectly identified as a processed pseudogene insertion">
##FILTER=<ID=COUNT,Description="Allele count < threshold">
##FILTER=<ID=MISSGT,Description="Genotype missing-ness rate > threshold ">
##FORMAT=<ID=RCP,Number=1,Type=Integer,Description="Count of positive cluster supporting reads">
##FORMAT=<ID=RCN,Number=1,Type=Integer,Description="Count of negative cluster supporting reads">
##FORMAT=<ID=SL,Number=1,Type=String,Description="List of samples where the variant was found (specially relevant for multi-tumor donors)">
##FORMAT=<ID=REPR,Number=1,Type=String,Description="Sample selected as representative among all the samples where the variant was found (specially relevant for multi-sample VCF).">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting the variant in this sample">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads covering variant location in this sample">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  
"""

        self.header = template.format(**context)


class VCFline():
    """

    """

    def __init__(self, VCFlineList):
        """

        """
        self.chrom = VCFlineList[0]
        self.pos = int(VCFlineList[1])
        self.id = VCFlineList[2]
        self.ref = VCFlineList[3]
        self.alt = VCFlineList[4]
        self.qual = VCFlineList[5]
        self.filter = VCFlineList[6]
        self.info = VCFlineList[7]
        self.format = VCFlineList[8]
        self.genotype = VCFlineList[9]
        self.genotypesDict = {}
        self.infoDict = self.read_info()

    def read_info(self):
        """

        """
        infoDict = {}

        infoList = self.info.split(';')

        # Iterate over info fields list.
        for field in infoList:

            fieldList = field.split('=')

            # A) Key-value pair
            if (len(fieldList) == 2):

                key = fieldList[0]
                value = fieldList[1]

            # B) Flag
            else:
                key = fieldList[0]
                value = "1"

            infoDict[key] = value

        return infoDict

    def make_info(self):
        """
        """

        ## Create list containing the order of info fields 
        infoOrder = [ "SVTYPE", "CLASS", "TYPE", "MECHANISM", "SCORE", "MANUAL", "BKPB", "CIPOS", "STRAND", "STRUCT", "LEN", "TSLEN", "SRCID", "SRCTYPE", "SRC", "TDC", "TDLEN", "TDLENR", "SRCGENE", "GR", "GERMDB", "POLYMORPHIC", "NOVEL", "REGION", "GENE", "ROLE", "COSMIC", "CPG", "REP", "DIV", "CONTIGA", "CONTIGB", "RP", "RN" ]

        flagList = ["POLYMORPHIC", "NOVEL", "COSMIC", "CPG" , "MANUAL"]

        ## Create info string in the correct order from dictionary
        infoList = []

        # Iterate over all possible info field keys
        for info in infoOrder:

            # Information about current info field available and applicable
            if (info in self.infoDict.keys()) and (self.infoDict[info] != "UNK") and (self.infoDict[info] != "NA"):

                # A) Flag
                if (info in flagList):

                    # Flag with positive value
                    if (self.infoDict[info] == "1"):
                        
                        infoField = info
                        infoList.append(infoField)

                # B) Key-value pair
                else:
                    infoField = info + "=" + str(self.infoDict[info])
                    infoList.append(infoField)

        # Concatenate elements in a single string separated by commas
        info = ';'.join(infoList)

        return(info)


    def nbPE(self):
        """
        Return the total number of paired-end reads supporting a MEI (sum + and - cluster supporting reads)
        """
        genotypeList = self.genotype.split(':')
        nbPE = int(genotypeList[0]) + int(genotypeList[1])

        return nbPE

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
        faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == ">"))
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


class annovar():
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
        self.lineList = []  # List of annovarLine objects

    #### METHODS ####
    def read_annovar(self, annovarPath):
        """

        """
        annovarFile = open(annovarPath, 'r')

        for line in annovarFile:

            line = line.rstrip('\n')
            annovarLineList = line.split('\t')
            annovarLineObj = annovarLine(annovarLineList)
            self.addLine(annovarLineObj)

    def addLine(self, annovarLineObj):
        """

        """
        self.lineList.append(annovarLineObj)

    def info2VCF(self, VCFObj):
        """

        """

        annotDictList = []

        for i in range(0, len(self.lineList), 1):

            annovarLineObj = self.lineList[i]
            VCFlineObj =  VCFObj.lineList[i]

            if (annovarLineObj.chrom == VCFlineObj.chrom) and (annovarLineObj.beg == VCFlineObj.pos):

                annotDict = {}
                annotDict["REGION"] = annovarLineObj.region

                if (annotDict["REGION"] != "intergenic"):
                    annotDict["GENE"] = annovarLineObj.genes

                VCFlineObj.infoDict.update(annotDict)
                VCFlineObj.info = VCFlineObj.make_info()

            else:
                print "ERROR. Unpaired VCF and annovar variants"


class annovarLine():
    """
        .....................

        Methods:
        -

    """

    def __init__(self, annovarLineList):
        """
            ...

            Output:
            -
        """
        self.region = annovarLineList[0]
        self.genes = re.sub('\(.*?\)', '', annovarLineList[1]) # remove not necesary information enclosed by slashes (.*)
        self.chrom = annovarLineList[2]
        self.beg = int(annovarLineList[3])
        self.end = int(annovarLineList[4])
        self.ref = annovarLineList[5]
        self.alt = annovarLineList[6]
