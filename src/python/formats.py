#!/usr/bin/env python
#coding: utf-8

import re
import time
import itertools 
import os 

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


    def create_header(self, genome, genomeObj):
        """
        """

        ### 1. Define variables
        ## 1.1 date
        date = time.strftime("%Y%m%d")

        ## 1.2 Extract species, build and version from genome file name
        species, build, version, nan = os.path.basename(genome).split("_")

        if species == "":
            species = "unknown"

        if build == "":
            build = "unknown"

        if version == "":
            version = "unknown"
       
        context = {
            "date": date,
            "source": "TraFiC-MEMv1.0",
            "reference": version
        }

        ### 2. Generate header
        ## 2.1 Header chunk 1
        template = """##fileformat=VCFv4.2
##fileDate={date}
##source={source}
##reference={reference}"""

        chunk1 = template.format(**context)
        
        ## 2.2 Header chunk 2
        chromList = genomeObj.fastaDict.keys()
        chromListSorted = sorted(chromList)

        chromList = []

        ## Create header chromosome rows one by one
        for chrom in chromListSorted:

            context = {
                "chrom": chrom,
                "build": build,
                "length": len(genomeObj.fastaDict[chrom]),
                "species": species
            }
            
            template = """
##contig=<ID={chrom},assembly={build},length={length},species={species}>"""
            chromRow = template.format(**context)
            chromList.append(chromRow)

        chunk2 = "".join(chromList)

        ## 2.3 Header chunk 3
        chunk3 = """
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Insertion type (TD0: solo, TD1: partnered-3'transduction, TD2: orphan-3'transduction, PSD: processed-pseudogene)>
##INFO=<ID=CLASS,Number=1,Type=String,Description="Mobile element family (L1, ALU, SVA or ERVK)">
##INFO=<ID=SUBFAMILY,Number=1,Type=String,Description="Mobile element subfamily">
##INFO=<ID=PDIV,Number=1,Type=Float,Description="Percentage divergence with respect to consensus mobile element subfamily">
##INFO=<ID=MECHANISM,Number=1,Type=String,Description="Insertion mechanism (TPRT: target primed reverse transcription, EI: endonuclease independent, DPA: double poly-A)>
##INFO=<ID=GR,Number=1,Type=String,Description="L1-mediated genomic rearrangement (DEL: deletion, DUP: duplication or TRANS: translocation)">
##INFO=<ID=MANUAL,Number=0,Type=Flag,Description="MEI manually verified through BAM inspection">
##INFO=<ID=SCORE,Number=1,Type=Integer,Description="Insertion score (5: 5' and 3' breakpoints (bkp) identified, 4: 3'bkp identified, 3: 5'bkp identified, 2: no bkp identified, 1: inconsistent (contradictory insertion features))">
##INFO=<ID=BKPB,Number=1,Type=Integer,Description="MEI right-most breakpoint position (bkpB). Left-most breakpoint position (bkpA) represented in the POS field">
##INFO=<ID=CIPOS,Number=1,Type=Integer,Description="Confidence interval around insertion breakpoint">
##INFO=<ID=CCO,Number=1,Type=String,Description="Clipped read clusters orientation: beg (clipping at the read begins) and end (clipping at the read ends). First orientation corresponds to the bkpA, while second to the bkpB. >
##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion DNA strand (+ or -)">
##INFO=<ID=STRUCT,Number=1,Type=String,Description="Mobile element structure (INV: 5'inverted, DEL: 5'deleted, FULL: full-length)">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Mobile element length">
##INFO=<ID=RANGE,Number=1,Type=Integer,Description="Begin and end position of the piece of mobile element integrated">
##INFO=<ID=TSLEN,Number=1,Type=Integer,Description="Target site duplication (+_value) or deletion (-_value) length">
##INFO=<ID=SRCID,Number=1,Type=String,Description="Germline source element cytoband identifier">
##INFO=<ID=SRCTYPE,Number=1,Type=String,Description="Source element type (GERMLINE or SOMATIC)">
##INFO=<ID=SRC,Number=1,Type=String,Description="Coordinates of the source element that mediated the transduction in the format: chrom_beg_end">
##INFO=<ID=TDC,Number=1,Type=String,Description="Begin and end coordinates of the integrated transduced or pseudogene sequence in the format: chrom_beg_end">
##INFO=<ID=TDLEN,Number=1,Type=Integer,Description="Transduced region length">
##INFO=<ID=TDLENR,Number=1,Type=Integer,Description="Transduced region length at RNA level">
##INFO=<ID=SRCGENE,Number=1,Type=String,Description="Source gene of the processed pseudogene insertion">
##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region where the mobile element is inserted (exonic, splicing, ncRNA, UTR5, UTR3, intronic, upstream, downstream, intergenic)">
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO gene symbol">
##INFO=<ID=ROLE,Number=1,Type=String,Description="Role in cancer (oncogene, TSG: tumor suppressor gene, oncogene/TSG: both roles)">
##INFO=<ID=COSMIC,Number=0,Type=Flag,Description="Reported as cancer driver in COSMIC cancer gene census database">
##INFO=<ID=CPG,Number=0,Type=Flag,Description="Reported as cancer predisposition gene in 10.1038/nature12981 (DOI).">
##INFO=<ID=REP,Number=1,Type=String,Description="Repetitive element overlapping the insertion breakpoints">
##INFO=<ID=DIV,Number=1,Type=Integer,Description="Millidivergence of the overlapping repetitive element with respect a consensus sequence">
##INFO=<ID=GERMDB,Number=.,Type=String,Description="List of sources that previously reported this event">
##INFO=<ID=IDS,Number=.,Type=String,Description="List of identifiers assigned to this event. Same order as in GERMDB">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Assembled contig sequence spanning bkpA">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Assembled contig sequence spanning bkpB">
##INFO=<ID=DP,Number=.,Type=String,Description="Discordant read-pairs Positive cluster">
##INFO=<ID=DN,Number=.,Type=String,Description="Discordant read-pairs Negative cluster">
##INFO=<ID=CA,Number=.,Type=String,Description="Clipped reads breakpoint A">
##INFO=<ID=CB,Number=.,Type=String,Description="Clipped reads breakpoint B">
##FILTER=<ID=SCORE,Description="Insertion with an score < threshold">
##FILTER=<ID=FPSOURCE,Description="Transduction mediated by a false somatic source element">
##FILTER=<ID=REP,Description="Insertion overlapping a satellite region or a repetitive element of the same class">
##FILTER=<ID=DUP,Description="Duplicated MEI call">
##FILTER=<ID=GERMLINE,Description="Germline MEI miscalled as somatic">
##FILTER=<ID=TD,Description="L1 transduction incorrectly identified as a processed pseudogene insertion">
##FILTER=<ID=COUNT,Description="Allele count < threshold">
##FILTER=<ID=MISSGT,Description="Genotype missing-ness rate > threshold ">
##FORMAT=<ID=NDP,Number=1,Type=Integer,Description="Number of Discordant read-pairs Positive cluster">
##FORMAT=<ID=NDN,Number=1,Type=Integer,Description="Number of Discordant read-pairs Negative cluster">
##FORMAT=<ID=NCA,Number=1,Type=Integer,Description="Number of Clipped reads breakpoint A">
##FORMAT=<ID=NCB,Number=1,Type=Integer,Description="Number of Clipped reads breakpoint B">
##FORMAT=<ID=SL,Number=1,Type=String,Description="List of samples where the variant was found">
##FORMAT=<ID=REPR,Number=1,Type=String,Description="Representative sample among all the samples where the variant is identified">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting the variant">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads supporting the reference allele">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
"""

        ## Create complete header:
        self.header = ''.join([chunk1, chunk2, chunk3])


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
        infoOrder = ['SVTYPE', 'TYPE', 'CLASS', 'SUBFAMILY', 'PDIV', 'MECHANISM', 'GR', 'MANUAL', 'SCORE', 'BKPB', 'CIPOS', 'CCO', 'STRAND', 'STRUCT', 'LEN', 'RANGE', 'TSLEN', 'SRCID', 'SRCTYPE', 'SRC', 'TDC', 'TDLEN', 'TDLENR', 'SRCGENE', 'REGION', 'GENE', 'ROLE', 'COSMIC', 'CPG', 'REP', 'DIV', 'GERMDB', 'IDS', 'AC', 'AN', 'AF', 'NS', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'CONTIGA', 'CONTIGB', 'DP', 'DN', 'CA', 'CB']

        flagList = ['COSMIC', 'CPG' , 'MANUAL']

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

    def sequences_lenght(self):
        """
        """
        

        

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
