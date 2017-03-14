#!/usr/bin/env python
#coding: utf-8

import re

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
        infoOrder = [ "SVTYPE", "CLASS", "TYPE", "SCORE", "MANUAL", "BKPB", "CIPOS", "STRAND", "STRUCT", "LEN", "TSLEN", "TSSEQ", "POLYA", "SRCID", "SRCTYPE", "SRC", "TDC", "TDLEN", "TDLENR", "SRCGENE", "GERMDB", "REGION", "GENE", "ROLE", "COSMIC", "CPG", "REP", "DIV", "CONTIGA", "CONTIGB", "RP", "RN" ]

        flagList = ["COSMIC", "CPG" , "MANUAL"]

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
