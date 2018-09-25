#!/usr/bin/env python
#coding: utf-8

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


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
from operator import itemgetter, attrgetter, methodcaller

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: tumorSpecimenId   icgcDonorId vcfPath')
parser.add_argument('metadata', help='PCAWG donor metadata')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
metadata = args.metadata
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "metadata: ", metadata
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## 1. READ METADATA FILE:
#########################
metadata = open(metadata, 'r')

conversionDict = {}

# Read file line by line
for line in metadata:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        fieldsList = line.split("\t")
        sampleIds = fieldsList[20].split(",")
        aliquotIds = fieldsList[21].split(",")
    
        for index, aliquotId in enumerate(aliquotIds):
            conversionDict[aliquotId] = sampleIds[index]

## 2. PROCESS VCF FILES
########################
outPath = outDir + '/' + fileName + '.tsv'

outFile = open(outPath, 'w')
         
## Write file header in the output file
fieldList = ["donorId", "sampleIds", "tumorType", "chrom", "firstBkp", "secondBkp", "insertionType", "family", "subfamily", "subfamDiv", "rg", "srcGene", "score", "discordantPlus", "discordantMinus", "clippedFirstBkp", "clippedSecondBkp", "structure", "rtLen", "strand", "tsLen", "srcId", "srcType", "srcCoord", "tdCoord", "tdLen", "rep", "div", "region", "gene", "role", "cosmic", "firstContig", "secondContig", "\n"]
row = "\t".join(fieldList)
outFile.write(row)

inputFile = open(inputPath, 'r')

# Per iteration, read a VCF
for line in inputFile:
    line = line.rstrip('\n')
    line = line.split("\t")

    icgcDonorId = line[0]
    tumorType = line[1]
    vcfPath = line[2]

    print "Processing: ", icgcDonorId, tumorType, vcfPath
    # Create VCF object
    VCFObj = formats.VCF()

    ## Raise error if the input VCF is not available
    if not os.path.isfile(vcfPath):

        print "[ERROR] Input file is not available"

    ## Input VCF available
    else:

        # Read VCF and add information to VCF object
        VCFObj.read_VCF(vcfPath)
        
        ## For each MEI
        for MEIObj in VCFObj.lineList:  
            
            ## Select only those MEI that pass all the filters
            if (MEIObj.filter == "PASS"):

                ### 1. Set variables of interest    
                ##################################    
                ## Basic
                chrom = MEIObj.chrom
                bkpA = str(MEIObj.pos)
                bkpB = MEIObj.infoDict["BKPB"] if "BKPB" in MEIObj.infoDict else 'NA'
                insertionType = MEIObj.infoDict["TYPE"] if "TYPE" in MEIObj.infoDict else 'NA'
                family = MEIObj.infoDict["CLASS"] if "CLASS" in MEIObj.infoDict else 'NA'
                subfamily = MEIObj.infoDict["SUBFAMILY"] if "SUBFAMILY" in MEIObj.infoDict else 'NA'
                subfamDiv = MEIObj.infoDict["PDIV"] if "PDIV" in MEIObj.infoDict else 'NA'
                rg = MEIObj.infoDict["GR"] if "GR" in MEIObj.infoDict else 'NA'
                srcGene = MEIObj.infoDict["SRCGENE"] if "SRCGENE" in MEIObj.infoDict else 'NA'

                ## Support
                score = MEIObj.infoDict["SCORE"] if "SCORE" in MEIObj.infoDict else 'NA'
                discordantPlus, discordantMinus, clippedA, clippedB, aliquotIds = MEIObj.genotype.split(":")    ## pending replace sampleIds terminology
                clippedA = "0" if clippedA == "." else clippedA
                clippedB = "0" if clippedB == "." else clippedB

                # convert aliquot into sample ids
                aliquotIdList = aliquotIds.split(";")
                sampleIdList = []

                for aliquotId in aliquotIdList:
                    sampleIdList.append(conversionDict[aliquotId])

                sampleIds = ",".join(sampleIdList)

                ## Insertion features
                structure = MEIObj.infoDict["STRUCT"] if "STRUCT" in MEIObj.infoDict else 'NA'
                rtLen = MEIObj.infoDict["LEN"] if "LEN" in MEIObj.infoDict else 'NA'
                strand = MEIObj.infoDict["STRAND"] if "STRAND" in MEIObj.infoDict else 'NA'
                tsLen = MEIObj.infoDict["TSLEN"] if "TSLEN" in MEIObj.infoDict else 'NA'
                srcId = MEIObj.infoDict["SRCID"] if "SRCID" in MEIObj.infoDict else 'NA'
                srcType = MEIObj.infoDict["SRCTYPE"] if "SRCTYPE" in MEIObj.infoDict else 'NA'
                srcCoord = MEIObj.infoDict["SRC"] if "SRC" in MEIObj.infoDict else 'NA'
                tdCoord = MEIObj.infoDict["TDC"] if "TDC" in MEIObj.infoDict else 'NA'
                tdLen = MEIObj.infoDict["TDLEN"] if "TDLEN" in MEIObj.infoDict else 'NA'
                mechanism = MEIObj.infoDict["MECHANISM"] if "MECHANISM" in MEIObj.infoDict else 'NA'

                ## MEI annotation
                rep = MEIObj.infoDict["REP"] if "REP" in MEIObj.infoDict else 'NA'
                div = MEIObj.infoDict["DIV"] if "DIV" in MEIObj.infoDict else 'NA'
                region = MEIObj.infoDict["REGION"] if "REGION" in MEIObj.infoDict else 'NA'
                gene = MEIObj.infoDict["GENE"] if "GENE" in MEIObj.infoDict else 'NA'
                role = MEIObj.infoDict["ROLE"] if "ROLE" in MEIObj.infoDict else 'NA'
                cosmic = MEIObj.infoDict["COSMIC"] if "COSMIC" in MEIObj.infoDict else 'NA'

                contigA = MEIObj.infoDict["CONTIGA"] if "CONTIGA" in MEIObj.infoDict else 'NA'
                contigB = MEIObj.infoDict["CONTIGB"] if "CONTIGB" in MEIObj.infoDict else 'NA'      



                ### 2. Reorder breakpoints and contigs to have them in the same order as in the tumour genome 
                ###############################################################################################
                ## a) Both breakpoints reconstructed and no target site or target site duplication
                if ((MEIObj.infoDict["SCORE"] == "1") or (MEIObj.infoDict["SCORE"] == "5")) and (int(tsLen) >= 0):
                    firstBkp = bkpB
                    secondBkp = bkpA
                    clippedFirstBkp = clippedB
                    clippedSecondBkp = clippedA
                    firstContig = contigB
                    secondContig = contigA
           
                ## b) At least one breakpoint not reconstructed or target site deletion
                else:
                    firstBkp = bkpA
                    secondBkp = bkpB
                    clippedFirstBkp = clippedA 
                    clippedSecondBkp = clippedB
                    firstContig = contigA
                    secondContig = contigB

                ### 3. Write MEI into the output file
                ######################################
                fieldList = [icgcDonorId, sampleIds, tumorType, chrom, firstBkp, secondBkp, insertionType, family, subfamily, subfamDiv, rg, srcGene, score, discordantPlus, discordantMinus, clippedFirstBkp, clippedSecondBkp, structure, rtLen, strand, tsLen, srcId, srcType, srcCoord, tdCoord, tdLen, rep, div, region, gene, role, cosmic, firstContig, secondContig, mechanism, "\n"]
                row = "\t".join(fieldList)
                outFile.write(row)


# Close output file:
outFile.close()

## Finish ##
print
print "***** Finished! *****"
print



