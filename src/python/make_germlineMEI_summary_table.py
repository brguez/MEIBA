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
parser.add_argument('vcf', help='')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
vcf = args.vcf
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", vcf
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## 2. PROCESS VCF FILES
########################
outPath = outDir + '/' + fileName + '.tsv'

outFile = open(outPath, 'w')
         
## Write file header in the output file
fieldList = ["chrom", "firstBkp", "secondBkp", "insertionType", "family", "score", "mechanism", "structure", "rtLen", "strand", "tsLen", "rep", "div", "region", "gene", "role", "cosmic", 'novelty', 'germDb', 'ids', 'AC', 'AN', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', "firstContig", "secondContig", "\n"]
row = "\t".join(fieldList)
outFile.write(row)


# Read VCF and add information to VCF object
VCFObj = formats.VCF()
VCFObj.read_VCF(vcf)
        
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
            
        ## Support
        score = MEIObj.infoDict["SCORE"] if "SCORE" in MEIObj.infoDict else 'NA'
            
        ## Insertion features
        mechanism = MEIObj.infoDict["MECHANISM"] if "MECHANISM" in MEIObj.infoDict else 'NA'
        structure = MEIObj.infoDict["STRUCT"] if "STRUCT" in MEIObj.infoDict else 'NA'
        rtLen = MEIObj.infoDict["LEN"] if "LEN" in MEIObj.infoDict else 'NA'
        strand = MEIObj.infoDict["STRAND"] if "STRAND" in MEIObj.infoDict else 'NA'
        tsLen = MEIObj.infoDict["TSLEN"] if "TSLEN" in MEIObj.infoDict else 'NA'

        ## MEI annotation
        rep = MEIObj.infoDict["REP"] if "REP" in MEIObj.infoDict else 'NA'
        div = MEIObj.infoDict["DIV"] if "DIV" in MEIObj.infoDict else 'NA'
        region = MEIObj.infoDict["REGION"] if "REGION" in MEIObj.infoDict else 'NA'
        gene = MEIObj.infoDict["GENE"] if "GENE" in MEIObj.infoDict else 'NA'
        role = MEIObj.infoDict["ROLE"] if "ROLE" in MEIObj.infoDict else 'NA'
        cosmic = MEIObj.infoDict["COSMIC"] if "COSMIC" in MEIObj.infoDict else 'NA'

        ## Population stats
        germDb = MEIObj.infoDict["GERMDB"] if "GERMDB" in MEIObj.infoDict else 'NA'
        novelty = 'NOVEL' if germDb == 'NA' else 'KNOWN'
        ids = MEIObj.infoDict["IDS"] if "IDS" in MEIObj.infoDict else 'NA'
        AC = MEIObj.infoDict["AC"] if "AC" in MEIObj.infoDict else 'NA'
        AN = MEIObj.infoDict["AN"] if "AN" in MEIObj.infoDict else 'NA'
        AF = MEIObj.infoDict["AF"] if "AF" in MEIObj.infoDict else 'NA'
        AFR_AF = MEIObj.infoDict["AFR_AF"] if "AFR_AF" in MEIObj.infoDict else 'NA'
        AMR_AF = MEIObj.infoDict["AMR_AF"] if "AMR_AF" in MEIObj.infoDict else 'NA'
        EAS_AF = MEIObj.infoDict["EAS_AF"] if "EAS_AF" in MEIObj.infoDict else 'NA'
        EUR_AF = MEIObj.infoDict["EUR_AF"] if "EUR_AF" in MEIObj.infoDict else 'NA'
        SAS_AF = MEIObj.infoDict["SAS_AF"] if "SAS_AF" in MEIObj.infoDict else 'NA'

        ## Contigs
        contigA = MEIObj.infoDict["CONTIGA"] if "CONTIGA" in MEIObj.infoDict else 'NA'
        contigB = MEIObj.infoDict["CONTIGB"] if "CONTIGB" in MEIObj.infoDict else 'NA'      

        ### 2. Reorder breakpoints and contigs to have them in the same order as in the tumour genome 
        ###############################################################################################
        ## a) Both breakpoints reconstructed and no target site or target site duplication
        if ((MEIObj.infoDict["SCORE"] == "1") or (MEIObj.infoDict["SCORE"] == "5")) and (int(tsLen) >= 0):
            firstBkp = bkpB
            secondBkp = bkpA
            firstContig = contigB
            secondContig = contigA
           
        ## b) At least one breakpoint not reconstructed or target site deletion
        else:
            firstBkp = bkpA
            secondBkp = bkpB
            firstContig = contigA
            secondContig = contigB

        ### 3. Write MEI into the output file
        ######################################
        fieldList = [chrom, firstBkp, secondBkp, insertionType, family, score, mechanism, structure, rtLen, strand, tsLen, rep, div, region, gene, role, cosmic, novelty, germDb, ids, AC, AN, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, firstContig, secondContig, "\n"]
        row = "\t".join(fieldList)
        outFile.write(row)

# Close output file:
outFile.close()

## Finish ##
print
print "***** Finished! *****"
print



