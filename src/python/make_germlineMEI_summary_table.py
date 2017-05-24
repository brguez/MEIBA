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
parser.add_argument('vcfPath', help='VCF with germline MEI calls')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
vcfPath = args.vcfPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", vcfPath
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print



## PROCESS VCF FILES
#####################
outPath = outDir + '/PCAWG_germlineRT_summaryTable_0.5.7.tsv'

outFile = open(outPath, 'w')

         
## Write file header in the output file
row = '#chrom' + '\t' + 'pos' + '\t' + 'class' + '\t' + 'type' + '\t' + 'score' + '\t' + 'cipos' + '\t' + 'length'  + '\t' + 'structure' + '\t' + 'strand' + '\t' + 'tsLen' + '\t' + 'tsSeq' + '\t' + 'polyA' + '\t' + 'region' + '\t' + 'gene' + '\n'
    
outFile.write(row)


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
            
        ## Select only those MEI that pass all the filters:
        if (MEIObj.filter == "PASS"):
            chrom = MEIObj.chrom
            pos = MEIObj.pos
            rtClass = MEIObj.infoDict['CLASS'] 
            rtType = MEIObj.infoDict['TYPE']  
            score = MEIObj.infoDict['SCORE']     
            cipos = MEIObj.infoDict['CIPOS']     
            length = MEIObj.infoDict['LEN'] if 'LEN' in MEIObj.infoDict else 'UNK' 
            structure = MEIObj.infoDict['STRUCT'] if 'STRUCT' in MEIObj.infoDict else 'UNK' 
            strand = MEIObj.infoDict['STRAND'] if 'STRAND' in MEIObj.infoDict else 'UNK' 
            tsLen = MEIObj.infoDict['TSLEN'] if 'TSLEN' in MEIObj.infoDict else 'UNK' 
            tsSeq = MEIObj.infoDict['TSSEQ'] if 'TSSEQ' in MEIObj.infoDict else 'UNK'
               
            # For those insertions with tsLen == 0, set tsSeq as NA
            if (tsLen == '0'):
                tsSeq = "NA"         
       
            polyA = MEIObj.infoDict['POLYA'] if 'POLYA' in MEIObj.infoDict else 'UNK' 
            region = MEIObj.infoDict['REGION'] 
            gene = MEIObj.infoDict['GENE'] if 'GENE' in MEIObj.infoDict else 'NA' 
                

            ## Write MEI into the output file
            row = chrom + '\t' + str(pos) + '\t' + rtClass + '\t' + rtType + '\t' + score + '\t' + cipos + '\t' + length  + '\t' + structure + '\t' + strand + '\t' + tsLen + '\t' + tsSeq + '\t' + polyA + '\t' + region + '\t' + gene + '\n'
            outFile.write(row)
 

## End ##
print
print "***** Finished! *****"
print
