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
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## PROCESS VCF FILES
#####################
outPath = outDir + '/' + fileName + '.tsv'

outFile = open(outPath, 'w')
         
## Write file header in the output file
row = 'chrom' + '\t' + 'pos' + '\t' + "icgc_donor_id" + "\t" + "histology_abbreviation" + "\t" + 'class' + '\t' + 'type' + '\t' + 'strand' + '\t' + 'region' + '\t' + 'gene' + '\t' + 'mechanism' + '\n'
    
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
            
            ## Select only those MEI with SCORE 5 or 4 that pass all the filters
            # Exclude pseudogenes and rearrangements
            if (MEIObj.filter == "PASS") and ((MEIObj.infoDict['SCORE'] == '4') or (MEIObj.infoDict['SCORE'] == '5')) and (MEIObj.infoDict['TYPE'] != "PSD") and ("GR" not in MEIObj.infoDict):
                chrom = MEIObj.chrom
                bkpA = int(MEIObj.pos)
                bkpB = int(MEIObj.infoDict["BKPB"]) if "BKPB" in MEIObj.infoDict else 'NA'
                contigA = MEIObj.infoDict['CONTIGA'] if 'CONTIGA' in MEIObj.infoDict else 'UNK' 
                contigB = MEIObj.infoDict['CONTIGB'] if 'CONTIGB' in MEIObj.infoDict else 'UNK' 
                rtClass = MEIObj.infoDict['CLASS'] 
                rtType = MEIObj.infoDict['TYPE']  
                strand = MEIObj.infoDict['STRAND'] if 'STRAND' in MEIObj.infoDict else 'UNK'               
                region = MEIObj.infoDict['REGION'] 
                gene = MEIObj.infoDict['GENE'] if 'GENE' in MEIObj.infoDict else 'NA' 
                tsLen = int(MEIObj.infoDict["TSLEN"]) if "TSLEN" in MEIObj.infoDict else 'NA'
                mechanism = MEIObj.infoDict['MECHANISM'] if 'MECHANISM' in MEIObj.infoDict else 'UNK' 

                ### Select 3' breakpoint position
                # A) Only 3' bkp reconstructed                
                if (MEIObj.infoDict['SCORE'] == '4'):
                    pos = bkpA
                    contig = contigA

                # B) Both bkp reconstructed
                else:
                    ## a) First bkp polyA/3': [Target site duplication + positive strand] OR [Target site deletion + negative strand]
                    if ((tsLen >= 0) and (strand == "+")) or ((tsLen < 0) and (strand == "-")):
                        pos = bkpA
                        contig = contigA
            
                    ## b) Second bkp polyA/3': [Target site duplication + negative strand] OR [Target site deletion + positive strand]
                    else:
                        pos = bkpB
                        contig = contigB
                
                ## Write MEI into the output file
                row = chrom + '\t' + str(pos) + '\t' + icgcDonorId + '\t' + tumorType + '\t' + rtClass + '\t' + rtType + '\t' + strand + '\t' + region + '\t' + gene + '\t' + mechanism + '\n'
                outFile.write(row)
 

## End ##
print
print "***** Finished! *****"
print
