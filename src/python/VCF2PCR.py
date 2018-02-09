#!/usr/bin/env python
#coding: utf-8

#### FUNCTIONS ####
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

def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string


#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import os.path
import formats
import os.path
from operator import attrgetter

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('VCF', help='VCF file to be filtered')
parser.add_argument('genome', help='reference genome fasta file')
parser.add_argument('donorId', help='donor identifier')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
genome = args.genome
donorId = args.donorId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", inputVCF
print "genome: ", genome
print "donorId: ", donorId
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### Create reference genome fasta object
header("Creating reference genome fasta object")
genomeObj = formats.fasta(genome)

#### Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

#### Produce output file with sequences for PCR validation
outFilePath = outDir + '/' + donorId + ".sequences4PCR.tsv"
outFile = open(outFilePath, 'w')

# Write header 
row = "#chrom" + "\t" + "bkp5prime" + "\t" + "newContig5prime" + "\t" + "bkp3prime" + "\t" + "newContig3prime" + "\t" + "rtClass" + "\t" + "insertionType" + "\t" + "score" + "\t" + "tsLen" + "\t" + "strand" + "\t" + "structure" + "\t" + "length" + "\n"
outFile.write(row)

for VCFlineObj in VCFObj.lineList:

    
    ### 1. Set variables of interest        
    chrom = VCFlineObj.chrom
    bkpA = int(VCFlineObj.pos)
    bkpB = int(VCFlineObj.infoDict["BKPB"]) if "BKPB" in VCFlineObj.infoDict else 'NA'
    rtClass = VCFlineObj.infoDict["CLASS"] 
    insertionType = VCFlineObj.infoDict["TYPE"]
    mechanism = VCFlineObj.infoDict["MECHANISM"] if "MECHANISM" in VCFlineObj.infoDict else 'NA'
    score = VCFlineObj.infoDict["SCORE"]
    structure = VCFlineObj.infoDict["STRUCT"] if "STRUCT" in VCFlineObj.infoDict else 'NA'
    length = VCFlineObj.infoDict["LEN"] if "LEN" in VCFlineObj.infoDict else 'NA'
    strand = VCFlineObj.infoDict["STRAND"] if "STRAND" in VCFlineObj.infoDict else 'NA'
    tsLen = int(VCFlineObj.infoDict["TSLEN"]) if "TSLEN" in VCFlineObj.infoDict else 'NA'
    contigA = VCFlineObj.infoDict["CONTIGA"] if "CONTIGA" in VCFlineObj.infoDict else 'NA'
    contigB = VCFlineObj.infoDict["CONTIGB"] if "CONTIGB" in VCFlineObj.infoDict else 'NA'      

    print "*** MEI: ", chrom, bkpA, bkpB, rtClass, score, structure, length, strand, tsLen, contigA, contigB

    ### 2. Reorder breakpoints and contigs to have the 5' before the 3' 
    ## A) Both breakpoints reconstructed and no target site or target site duplication
    if ((VCFlineObj.infoDict["SCORE"] == "1") or (VCFlineObj.infoDict["SCORE"] == "5")) and (tsLen >= 0):
        bkp5prime = bkpB
        bkp3prime = bkpA
        contig5prime = contigB
        contig3prime = contigA
            
    ## B) At least one breakpoint not reconstructed or target site deletion
    else:
        bkp5prime = bkpA
        bkp3prime = bkpB
        contig5prime = contigA
        contig3prime = contigB
            
    ### 3. Expand target region sequence for designing the primers
    # A) No breakpoint reconstructed
    if (VCFlineObj.infoDict["SCORE"] == "2"):
        newContig5prime = "NA"
        newContig3prime = "NA"

    # B) Both breakpoints reconstructed
    elif ((VCFlineObj.infoDict["SCORE"] == "1") or (VCFlineObj.infoDict["SCORE"] == "5")):
        
        ## 5 prime
        seqA5prime = genomeObj.fastaDict[chrom][bkp5prime-500 : bkp5prime]
        seqB5prime = contig5prime.split("[MEI]>")[1]
        newContig5prime = seqA5prime + '[MEI]>' + seqB5prime
    
        ## 3 prime
        seqB3prime = genomeObj.fastaDict[chrom][bkp3prime : bkp3prime+500]    
        seqA3prime = contig3prime.split("<[MEI]")[0]
        newContig3prime = seqA3prime + '<[MEI]' + seqB3prime

    # C) One breakpoint reconstructed 
    else:

        ## a) 5 prime bkp reconstructed
        if len(contig5prime.split("[MEI]>")) == 2:

            seqA5prime = genomeObj.fastaDict[chrom][bkp5prime-500 : bkp5prime]
            seqB5prime = contig5prime.split("[MEI]>")[1]
            newContig5prime = seqA5prime + '[MEI]>' + seqB5prime
            newContig3prime = "NA"

        ## b) 3 prime bkp reconstructed
        else: 

            ## Switch the breakpoints 
            bkp3prime = bkp5prime
            bkp5prime = "NA"
            contig3prime = contig5prime
            contig5prime = "NA"

            ## Extend 3 prime contig genomic sequence
            seqB3prime = genomeObj.fastaDict[chrom][bkp3prime : bkp3prime+500]    
            seqA3prime = contig3prime.split("<[MEI]")[0]
            newContig3prime = seqA3prime + '<[MEI]' + seqB3prime
            newContig5prime = "NA"            

    ## Write row into the output file
    row = contigA + "\t" + contigB + "\t" + chrom + "\t" + str(bkp5prime) + "\t" + newContig5prime + "\t" + str(bkp3prime) + "\t" + newContig3prime + "\t" + insertionType + "\t" + rtClass + "\t" + score + "\t" + str(tsLen) + "\t" + strand + "\t" + structure + "\t" + length + "\n"
    outFile.write(row)


# Close output file:
outFile.close()

## Finish ##
print
print "***** Finished! *****"
print

  



