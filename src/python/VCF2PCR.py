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


def rev_complement(seq):
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
row = "#chrom" + "\t" + "firstBkp" + "\t" + "firstContig" + "\t" + "secondBkp" + "\t" + "secondContig" + "\t" + "rtClass" + "\t" + "insertionType" + "\t" + "mechanism" + "\t" + "filter" + "\t" + "score" + "\t" + "tsLen" + "\t" + "strand" + "\t" + "structure" + "\t" + "length" + "\n"
outFile.write(row)

for VCFlineObj in VCFObj.lineList:

    ### 1. Set variables of interest    
    ##################################    
    chrom = VCFlineObj.chrom
    bkpA = int(VCFlineObj.pos)
    bkpB = int(VCFlineObj.infoDict["BKPB"]) if "BKPB" in VCFlineObj.infoDict else 'NA'
    rtClass = VCFlineObj.infoDict["CLASS"] 
    insertionType = VCFlineObj.infoDict["TYPE"]
    mechanism = VCFlineObj.infoDict["MECHANISM"] if "MECHANISM" in VCFlineObj.infoDict else 'NA'
    filterField = VCFlineObj.filter
    score = VCFlineObj.infoDict["SCORE"]
    structure = VCFlineObj.infoDict["STRUCT"] if "STRUCT" in VCFlineObj.infoDict else 'NA'
    length = VCFlineObj.infoDict["LEN"] if "LEN" in VCFlineObj.infoDict else 'NA'
    strand = VCFlineObj.infoDict["STRAND"] if "STRAND" in VCFlineObj.infoDict else 'NA'
    tsLen = int(VCFlineObj.infoDict["TSLEN"]) if "TSLEN" in VCFlineObj.infoDict else 'NA'
    contigA = VCFlineObj.infoDict["CONTIGA"] if "CONTIGA" in VCFlineObj.infoDict else 'NA'
    contigB = VCFlineObj.infoDict["CONTIGB"] if "CONTIGB" in VCFlineObj.infoDict else 'NA'      
    tdCoord = VCFlineObj.infoDict["TDC"] if "TDC" in VCFlineObj.infoDict else 'NA'

    print "*** MEI: ", chrom, bkpA, bkpB, rtClass, insertionType, mechanism, score, structure, length, strand, tsLen, contigA, contigB

    ### 2. Reorder breakpoints and contigs to have them in the same order as in the tumour genome 
    ###############################################################################################
    ## A) Both breakpoints reconstructed and no target site or target site duplication
    if ((VCFlineObj.infoDict["SCORE"] == "1") or (VCFlineObj.infoDict["SCORE"] == "5")) and (tsLen >= 0):
        firstBkp = bkpB
        secondBkp = bkpA
        firstContig = contigB
        secondContig = contigA
            
    ## B) At least one breakpoint not reconstructed or target site deletion
    else:
        firstBkp = bkpA
        secondBkp = bkpB
        firstContig = contigA
        secondContig = contigB
            
    ### 3. Expand bkp contig sequences for designing the primers
    ###############################################################
    # A) No breakpoint reconstructed
    if (VCFlineObj.infoDict["SCORE"] == "2"):
        newFirstContig = "NA"
        newSecondContig = "NA"

    # B) Both breakpoints reconstructed
    elif ((VCFlineObj.infoDict["SCORE"] == "1") or (VCFlineObj.infoDict["SCORE"] == "5")):
        
        ### 3.1 Extend external contig fragments
        firstContigSeqA = genomeObj.fastaDict[chrom][firstBkp-500 : firstBkp]
        secondContigSeqB = genomeObj.fastaDict[chrom][secondBkp : secondBkp+500]    

        ### 3.2 Extend internal fragments
        # a) Insertion in the plus strand:  --------##########AAAAAA----------
        if (strand == "+"):
            firstContigSeqB = firstContig.split("[MEI]>")[1]

            # Solo 
            if (insertionType == "TD0"):
                secondContigSeqA = secondContig.split("<[MEI]")[0]

            # transduction
            else:
                chromTd = tdCoord.split("_")[0]
                begTd = int(tdCoord.split("_")[1]) 
                endTd = int(tdCoord.split("_")[2])
                secondContigSeqA = genomeObj.fastaDict[chromTd][begTd : endTd] + "AAAAAAAAAAAAAAAAAAAAAAAA"

        # b) Insertion in the minus strand: --------TTTTTT##########----------
        else:
            secondContigSeqA = secondContig.split("<[MEI]")[0]

            # Solo 
            if (insertionType == "TD0"):
                firstContigSeqB = firstContig.split("[MEI]>")[1]

            # transduction
            else:
                chromTd = tdCoord.split("_")[0]
                begTd = int(tdCoord.split("_")[1]) 
                endTd = int(tdCoord.split("_")[2])
                firstContigSeqB = rev_complement(genomeObj.fastaDict[chromTd][begTd : endTd] + "AAAAAAAAAAAAAAAAAAAAAAAA")
 
        ## 3.3 Join start and end
        newFirstContig = firstContigSeqA + '[MEI]>' + firstContigSeqB
        newSecondContig = secondContigSeqA + '<[MEI]' + secondContigSeqB

    # C) One breakpoint reconstructed 
    else:

        ## a) First bkp reconstructed
        # ----------[MEI]
        if len(firstContig.split("[MEI]>")) == 2:

            ### 3.1. Extend external contig fragment
            firstContigSeqA = genomeObj.fastaDict[chrom][firstBkp-500 : firstBkp]

            ### 3.2 Extend internal fragment
            ## Solo integration or no poly-A bkp  
            # -----#####            
            if (insertionType == "TD0") or (score == "3"):
                firstContigSeqB = firstContig.split("[MEI]>")[1]

            ## Transduction and poly-A bkp
            # --------TTTTTT###         * I will assume the insertion has - orientation  
            else:
                chromTd = tdCoord.split("_")[0]
                begTd = int(tdCoord.split("_")[1]) 
                endTd = int(tdCoord.split("_")[2])
                firstContigSeqB = rev_complement(genomeObj.fastaDict[chromTd][begTd : endTd] + "AAAAAAAAAAAAAAAAAAAAAAAA")
    
            ## 3.3 Join start and end
            newFirstContig = firstContigSeqA + '[MEI]>' + firstContigSeqB
            newSecondContig = "NA"

        ## b) Second bkp reconstructed
        # [MEI]----------
        else: 

            ## 3.0 Switch the breakpoints 
            secondBkp = firstBkp
            firstBkp = "NA"
            secondContig = firstContig
            firstContig = "NA"

            ### 3.1. Extend external contig fragment
            secondContigSeqB = genomeObj.fastaDict[chrom][secondBkp : secondBkp+500]    

            ### 3.2 Extend internal fragment
            ## Solo integration or no poly-A bkp  
            # #####------            
            if (insertionType == "TD0") or (score == "3"):
                secondContigSeqA = secondContig.split("<[MEI]")[0]

            ## Transduction and poly-A bkp
            # #####AAAAAA--------    * I will assume the insertion has + orientation  
            else:
                chromTd = tdCoord.split("_")[0]
                begTd = int(tdCoord.split("_")[1]) 
                endTd = int(tdCoord.split("_")[2])
                secondContigSeqA = genomeObj.fastaDict[chromTd][begTd : endTd] + "AAAAAAAAAAAAAAAAAAAAAAAA"
    
            ## 3.3 Join start and end
            newFirstContig = "NA" 
            newSecondContig = secondContigSeqA + '<[MEI]' + secondContigSeqB


    ## 4. Write row into the output file
    ######################################
    row = chrom + "\t" + str(firstBkp) + "\t" + newFirstContig + "\t" + str(secondBkp) + "\t" + newSecondContig + "\t" + rtClass + "\t" + insertionType + "\t" + mechanism + "\t" + filterField + "\t" + score + "\t" + str(tsLen) + "\t" + strand + "\t" + structure + "\t" + length + "\n"
    outFile.write(row)


# Close output file:
outFile.close()

## Finish ##
print
print "***** Finished! *****"
print

  



