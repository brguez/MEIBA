#!/usr/bin/env python
#coding: utf-8 

def header(string):
    """ 
        Display  header
    """ 
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats

## Get user's input ## 
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('--min-score', default=5, dest='minScore', type=int, help='Minimum insertion score for L1, Alu and SVA. Default: 5, both breakpoints expected to be assembled (5-prime and 3-prime).' )
parser.add_argument('--min-score-ERVK', default=3, dest='minScoreERVK', type=int, help='Minimum insertion score for ERVK. Note that mechanism of retrotransposition different, no polyA, so an ERVK insertion can not have an score > 3. Default: 3.' )
parser.add_argument('--max-divergence', default=300, dest='maxDiv', type=int, help='Maximum millidivergence. Default: 300.' )
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
donorId = args.donorId
minScore = args.minScore
minScoreERVK = args.minScoreERVK
maxDiv = args.maxDiv
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "donorId: ", donorId
print "minScore: ", minScore
print "minScoreERVK: ", minScoreERVK
print "maxDiv: ", maxDiv
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 

## Start ## 

outFilePath = outDir + '/' + donorId + ".filtered.vcf"

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## 2. Filter MEI

# Iterate over each MEI in the VCF
for VCFlineObj in VCFObj.lineList:

    ## 2.1 Apply score filter:
    # score < minScore_threshold -> filter out
    
    # A) L1, Alu or SVA insertion
    if (VCFlineObj.infoDict["CLASS"] != "ERVK") and (int(VCFlineObj.infoDict["SCORE"]) < minScore):
	VCFlineObj.filter = "SCORE"

    # B) ERVK insertion
    elif (VCFlineObj.infoDict["CLASS"] == "ERVK") and (int(VCFlineObj.infoDict["SCORE"]) < minScoreERVK):
	VCFlineObj.filter = "SCORE"

    ## 2.2 Repeats filter:
    # Filter out those MEI that overlap a repetitive element
    # Plus at least one of these two conditions:
    # 1. Insertion overlapping a repetitive element of the same class with a millidivergence < maxDiv_threshold. Notes: 
    #    - Millidivergence is a measure of the degree betweene a given repetitive element and a consensus reference sequence.
    #    - Goes from 1 (0.1%) to 1000 (100%) 
    #    - It is the opposite of sequence identity. 
    #    - Elements with low divergence to consensus sequence are highly repetitive in our genome and thus, prone to missalignments and 
    #     ,in last term, false positive MEI calls	
    # 2. Insertion overlapping one of the possible satellite regions in the database (ALR/Alpha, BSR/Beta, HSATII)

    # MEI overlaps a repetitive element
    if ('REP' in VCFlineObj.infoDict):
	
	## Condition 1 or 2 fulfilled:
    	if ((VCFlineObj.infoDict["CLASS"] == VCFlineObj.infoDict["REP"]) and (int(VCFlineObj.infoDict["DIV"]) < maxDiv)) or ((VCFlineObj.infoDict["REP"] == "ALR/Alpha") or (VCFlineObj.infoDict["REP"] == "BSR/Beta") or (VCFlineObj.infoDict["REP"] == "HSATII")):

	    # a) First filter failed  -> substitute . by filtering reason 
	    if (VCFlineObj.filter == "."):
		VCFlineObj.filter = "REP"
	    
	    # b) Already filtered by score -> append filtering reason 
	    else:
		VCFlineObj.filter = VCFlineObj.filter + ";REP"

    ## 2.3 MEI passed all the filters: 
    if (VCFlineObj.filter == "."):
	VCFlineObj.filter = "PASS"
	
## 3. Make output VCF 

# 3.1 Write header
VCFObj.write_header(outFilePath)

# 3.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print 
print "***** Finished! *****"
print 


