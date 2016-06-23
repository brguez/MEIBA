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
parser.add_argument('repeatAnnot', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
repeatAnnot = args.repeatAnnot
donorId = args.donorId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "repeat-annotation: ", repeatAnnot
print "donorId: ", donorId
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 

## Start ## 

outFilePath = outDir + '/' + donorId + ".annotated.vcf"

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## 2. Add repeat annotation to VCF file
repeatAnnotFile = open(repeatAnnot, 'r')

counter = 0

for line in repeatAnnotFile:
	line = line.rstrip('\n')
	lineList = line.split('\t')
	chrom = lineList[0] 
	beg = int(lineList[1])
	end = int(lineList[2])
	category = lineList[3]
	REP = lineList[4]
	DIV = lineList[5]
	

        VCFlineObj =  VCFObj.lineList[counter]

	# A) Target site identified
	if ( VCFlineObj.infoDict["SCORE"] == "1"):

		# Compute expected beg and end
		exBeg = VCFlineObj.pos - 1
		exEnd = VCFlineObj.pos + abs(int(VCFlineObj.infoDict["TSLEN"]))

	# B) Target site not identified
	else:

		# Compute expected beg and end
		exBeg = VCFlineObj.pos - int(VCFlineObj.infoDict["CIPOS"]) - 1 
		exEnd = VCFlineObj.pos + int(VCFlineObj.infoDict["CIPOS"])

	## Insertion overlapping a repetitive sequence
	if (REP != "."):
	
		# A) VCF line correctly paired with repeat annotation file
		if (VCFlineObj.chrom == chrom) and (beg == exBeg) and (end == exEnd):
		
			annotDict = {}
			annotDict["REP"] = REP

			# Divergence applicable, not applicable (na) only for satellite repeats
			if (DIV != "na"):		
		
				annotDict["DIV"] = DIV		
			
			VCFlineObj.infoDict.update(annotDict)
			VCFlineObj.info = VCFlineObj.make_info()

		# A) VCF line unpaired with repeat annotation file
	    	else:
			print "ERROR. Unpaired VCF and repetitive sequences"

	counter += 1	    


repeatAnnotFile.close()

## 3. Make output VCF 

# 3.1 Write header
VCFObj.write_header(outFilePath)

# 3.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print 
print "***** Finished! *****"
print 


