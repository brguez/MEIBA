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

outFilePath = outDir + '/' + donorId + ".repeatAnnot.vcf"

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## 2. Remove duplicities. For each insertion select a single overlapping repeat
repeatAnnotFile = open(repeatAnnot, 'r')
counter = 0

repeatDict = {}

for line in repeatAnnotFile:

    line = line.rstrip('\n')

    lineList = line.split('\t')

    chrom = lineList[0]
    beg = lineList[1]
    end = lineList[2]
    category = lineList[3]
    REP = lineList[4]
    DIV = lineList[5]

    insertionId = category + ":" + chrom + "_" + beg + "_" + end

    annotDict = {}
    annotDict["REP"] = REP

    # Divergence applicable, not applicable (na) only for satellite repeats
    if (DIV != "na"):
        annotDict["DIV"] = DIV

    # A) First annotation appearance for a given insertion
    if insertionId not in repeatDict:

        repeatDict[insertionId] = annotDict

    # B) Already annotation appearance for a given insertion
    else:

        # B.1) Previous satellite region
        if (repeatDict[insertionId]["REP"] in ["ALR/Alpha", "BSR/Beta", "HSATII"]):
            repeatDict[insertionId]["REP"] = repeatDict[insertionId]["REP"]

        # B.2) Overlapping satellite region or repeat of the same family while previus not
        elif (DIV == "na") or ((category == REP) and (category != repeatDict[insertionId]["REP"])):
            repeatDict[insertionId] = annotDict

        # B.3) Overlapping repeat with lower millidivergence than previous one
        elif (int(DIV) < int(repeatDict[insertionId]["DIV"])):
            repeatDict[insertionId] = annotDict



## 3. Add repeat annotation to VCF file
for VCFlineObj in VCFObj.lineList:

    # Skip repeat annotation if pseudogene insertion. Repeat annotation only applicable to L1, Alu, SVA and ERVK? insertions
    if (VCFlineObj.infoDict["TYPE"] == "PSD") or (VCFlineObj.infoDict["TYPE"] == "DEL"):
        print "[WARNING] Skip repeat annotation since " + VCFlineObj.infoDict["TYPE"] + " insertion"
        continue

    # a) bkpB identified
    if ("BKPB" in VCFlineObj.infoDict):

        # Compute expected beg and end
        exBeg = VCFlineObj.pos - 150
        exEnd = int(VCFlineObj.infoDict["BKPB"]) + 150

    # b) bkpB not identified
    else:

        # Compute expected beg and end
        exBeg = VCFlineObj.pos - 150
        exEnd = VCFlineObj.pos + 150

    insertionId = VCFlineObj.infoDict["CLASS"] + ":" + VCFlineObj.chrom + "_" + str(exBeg) + "_" + str(exEnd)

    ## Insertion overlapping a repetitive sequence
    if (repeatDict[insertionId]["REP"] != "."):

        # Add overlapping repeat information to info
        VCFlineObj.infoDict.update(repeatDict[insertionId])
        VCFlineObj.info = VCFlineObj.make_info()

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
