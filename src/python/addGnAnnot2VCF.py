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
parser.add_argument('annovarOut', help='...')
parser.add_argument('donorId', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputVCF = args.VCF
annovarOut = args.annovarOut
donorId = args.donorId
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "annovar-output: ", annovarOut
print "donorId: ", donorId
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 


## Start ## 

outFilePath = outDir + '/' + donorId + ".gnAnnot.vcf"

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## 2. Create Annovar object and read input annovar file
annovarObj = formats.annovar()
annovarObj.read_annovar(annovarOut)

## 3. Add annovar annotation information to output VCF file
annovarObj.info2VCF(VCFObj)

## 4. Make output VCF 

# 4.1 Write header
VCFObj.write_header(outFilePath)

# 4.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print 
print "***** Finished! *****"
print 


