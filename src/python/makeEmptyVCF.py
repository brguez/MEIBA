#!/usr/bin/env python
#coding: utf-8

## Import modules ##
import argparse
import time
import sys
import os
import formats

#### FUNCTIONS ####

def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"

def parse_args():
    """Define and parse command line parameters."""
    parser = argparse.ArgumentParser(description="Generate MEI VCF file only containing the header")
    parser.add_argument('donorId', help='Donor identifier. The output vcf will be named accordingly')
    parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

    args = parser.parse_args()
    return args

#### MAIN ####
if __name__ == "__main__":

    ## Get user's input ##
    args = parse_args()
    donorId = args.donorId
    outDir = args.outDir

    scriptName = os.path.basename(sys.argv[0])

    ## Display configuration to standard output ##
    print
    print "***** ", scriptName, " configuration *****"
    print "donorId: ", donorId
    print "outDir: ", outDir
    print
    print "***** Executing ", scriptName, " *****"
    print


    ## Start ## 
    outFilePath = outDir + '/' + donorId + '.vcf'

    ## 1. Create VCF object and print VCF header    
    header("Creating VCF object and printing VCF header into the output file")
    VCFObj = formats.VCF()
    VCFObj.create_header()

    VCFObj.write_header(outFilePath)

    ## Finish ##
    print
    print "***** Finished! *****"
    print
