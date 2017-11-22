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


def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


#### MAIN ####
## Import modules ##
import argparse
import sys
import os.path
import time

## Get user's input ##
parser = argparse.ArgumentParser(description="Add genomic coordinates to COSMIC cancer genes")
parser.add_argument('cosmic', help='Cosmic cancer genes')
parser.add_argument('annot', help='Bed file containing protein coding gene coordinates')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
cosmic = args.cosmic
annot = args.annot
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "cosmic: ", cosmic
print "annot: ", annot
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

## 1. Read annotation file
###########################
# Initialize a dictionary with the following structure:
# - dict: key(geneName) -> geneCoordinates

header("1. Read annotation file")

annotFile = open(annot, 'r')
geneCoordDict = {}

for line in annotFile:

    # Skip header
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        chrom = line[0]
        beg = line[1]
        end = line[2]
        geneName = line[3]        
        geneCoordDict[geneName] = chrom + "_" + beg + "_" + end 


## 2. Incorporate gene coordinates into the cosmic file
########################################################

header("2. Incorporate gene coordinates into the cosmic file")

cosmicFile = open(cosmic, 'r')
outPath = outDir + '/cancerGenes_COSMIC_CPG.gnCoord.tsv'
outFile = open(outPath, 'w')

## Write file header in the output file
row = "## 20171115" + "\n" + "## Database of cancer genes. Composed by:" + "\n" + "## 1. Cancer driver genes from cosmic cancer gene census database (27/06/2016, COSMIC v77)" + "\n" + "## access_info: http://cancer.sanger.ac.uk/census/" + "\n" + "## 2. Cancer predisposition genes (CPG) from Rahman N. (doi:10.1038/nature12981)" + "\n" + "## access_info: Supplementary table 1" + "\n" + "## <CHROM>  <BEG>   <END>   <GENE_SYMBOL>	    <ROLE>  <COSMIC>    	<CPG>" + "\n"	
outFile.write(row)

for line in cosmicFile:

    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split('\t')

        geneName, role, COSMIC, CPG = line
        
        ## Coordinates available for the gene
        if geneName in geneCoordDict:
            chrom, beg, end = geneCoordDict[geneName].split("_")

        ## Coordinates NOT available for the gene
        else:
            print "Coordinates NOT available for: ",  geneName
            chrom, beg, end = ["UNK", "UNK", "UNK"]

        ## write line into the output file
        row = chrom + "\t" + beg + "\t" + end + "\t" + geneName + "\t" + role + "\t" + COSMIC + "\t" + CPG + "\n"    
        outFile.write(row)

####
header("Finished")

