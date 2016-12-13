#!/usr/bin/env python
#coding: utf-8



#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import time
import numpy as np
import pandas as pd

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputFile', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputFile = args.inputFile
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputFile: ", inputFile
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

# Load table into dataframe
df = pd.read_csv(inputFile, header=0, sep='\t')

df = df.set_index('donor_id')

# Traspose dataframe
df=df.T 

print df

# Save output into tsv
outFilePath = outDir + '/germline_source_element_transductions_perDonor.tsv'
df.to_csv(outFilePath, sep='\t') 

####
