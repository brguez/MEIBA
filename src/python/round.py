
#!/usr/bin/env python
#coding: utf-8

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import scipy.stats as stats
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputFile', help='tsv')
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


#### 1. Load input table:
##########################
dataframe = pd.read_csv(inputFile, header=0, sep='\t')

dataframe = dataframe.round(3)

dataframe.fillna("UNK", inplace=True)

#### 2. Write output table:
############################
outFilePath = outDir + '/supplementaryTable2.tsv'
dataframe.to_csv(outFilePath, index=False, sep='\t') 


print "***** Finished! *****"
print

