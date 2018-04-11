#!/usr/bin/env python
#coding: utf-8

def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"

def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import time
import pandas as pd
import numpy as np


## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('rtCounts', help='')
parser.add_argument('scoreStats', help='')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
rtCounts = args.rtCounts
scoreStats = args.scoreStats
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "rtCounts: ", rtCounts
print "scoreStats: ", scoreStats
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ##

## 1. Read input files
###########################
rtCountsDf = pd.read_csv(rtCounts, header=0, index_col=0, sep='\t')
scoreStatsDf = pd.read_csv(scoreStats, header=0, index_col=0, sep='\t')


## 2. Incorporate the fraction of MEI with polyA bkp per sample
#################################################################
rtCountsDf["fractionPolyA"] = scoreStatsDf["fractionPolyA"]


## 3. Compute for each sample the ratio between the number of MEI 
#################################################################
# in the sample and the median number of MEI in the tumor type
#############################################################
tumorTypes = set(rtCountsDf['tumorType'].tolist())

nbTotalFcSeries = pd.Series()
 
## Per tumor type
for tumorType in tumorTypes:

    ## compute the median number of retrotransposition events per sample
    medianRt = rtCountsDf[rtCountsDf['tumorType'] == tumorType]['nbTotal'].median()
    medianRt = medianRt + 1 # Add pseudocount of 1

    ## compute for each sample of the current tumor type the ratio of nbMEI/medianNbMEI
    nbTotalFcSeries = pd.concat([nbTotalFcSeries, rtCountsDf[rtCountsDf['tumorType'] == tumorType]['nbTotal']/medianRt])

## Incorporate the fold change ratio into the dataframe
rtCountsDf["nbTotalFC"] = nbTotalFcSeries


## 4. Exclude outlier samples with a high percentage of 
########################################################
# likely false positive calls
##############################
## Exluded samples will fulfil both criteria:
# - Outlier: more than 5 times more events than the median number of events in the tumor type
# - Bad quality: less than 20% of the sample with polyA identified
rtCountsDf["Excluded"] = (rtCountsDf["fractionPolyA"] < 0.2) & (rtCountsDf["nbTotalFC"] > 5)


## 5. Report the number of insertions per category before  
##########################################################
# and after the excluding samples
#####################################

totalCounts = rtCountsDf.sum(axis=0)

print "*** Counts before excluding samples"
print "nbTotal: ", totalCounts["nbTotal"]
print "nbL1Solo: ", totalCounts["nbL1Solo"]
print "nbAlu: ", totalCounts["nbAlu"]
print "nbL1TD: ", totalCounts["nbL1TD"]
print "nbSVA: ", totalCounts["nbSVA"]
print "nbERVK: ", totalCounts["nbERVK"]
print "nbPSD: ", totalCounts["nbPSD"]
print "nbL1DEL: ", totalCounts["nbL1DEL"]
print "nbL1DUP: ", totalCounts["nbL1DUP"]
print 

rtCountsFilteredDf = rtCountsDf[rtCountsDf["Excluded"] == False] 
totalCounts = rtCountsFilteredDf.sum(axis=0)

print "*** Counts after excluding samples"
print "nbTotal: ", totalCounts["nbTotal"]
print "nbL1Solo: ", totalCounts["nbL1Solo"]
print "nbAlu: ", totalCounts["nbAlu"]
print "nbL1TD: ", totalCounts["nbL1TD"]
print "nbSVA: ", totalCounts["nbSVA"]
print "nbERVK: ", totalCounts["nbERVK"]
print "nbPSD: ", totalCounts["nbPSD"]
print "nbL1DEL: ", totalCounts["nbL1DEL"]
print "nbL1DUP: ", totalCounts["nbL1DUP"]
print 


## 5. Write dataframe into the output file
############################################
outPath = outDir + '/' + fileName + '.tsv'
rtCountsDf.to_csv(outPath, sep='\t') 

## End ##
print
print "***** Finished! *****"
print
