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
import formats
import time
from operator import itemgetter, attrgetter, methodcaller
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import seaborn as sns

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('metadata', help='PCAWG donors metadata')
parser.add_argument('expr', help='Expression matrix')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
metadata = args.metadata
expr = args.expr
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "metadata: ", metadata
print "expr: ", expr
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ##

#### 1. Read metadata file and create dictionary with ...
info("1. Read metadata file and create dictionary with ...")

## Create dictionary with RNA-seq aliquot id and tumour histology equivalences
metadata = open(metadata, 'r')

histologyDict = {}

## For each donor:
for line in metadata:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")

        tumorHistology = line[12]
        tumor_wgs_has_matched_rna_seq = line[24]
        tumor_rna_seq_aliquot_ids = line[30]

        if (tumor_wgs_has_matched_rna_seq == "True"):
            aliquotIdList = tumor_rna_seq_aliquot_ids.split(",")
            
            for aliquotId in aliquotIdList:

                histologyDict[aliquotId] = tumorHistology 


## 1078 tumor aliquots with RNA-seq available
histologySeries =  pd.Series(histologyDict)
print "histologySeries: ", histologySeries.shape

tumorAliquotIdList =  list(histologySeries.index)

#### 2. 
info("2. ")

exprDf = pd.read_csv(expr, header=0, index_col=0, sep='\t')

exprDf = exprDf.T

## 1521 total aliquots (tumor+normal) with RNA-seq available.
# I guess ~1000 are tumor and 500 are normal
print "exprDf: ", exprDf.shape

## Select tumor aliquouts with tumor histology available
exprFilteredDf = exprDf.loc[tumorAliquotIdList, :]

## Add the tumor histology information to the histology dataframe
exprFilteredDf["tumorHistology"] = histologySeries

print "exprFilteredDf: ", exprFilteredDf

## Select expression values for Lung-SCC aliquots:
exprFilteredLungSCCDf = exprFilteredDf[exprFilteredDf["tumorHistology"] == "Lung-SCC"]

print exprFilteredLungSCCDf.shape

exprFilteredLungTargetSampleDf = exprFilteredLungSCCDf.loc["78ed13d4-1578-46b8-81fa-b4395643691d"]

exprTargetSample = exprFilteredLungTargetSampleDf.loc["CCND1"]

print "targetSample: ", exprTargetSample

exprAllLungSamples = exprFilteredLungSCCDf["CCND1"].values

print "exprAllLungSamples: ", exprAllLungSamples


## Make statistical test
###########################

## Make boxplot:
#################
sns.set_style("whitegrid")

#title = "U: " + str(Ustat) + ", p-value: " + str(pvalue)

fig = plt.figure(figsize=(4,5))
#fig.suptitle(title, fontsize=14)

ax = sns.boxplot(y="CCND1",data=exprFilteredLungSCCDf, width=0.5, showfliers=False)
ax = sns.stripplot(y="CCND1", data=exprFilteredLungSCCDf, jitter=True, color=".3")
ax.plot(exprTargetSample, 'or', markeredgecolor='black', markeredgewidth='1', markersize='5')

ax.set_ylabel('CCND1 (FPKM)')
ax.set_xlabel('')


# Add mann whitney U statistic and p-value to the plot:

## Save figure 
fileName = outDir + "/CCND1_expr_lungSCC.pdf"
plt.savefig(fileName)


## End ##
print
print "***** Finished! *****"
print
