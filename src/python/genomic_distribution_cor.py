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

def scatterCorr(arrayA, arrayB, threshold, outPath):
    """

        Interpretation of strength of correlation

        very weak: < 0,15 
        weak: 0,15-0,25  
        moderate: 0,25-0,40 
        strong: 0,40-0,75
        very strong: >0,75

    """
    corr = stats.spearmanr(arrayA, arrayB)
    coefficient = float(format(corr[0], '.3f'))
    pvalue = float(corr[1])
    print "pvalue: ", pvalue

    ## Make scatterplot if rho >= threshold or <= -theshold
    if (coefficient >= threshold) or (coefficient <= -threshold):
    
        # Make scatterplot
        fig = plt.figure(figsize=(6,6))
        ax1 = fig.add_subplot(1, 1, 1)
        #plot = sns.jointplot(x=arrayA, y=arrayB, kind="hex", xlim=(0,40), gridsize=50, dropna=True, cmap="Blues", stat_func=spearmanr)
        plot = sns.jointplot(x=arrayA, y=arrayB, kind="kde", space=0, xlim=(0,30), gridsize=50, dropna=True, cmap="Blues", stat_func=spearmanr)
        plt.xlabel('# L1', fontsize=12)
        plt.ylabel('Replication time', fontsize=12)

#        sns.plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  # shrink fig so cbar is visible
#        cax = plot.fig.add_axes([.85, .25, .05, .4])  # x, y, width, height
#        sns.plt.colorbar(cax=cax)

        #sns.jointplot(x=arrayA, y=arrayB, kind="kde", space=0, color="b", xlim=(0,30))



        ## Save figure
        fileName = outPath + '_' + str(coefficient) + '_correlation.pdf' 

        plt.savefig(fileName)

    return coefficient, pvalue


#### MAIN ####

## Import modules ##
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
from scipy.stats import spearmanr

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('input', help='tsv...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
inputFile = args.input
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

#### 1. Load input tables:
##########################
inputDf = pd.read_csv(inputFile, header=0, sep='\t')

print "HOLA: ", inputDf

fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("nbL1", "medianRT", data=inputDf, xlim=(0,30), kind="kde", space=0, dropna=True, cmap="Blues", stat_func=spearmanr)

## Save figure
outPath = outDir + '/nbL1_medianRT_correlation.pdf'
plt.savefig(outPath)

fig = plt.figure(figsize=(6,6))
ax2 = fig.add_subplot(1, 1, 1)
plot = sns.jointplot("nbL1", "medianExpr", data=inputDf, xlim=(0,30), ylim=(0,700000), shade_lowest=True, kind="kde", space=0, dropna=True, cmap="Blues", stat_func=spearmanr)

## Save figure
outPath = outDir + '/nbL1_medianExpression_correlation.pdf'
plt.savefig(outPath)

sys.exit(1)
signaturesDf = pd.read_csv(signatures, header=0, index_col=0, sep='\t')
signaturesList = signaturesDf.columns.values.tolist()

## Add the tumor type and the number of somatic L1 insertions to the signatures dataframe
signaturesDf.insert(0, 'nbL1', rtCountsDf['nbL1'])
signaturesDf.insert(0, 'tumorType', rtCountsDf['dcc_project_code'])

## Remove those samples excluded from RT analysis (14 specimens are removed)
signaturesFilteredDf = signaturesDf.dropna()


#### 2. Compute nbL1/signatures correlation for each tumor type
################################################################
header("2. Compute nbL1/signatures correlation for each tumor type")

## Make list with tumor types and signatures:

tumorTypesList = set(signaturesFilteredDf['tumorType'].tolist())

coeffDict = {}
pvalueDict = {}

## For each tumor type
for tumorType in tumorTypesList:
    
    tumorTypeSignaturesDf = signaturesFilteredDf[signaturesFilteredDf['tumorType']  == tumorType]
    
    coeffDict[tumorType] = {}
    pvalueDict[tumorType] = {}

    ## For each signature
    for signature in signaturesList:
        
        nbL1 = tumorTypeSignaturesDf['nbL1'].tolist()
        exposures = tumorTypeSignaturesDf[signature].tolist()
       
        # a) All values equal to 0 either for L1 events or signatures exposure
        if all(v == 0 for v in nbL1) or all(v == 0 for v in exposures):
            coefficient = "NA"
            pvalue = "NA"
        
        # b) Not all values equal to 0. 
        else:
            fileName = outDir + '/Pictures/' + tumorType + '_' + signature 
            coefficient, pvalue = scatterCorr(nbL1, exposures, 0.25, fileName)

        ## Save rho and pvalues into dictionary
        coeffDict[tumorType][signature] = coefficient
        pvalueDict[tumorType][signature] = pvalue



## 3. Convert into dataframes and print output
#### Coefficients
# Create pandas dataframe from dictionary
coeffDf = pd.DataFrame(coeffDict) 

# transpose dataframe
coeffDf = coeffDf.T 

# Reorder signatures
coeffDf = coeffDf[signaturesList]

# Save output into tsv
outFilePath = outDir + '/RT_tumortypes_signatures_corr.coefficients.tsv'
coeffDf.to_csv(outFilePath, sep='\t') 

#### P-values
# Create pandas dataframe from dictionary
pvalueDf = pd.DataFrame(pvalueDict) 

# transpose dataframe
pvalueDf = pvalueDf.T 

# Reorder signatures
pvalueDf = pvalueDf[signaturesList]

# Save output into tsv
outFilePath = outDir + '/RT_tumortypes_signatures_corr.pvalues.tsv'
pvalueDf.to_csv(outFilePath, sep='\t') 

####
header("Finished")

