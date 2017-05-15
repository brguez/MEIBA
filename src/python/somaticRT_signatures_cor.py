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

def scatterCorr(listA, listB, threshold, fileName):
    """
        listA - number of somatic L1 insertions per individual
        listB - Exposure to a given signature

        Interpretation of strength of correlation

        very weak: < 0,15 
        weak: 0,15-0,25  
        moderate: 0,25-0,40 
        strong: 0,40-0,75
        very strong: >0,75

    """
    corr = stats.spearmanr(listA, listB)
    coefficient = float(format(corr[0], '.3f'))
    pvalue = float(corr[1])
    

    ## Make scatterplot if rho >= threshold or <= -theshold
    if (coefficient >= threshold) or (coefficient <= -threshold):
                
        text = 'spearman_corr: ' + str(coefficient) + "; p_value: " + str(pvalue)  
    
        # Make scatterplot
        fig = plt.figure(figsize=(6,6))

        ax1 = fig.add_subplot(1, 1, 1)
        plt.scatter(listA, listB, color='#008000', alpha=.4)
        plt.xlim((-5, (max(listA) + 10)))
        plt.ylim((-5, (max(listB) + 10)))
        plt.xlabel('# L1', fontsize=12)
        plt.ylabel('Signature exposure', fontsize=12)
        ax1.text(0, 1.02, text, transform = ax1.transAxes)

        ## Add trend line:
        plt.plot(np.unique(listA), np.poly1d(np.polyfit(listA, listB, 1))(np.unique(listA)), color='#000000',)
    
        ## Save figure
        fileName = fileName + '_' + str(coefficient) + '_correlation.pdf' 

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

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('rtCounts', help='tsv with MEI allele counts')
parser.add_argument('signatures', help='tsv with MEI allele freqs')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
rtCounts = args.rtCounts
signatures = args.signatures
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "rtCounts: ", rtCounts
print "signatures: ", signatures
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### 1. Load input tables:
##########################
rtCountsDf = pd.read_csv(rtCounts, header=0, index_col=0, sep='\t')
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

