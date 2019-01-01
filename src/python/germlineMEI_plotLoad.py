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

def pointSize(alleleFreq):
    """
    """
    pointSizes = [20*1.25**n for n in range(0,11)]
    index = int(alleleFreq / 10)
    pointSize = pointSizes[index]

    return pointSize

def classify(cluster, percSamples):
    """
    """
    if cluster !=  -1:
        activityType = 'no-hot' 

    elif percSamples >= 2:
        activityType = 'strombolian'
    
    else:
        activityType = 'plinian'

    return activityType


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats
import formats

## Graphic style ##
sns.set_style("white")
sns.set_style("ticks")

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('rtCounts', help='')
parser.add_argument('metrics', help='')
parser.add_argument('vcf', help='')
parser.add_argument('metadata', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
rtCounts = args.rtCounts
metrics = args.metrics
vcf = args.vcf
metadata = args.metadata

outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "rtCounts: ", rtCounts
print "metrics: ", metrics
print "vcf: ", vcf
print "metadata: ", metadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ## 

#### Read input 
################
rtCountsDf = pd.read_csv(rtCounts, header=0, index_col=0, sep='\t')
metricsDf = pd.read_csv(metrics, header=0, index_col=0, sep='\t')

## Read VCF file
VCFObj = formats.VCF()
donorIdList = VCFObj.read_VCF_multiSample(vcf)

## Make dict with donor id equivalencies
metadata = open(metadata, 'r')
donorIdEq = {}

## For sample's metadata:
for line in metadata:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
       
        submitterDonorId = line[0]
        donorId = line[1]
        histology = line[12]
        donorIdEq[submitterDonorId] = (donorId, histology)

#### Initialize donor load dictionary
######################################
donorLoad = {}

for donorId in donorIdList:
    donorId = donorIdEq[donorId][0]
    donorLoad[donorId] = {}
    donorLoad[donorId]['src'] = 0
    donorLoad[donorId]['hot'] = 0
    donorLoad[donorId]['plinian'] = 0
    donorLoad[donorId]['strombolian'] = 0
    donorLoad[donorId]['histology'] = ''

print 'donorLoad: ', donorLoad

#### Fill donor load dictionary
#################################
for MEIObj in VCFObj.lineList:

    SRCID = MEIObj.infoDict['SRCID']

    ## Select only those MEI that passes all the filters and that are active in PCAWG cohort
    if (MEIObj.filter == "PASS") and (SRCID in metricsDf.index):

        activityType = metricsDf.loc[SRCID]['activityType']

        # A) MEI absent in reference genome
        if (MEIObj.alt == "<MEI>"):

            # For each donor and genotype
            for donorId, genotypeField in MEIObj.genotypesDict.iteritems():
                donorId, histology = donorIdEq[donorId]
                genotypeFieldList = genotypeField.split(":")
                genotype = genotypeFieldList[0]   

                if donorLoad[donorId]['histology'] == '':
                    donorLoad[donorId]['histology'] = histology             

                # a) Carrier
                if (genotype == '1/1') or (genotype == '0/1') or (genotype == '1'):
                    donorLoad[donorId]['src'] += 1

                    if activityType == 'strombolian':
                        donorLoad[donorId]['hot'] += 1
                        donorLoad[donorId]['strombolian'] += 1

                    elif activityType == 'plinian':
                        donorLoad[donorId]['hot'] += 1
                        donorLoad[donorId]['plinian'] += 1

                # c) Not carrier

        ## B) MEI in the reference genome 
        elif (MEIObj.ref == "<MEI>"):

            # For each donor and genotype
            for donorId, genotypeField in MEIObj.genotypesDict.iteritems():
                donorId, histology = donorIdEq[donorId]
                genotypeFieldList = genotypeField.split(":")
                genotype = genotypeFieldList[0]

                if donorLoad[donorId]['histology'] == '':
                    donorLoad[donorId]['histology'] = histology  

                # a) Carrier
                if (genotype == '0/0') or (genotype == '0/1') or (genotype == '0'):

                    donorLoad[donorId]['src'] += 1

                    if activityType == 'strombolian':
                        donorLoad[donorId]['hot'] += 1
                        donorLoad[donorId]['strombolian'] += 1

                    elif activityType == 'plinian':
                        donorLoad[donorId]['hot'] += 1
                        donorLoad[donorId]['plinian'] += 1

                # c) Not carrier

    print '-------------------------'


print 'donorLoad: ', donorLoad


#### Convert dictionary into dataframe and incorporate RT counts
#######################################

## Create list of tuples
tuples = []
for donorId in donorLoad:
    nbL1 = rtCountsDf.loc[donorId]['nbL1']
    nbSrc = donorLoad[donorId]['src']
    nbHot = donorLoad[donorId]['hot']
    nbStrombolian = donorLoad[donorId]['strombolian']
    nbPlinian = donorLoad[donorId]['plinian']
    histology = donorLoad[donorId]['histology']


    tuples.append((donorId, histology, nbL1, nbSrc, nbHot, nbStrombolian, nbPlinian))

## Create dataframe
loadDf = pd.DataFrame(tuples, columns=['donorId', 'histology', 'nbL1', 'nbSrc', 'nbHot', 'nbStrombolian', 'nbPlinian'])
loadDf.set_index('donorId', inplace=True)

## Save metrics into file
outFilePath = outDir + '/source_elements_load.tsv'
loadDf.to_csv(outFilePath, sep='\t') 

print 'loadDf: ', loadDf

####
header("Finished")

