#!/usr/bin/env python
#coding: utf-8

### Functions ###
def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, "****", string, "****"

def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string

def overlap(begA, endA, begB, endB):

    """
    Check if both ranges overlap. 2 criteria for defining overlap: 
    ## A) Begin of the range A within the range B         
    #       *beg* <---------range_A---------->                         
    # <---------range_B----------> 
                
    #    *beg* <-------range_A----->
    # <-------------range_B------------------>

    ## B) Begin of the range B within the range A     
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->
          
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>

    """    
       
    # a) Begin of the range A within the range B   
    if ((begA >= begB) and (begA <= endB)):
        overlap = True
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True

    # c) Ranges do not overlapping
    else:
        overlap = False

    return overlap


def computeMetrics(VCFObj, MEIDbObj):
    """
    Precision. How many of the positively classified were relevant: 
        
        TP / P = TP / (TP + FP)
    
    Recall. How good a test is at detecting the positives: 

        TP / T = TP / (TP + FN)

    """   

    ### 1. Initialize counters
    insertionLenTupleList = []
    countsDict = {}
    strandDict = {}
    strandDict["consistent"] = 0
    strandDict["inconsistent"] = 0

    for Type in ["TD0", "TD1", "TD2"]:
    
        countsDict[Type] = {}
        countsDict[Type]["TP"] = 0
        countsDict[Type]["FP"] = 0
        countsDict[Type]["TN"] = 0
        countsDict[Type]["FN"] = 0


    ### 2. Intersect each detected MEI with the simulated list and count the number of false positives
    # A false positive will be a call not included in the simulated list
    MEIObjList = VCFObj.lineList

    # For each called MEI
    for MEIObj in MEIObjList:

        ## Select only those MEIs that passes the filters
        if (MEIObj.filter == "PASS"):

            length = MEIObj.infoDict["LEN"] if "LEN" in MEIObj.infoDict else "NA"
            strand = MEIObj.infoDict["STRAND"] if "STRAND" in MEIObj.infoDict else "NA"
            Type = MEIObj.infoDict["TYPE"]

            #info("Checking if " + MEIObj.chrom + ":" + str(MEIObj.pos) + " MEI is a TP..." )
         
            falsePositive = True # Assume it's a false positive till match 

            ### Assess if it overlaps with a simulated MEI 
            if (MEIObj.chrom in MEIDbObj.MEIDict):

                # For each simulated MEI in the chromosome
                for simMEIpos in MEIDbObj.MEIDict[MEIObj.chrom]:

                    simMEIObj = MEIDbObj.MEIDict[MEIObj.chrom][simMEIpos]

                    ## Define identified insertion coordinates range
                    begA = int(MEIObj.pos) - 50
                    endA = int(MEIObj.pos) + 50

                    ## Define simulated insertion coordinates range                    
                    begB = int(simMEIpos) - 50
                    endB = int(simMEIpos) + 50
    
                    ## True positive as MEI was simulated          
                    if overlap(int(begA), int(endA), int(begB), int(endB)):

                        falsePositive = False
                        simMEIObj.detected = True

                        ## Solo insertions length
                        if (Type == "TD0") and (length != "NA") and (int(length) < 6500) and (simMEIObj.Type == "TD0"):
                            insertionLenTupleList.append((int(length), int(simMEIObj.length)))

                        ## Insertions strand
                        # Strand available
                        if strand != "NA":

                            # a) Consistent strand
                            if (strand == simMEIObj.strand):
                                strandDict["consistent"] += 1
    
                            # b) Inconsistent strand
                            else:
                                strandDict["inconsistent"] += 1                              

                        break

            # False positive as not included 
            if falsePositive:
                countsDict[Type]["FP"] += 1

    ### 3. Count the number of true positives and false negatives
    for chrom in MEIDbObj.MEIDict:
        for pos in MEIDbObj.MEIDict[chrom]:
            simMEIObj = MEIDbObj.MEIDict[chrom][pos]
        
            ## a) True positive
            if simMEIObj.detected:
                countsDict[simMEIObj.Type]["TP"] += 1

            ## b) False negative
            else:   
                countsDict[simMEIObj.Type]["FN"] += 1

    ### 4. Compute precision and recall 
    metricsDict = {}

    for Type in ["TD0", "TD1", "TD2"]:
    
        TP = float(countsDict[Type]["TP"])
        FP = float(countsDict[Type]["FP"])
        FN = float(countsDict[Type]["FN"])

        precision = TP / (TP + FP)
        recall = TP / (TP + FN)

        ## Add to the dictionary
        metricsDict[Type] = {}
        metricsDict[Type]["Precision"] = precision
        metricsDict[Type]["Recall"] = recall             
     
    return metricsDict, insertionLenTupleList, strandDict




#### CLASSES ####
class MEIDb():
    """
            .....................

            Methods:
            -

    """
    def __init__(self):
        """

        """
        self.MEIDict = {}

    def read_bed(self, inputBed):
        """

        """
        # Open germline MEI database input file
        inputBedFile = open(inputBed, 'r')

        ## Read first input file line by line
        for line in inputBedFile:
            line = line.rstrip('\n')

            ## Skip header
            if not line.startswith("#"):
                line = line.split('\t')

                ## Create germline MEI object: 
                chrom = line[0]
                beg = line[1]
                end = line[2]
                Type = line[3]
                length = line[4]
                strand = line[5]
                inGap = line[6]                

                MEIObj = MEI(chrom, beg, end, Type, length, strand)

                ### Create a nested dictionary containing for each chromosome a dictionary of
                ## germline MEI objects organized by coordinates
                # key1(chrId) -> value(dict2) -> key2(coordinate) -> value(MEIObj)

                ## Create nested dictionaries as needed
                # B) First MEI in a given chromosome
                if MEIObj.chrom not in self.MEIDict:
                    self.MEIDict[MEIObj.chrom] = {}

                ## Add germline MEI object to the dictionary
                self.MEIDict[MEIObj.chrom][MEIObj.beg] = MEIObj

        # Close germline MEI database file
        inputBedFile.close()

class MEI():
    """
            .....................

            Methods:
            -

    """
    def __init__(self, chrom, beg, end, Type, length, strand):
        """

        """
        self.chrom = chrom
        self.beg = beg
        self.end = end
        self.Type = Type
        self.length = length
        self.strand = strand
        self.detected = False


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import scipy.stats as stats
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from operator import itemgetter, attrgetter, methodcaller
from collections import Counter
import seaborn as sns



## Get user's input ##
parser = argparse.ArgumentParser(description= "Compare PCAWG with 1KGP MEI and obtain different statistics")
parser.add_argument('paths', help='Path to output VCFs')
parser.add_argument('simMEI', help='BED file containing the simulated MEI')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
paths = args.paths
simMEI = args.simMEI
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "paths: ", paths
print "simMEI: ", simMEI
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

###########
## Start ## 
###########


##################################################################
## 1. Create database object containing the simulated insertions #
##################################################################
header("1. Create database object containing the simulated insertions")
MEIDbObj = MEIDb()
MEIDbObj.read_bed(simMEI)

###################################################
## 2. Assess PCAWG RT pipeline for each clonality #
###################################################
header("2. Assess PCAWG RT pipeline for each clonality")

paths = open(paths, 'r')

metricsDict = {}
lengthsDict = {}
strandDict = {}

# Read file line by line
for line in paths:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        vaf = fieldsList[0]
        VCF = fieldsList[1]

        ## Read VCF
        VCFObj = formats.VCF()
        VCFObj.read_VCF(VCF)

        ## Compute metrics
        metricsDict[vaf], lengthsDict[vaf], strandDict[vaf] = computeMetrics(VCFObj, MEIDbObj)

print "metricsDict: ", metricsDict

print "lengthsDict: ", lengthsDict

print "strandDict: ", strandDict

# Convert metrics dictionary into dataframe.
metricsList = []

for vaf in metricsDict:
    for Type in metricsDict[vaf]:
        tmpDict = metricsDict[vaf][Type]    
        tmpDict["VAF"] = vaf
        tmpDict["MEI"] = Type
        metricsList.append(tmpDict)

metricsDf = pd.DataFrame(metricsList)


#######################################
## 3. Plot information gathered in 2) #
#######################################
sns.set_style("whitegrid")

## 3.1 Recall and precision
############################
fig = plt.figure(figsize=(4,4))

g = sns.lmplot(x="Recall", y="Precision", data=metricsDf, fit_reg=False, legend=False)
 
#Clear the axes containing the scatter plot
g.ax.cla()

# Define colors and markers
markersDict = {}
markersDict["100"] = "o"
markersDict["75"] = "^"
markersDict["50"] = "s"
markersDict["25"] = "p"

colorsDict = {}
colorsDict["TD0"] = "black"
colorsDict["TD1"] = "blue"
colorsDict["TD2"] = "green"

## Plot each individual point separately
for i,row in enumerate(metricsDf.values):
    
    Type, precision, recall, vaf = row
    g.ax.plot(recall, precision, color=colorsDict[Type], marker=markersDict[vaf])

## Axis limits
sns.plt.ylim(0.9, 1.05)
sns.plt.xlim(0.9, 1)
g.set_axis_labels('Recall', 'Precision')

## Save figure 
fileName = outDir + "/Pictures/precision_recall.pdf"
plt.savefig(fileName)


## 3.2 MEI length
##################

## Prepare data
tmpList = map(list, zip(*lengthsDict["50"]))
predictedLenList = tmpList[0]
simLenList = tmpList[1]

## Compute correlation
corr = stats.pearsonr(predictedLenList, simLenList)
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

print "CORRELATION: ", pvalue,  text

# Make scatterplot
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)
plt.scatter(predictedLenList, simLenList, color='#6497b1', alpha=.4, s=30)
plt.xlim((0, (max(predictedLenList) + 500)))
plt.ylim((0, (max(simLenList) + 500)))
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('Simulated', fontsize=12)
ax.text(0.5, 0.1, text, transform = ax.transAxes)


## Save figure 
fileName = outDir + "/Pictures/length_corr.svg"
plt.savefig(fileName)

fileName = outDir + "/Pictures/length_corr.pdf"
plt.savefig(fileName)


## 3.3 MEI strand
##################
## Prepare data
consistentList = []
inconsistentList = []

for vaf in ["100", "75", "50", "25"]:
    total = float(strandDict[vaf]['consistent'] + strandDict[vaf]['inconsistent'])
    percConsistent = strandDict[vaf]['consistent'] / total * 100
    percInconsistent = strandDict[vaf]['inconsistent'] / total * 100

    consistentList.append(percConsistent)
    inconsistentList.append(percInconsistent)

print "consistentList: ", consistentList
print "inconsistentList: ", inconsistentList

##  Make figure
# Make bar plot
xpos = np.arange(4)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(5,5))
plt.ylabel('%', fontsize=12)
ax = fig.add_subplot(111)
p1 = ax.bar(xpos, consistentList, color='#008000', alpha=0.75, edgecolor='#000000', width=width, align='center')
p2 = ax.bar(xpos, inconsistentList, color='#ff0000', alpha=0.75, edgecolor='#000000', width=width, align='center',
             bottom=consistentList)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize ticks
plt.yticks(np.arange(0, 101, 10))
plt.xticks(xpos, ("100", "75", "50", "25"))


## Make legend
circle1 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.75)
circle2 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.75)
plt.figlegend((circle1, circle2), ('Consistent', 'Inconsistent'), loc = 'lower center', ncol=2, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/Pictures/strand_consistency.pdf"
plt.savefig(fileName)

## End ##
print
print "***** Finished! *****"
print
