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

####### CLASSES #######
class cohort():
    """
        .....................

        Methods:
        -
    """
    def __init__(self):
        """
        """
        self.VCFdict = {}


    def read_VCFs(self, inputPath):
        """
        """
        inputFile = open(inputPath, 'r')

        info("Read input VCFs ")

        # Per iteration, read a VCF, generate a VCF object and add it to the cohort
        for line in inputFile:
            line = line.rstrip('\n')
            line = line.split("\t")

            projectCode = line[0]
            sampleId = line[1]
            VCFfile = line[2]

            # Create VCF object
            VCFObj = formats.VCF()

            info("Reading " + VCFfile + "...")

            # Input VCF available
            if os.path.isfile(VCFfile):

                # Read VCF and add information to VCF object
                VCFObj.read_VCF(VCFfile)

                # Add projectCode and sampleId information to the genotype field in each MEI object
                for MEIObject in VCFObj.lineList:

                    MEIObject.format = MEIObject.format + ':SAMPLEID'
                    MEIObject.genotype = MEIObject.genotype + ':' + sampleId

                # Initialize the donor list for a given project if needed 
                if projectCode not in self.VCFdict:

                    self.VCFdict[projectCode] = []   

                # Add donor VCF to cohort
                self.VCFdict[projectCode].append(VCFObj)

            else:
                print "[ERROR] Input file does not exist"


####### FUNCTIONS #######
def autolabel(rects, ax, valuesList):
    # Get y-axis height to calculate label position from.
    (y_bottom, y_top) = ax.get_ylim()
    y_height = y_top - y_bottom

    index = 0
    for rect in rects:
        value = valuesList[index]

        ax.text(rect.get_x() + rect.get_width()/2., 1.01*y_height,
                '%d' % int(value),
                ha='center', va='bottom', rotation=90)    


        index += 1

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
from operator import itemgetter, attrgetter, methodcaller
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: projectCode sampleId vcf_path')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ##

### 1. Initialize cohort object
cohortObj = cohort()

### 2. Read VCF files, create VCF objects and organize them
cohortObj.read_VCFs(inputPath)

print "test: ", cohortObj.VCFdict

### 3. Make a dictionary containing per tumor type the number of samples with the following retrotransposition activity categories:
# None: 0 insertions
# Low: 1-10 insertions
# Moderate: 10-100 insertions
# High: 100-500 insertions
# Very-high: > 500 insertions


categoryCountsDict = {}

for projectCode in cohortObj.VCFdict:

    ## Initialize category counts    
    categoryCountsDict[projectCode] = {}
    categoryCountsDict[projectCode]["None"] = 0
    categoryCountsDict[projectCode]["Low"] = 0
    categoryCountsDict[projectCode]["Moderate"] = 0
    categoryCountsDict[projectCode]["High"] = 0
    categoryCountsDict[projectCode]["Very-high"] = 0


    ## Count total number of donors per tumor type and number of donor per category
    for VCFObj in cohortObj.VCFdict[projectCode]:

        nbInsertions = int(len(VCFObj.lineList))
        print "nbInsertions: ", nbInsertions

        ## a) None:
        if (nbInsertions == 0):
            print "None"
            categoryCountsDict[projectCode]["None"] += 1
        ## b) Low:
        elif (nbInsertions > 0) and (nbInsertions <= 10):
            print "Low"
            categoryCountsDict[projectCode]["Low"] += 1

        ## c) Moderate
        elif (nbInsertions > 10) and (nbInsertions <= 100):
            print "Moderate"
            categoryCountsDict[projectCode]["Moderate"] += 1

        ## d) High
        elif (nbInsertions > 100) and (nbInsertions <= 500):
            print "High"
            categoryCountsDict[projectCode]["High"] += 1                

        ## e) Very-high  
        elif (nbInsertions > 500):
            print "Very-high"
            categoryCountsDict[projectCode]["Very-high"] += 1
        ## d) Error
        else:
            print "ERROR..."

### 4. Make dataframe with the number of samples per category for each tumor type
# Project codes: columns
# Categories number of samples:
#                   ProjectCode1 ProjectCode2 ProjectCode3....
# 0_insertions      X1           Y1           Z1     
# 1-10_insertions   X2           Y2           Z2
# ...

categoryCountsDataframe =  pd.DataFrame(categoryCountsDict)


### 5. Make dataframe with the percentage of samples per category for each tumor type
# Project codes: columns
# Categories % samples:
#                   ProjectCode1 ProjectCode2 ProjectCode3....
# 0_insertions      X1%          Y1%          Z1%     
# 1-10_insertions   X2%          Y2%          Z2%

donorsPerTumorTypeSerie = categoryCountsDataframe.sum(axis=0)

categories = categoryCountsDataframe.index
projecCodes = categoryCountsDataframe.columns

categoryPercDataframe = pd.DataFrame(index=categories, columns=projecCodes)

# Iterate over row index labels (activity categories)
for category in categories:

    # Iterate over column index labels (project codes)
    for projectCode in projecCodes:
    
        categoryCountProjectCode = categoryCountsDataframe.loc[category, projectCode]
        nbDonorsInTumortype = donorsPerTumorTypeSerie.loc[projectCode]

        # Compute the percentage
        donorPercTumorType = float(categoryCountProjectCode)/float(nbDonorsInTumortype) * 100

        # Add source element contribution to dataframe
        categoryPercDataframe.loc[category, projectCode] = donorPercTumorType

# Order dataframe columns (tumor types) first according to "High" and then to "None" category
transposedDf = categoryPercDataframe.transpose()
sortedDf = transposedDf.sort_values(['High', 'None'], ascending=[False, True])
categoryPercSortedDataframe = sortedDf.transpose()


### 6. Make list per category containing the percentage of samples belonging to a given category in each tumor type
# list 0 insertions [%ProjectCode1, %ProjectCode2, ... ]
# list 1-10 insertions [%ProjectCode1, %ProjectCode2, ... ]
# ...

percHighList, percLowList, percModerateList, percNoneList, percVeryHigh = categoryPercSortedDataframe.values.tolist()

### 7. Make ordered list containing the number of donors per tumor type

tumorTypeList = categoryPercSortedDataframe.columns.values.tolist()

donorCountsSortedSeries = donorsPerTumorTypeSerie.reindex(tumorTypeList)


###  8. Make bar plot
xpos = np.arange(2, len(percHighList) + 2)    # the x locations for the groups

width = 0.75      # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(14,8))
fig.suptitle('Number of samples', fontsize=12)
plt.ylabel('% Samples', fontsize=12)
ax = fig.add_subplot(111)

p1 = ax.bar(xpos, percVeryHigh, color='#ff0000', alpha=1, edgecolor='#000000', width=width, align='center')

p2 = ax.bar(xpos, percHighList, color='#f15a22', alpha=0.80, edgecolor='#000000', width=width, align='center',
             bottom=[i for i in percVeryHigh])

p3 = ax.bar(xpos, percModerateList, color='#fdb913', alpha=0.90, edgecolor='#000000', width=width, align='center',
             bottom=[i+j for i,j in zip(percVeryHigh, percHighList)])

p4 = ax.bar(xpos, percLowList, color='#fff68f', alpha=0.90, edgecolor='#000000', width=width, align='center',
             bottom=[i+j+x for i,j,x in zip(percVeryHigh, percHighList, percModerateList)])

p5 = ax.bar(xpos, percNoneList, color='#DCDCDC', alpha=0.80, edgecolor='#000000', width=width, align='center',
             bottom=[i+j+x+z for i,j,x,z in zip(percVeryHigh, percHighList, percModerateList, percLowList)])

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize axis
plot_xmargin = 2
plot_ymargin = 20

x0, x1, y0, y1 = plt.axis()
plt.axis((x0,
          x1 - plot_xmargin,
          y0,
          y1 - plot_ymargin))

## Customize ticks
plt.yticks(np.arange(0, 100.001, 10))
tumorTypesList = categoryPercSortedDataframe.columns.values.tolist()
plt.xticks(xpos, tumorTypesList)

# Rotate them
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)

## Add the number of samples per tumor type on the top of each bar
donorCountsList = donorCountsSortedSeries.values.tolist()
autolabel(p5, ax, donorCountsList) ## autolabel function

## Make legend
circle1 = mpatches.Circle((0, 0), 5, color='#DCDCDC', alpha=0.90)
circle2 = mpatches.Circle((0, 0), 5, color='#fff68f', alpha=0.90)
circle3 = mpatches.Circle((0, 0), 5, color='#fdb913', alpha=0.90)
circle4 = mpatches.Circle((0, 0), 5, color='#f15a22', alpha=0.90)
circle5 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.90)
a
l = plt.figlegend((circle1, circle2, circle3, circle4, circle5), ('0', '[1-10]', '(10-100]', '(100-500]', '>500'), loc = 'center', ncol=2, labelspacing=0.75, title="Number of retrotransposition events", fontsize=11, fancybox=True, bbox_to_anchor=(0.75, 0.75))

## Save figure
fileName = outDir + "/PCAWG_retrotranspositionRates_barPlot2.pdf"

plt.savefig(fileName)

## End ##
print
print "***** Finished! *****"
print
