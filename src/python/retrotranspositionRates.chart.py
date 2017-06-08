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

            donorId = line[0]
            projectCode = line[1]
            VCFfile = line[2]

            #print "tiooo: ", donorId, projectCode, VCFfile

            # Create VCF object
            VCFObj = formats.VCF()

            info("Reading " + VCFfile + "...")

            # Input VCF available
            if os.path.isfile(VCFfile):

                # Read VCF and add information to VCF object
                VCFObj.read_VCF(VCFfile)

                # Initialize the donor list for a given project if needed 
                if projectCode not in self.VCFdict:

                    self.VCFdict[projectCode] = []   

                # Add donor VCF to cohort
                self.VCFdict[projectCode].append(VCFObj)

            else:
                print "[ERROR] Input file does not exist"


####### FUNCTIONS #######
def autolabel(rects, ax, valuesList):
    # Get x-axis height to calculate label position from.
    (x_left, x_right) = ax.get_xlim()    
    x_length = x_right - x_left

    index = 0
    for rect in rects:
        value = valuesList[index]
        ax.text(1.05*x_length, rect.get_y(), 
                '%d' % int(value),
                ha='center', va='bottom', fontsize=8)    

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
parser.add_argument('inputPath', help='Tabular text file containing one row per donor with the following consecutive fields: projectCode donorId vcf_path')
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

### 3. Make a dictionary containing per tumor type the number of samples with the following retrotransposition activity categories:
# None: 0 insertions
# Low: 1-10 insertions
# Moderate: 10-100 insertions
# High: > 100 insertions

categoryCountsDict = {}

for projectCode in cohortObj.VCFdict:

    ## Initialize category counts    
    categoryCountsDict[projectCode] = {}
    categoryCountsDict[projectCode]["None"] = 0
    categoryCountsDict[projectCode]["Low"] = 0
    categoryCountsDict[projectCode]["Moderate"] = 0
    categoryCountsDict[projectCode]["High"] = 0

    ## Count total number of donors per tumor type and number of donor per category
    ## For each donor    
    for VCFObj in cohortObj.VCFdict[projectCode]:

        nbInsertions = 0

        # For each MEI in a given donor
        for MEIObj in VCFObj.lineList:
        
            # Only count MEI that passes the filters:
            if (MEIObj.filter == "PASS"):
                nbInsertions += 1
        
        ## a) None (0 insertions):
        if (nbInsertions == 0):
            categoryCountsDict[projectCode]["None"] += 1

        ## b) Low ([1-10] insertions):
        elif (nbInsertions > 0) and (nbInsertions <= 10):
            categoryCountsDict[projectCode]["Low"] += 1

        ## c) Moderate ((10-100] insertions):
        elif (nbInsertions > 10) and (nbInsertions <= 100):
            categoryCountsDict[projectCode]["Moderate"] += 1

        ## d) High (>100 insertions):
        elif (nbInsertions > 100):
            categoryCountsDict[projectCode]["High"] += 1                

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

print "categoryCountsDataframe: ", categoryCountsDataframe

## Filter out those tumor types with 0 donors with high retrotransposition activity. 
categoryCountsFinalDataframe = categoryCountsDataframe
#tumorTypesBoolSerie =  categoryCountsDataframe.loc['High'] >= 1 # Enable this for selecting a subset of tumor types
#categoryCountsFinalDataframe = categoryCountsDataframe[tumorTypesBoolSerie.index[tumorTypesBoolSerie]]

print "categoryCountsFinalDataframe: ", categoryCountsFinalDataframe

### Group those elements contributing less than 1% in the category 'Other'
tumorTypesBoolSerie =  categoryCountsDataframe.loc['High'] < 1

categoryCountsOtherDataframe = categoryCountsDataframe[tumorTypesBoolSerie.index[tumorTypesBoolSerie]]

print "categoryCountsOtherDataframe: ", categoryCountsOtherDataframe

categoryCountsSerieOther = categoryCountsOtherDataframe.sum(axis=1)

print "categoryCountsSerieOther: ", categoryCountsSerieOther

### Make final dataframe and add the other category as an additional row
colName = "Other" 

categoryCountsFinalDataframe[colName] = categoryCountsSerieOther
print "final-counts: ", categoryCountsFinalDataframe

### 5. Make dataframe with the percentage of samples per category for each tumor type
# Project codes: columns
# Categories % samples:
#                   ProjectCode1 ProjectCode2 ProjectCode3....
# 0_insertions      X1%          Y1%          Z1%     
# 1-10_insertions   X2%          Y2%          Z2%

donorsPerTumorTypeSerie = categoryCountsFinalDataframe.sum(axis=0)
#donorsPerTumorTypeSerie = categoryCountsDataframe.sum(axis=0)

categories = categoryCountsFinalDataframe.index
projecCodes = categoryCountsFinalDataframe.columns

categoryPercDataframe = pd.DataFrame(index=categories, columns=projecCodes)

# Iterate over row index labels (activity categories)
for category in categories:

    # Iterate over column index labels (project codes)
    for projectCode in projecCodes:
    
        categoryCountProjectCode = categoryCountsFinalDataframe.loc[category, projectCode]
        nbDonorsInTumortype = donorsPerTumorTypeSerie.loc[projectCode]

        # Compute the percentage
        donorPercTumorType = float(categoryCountProjectCode)/float(nbDonorsInTumortype) * 100

        # Add source element contribution to dataframe
        categoryPercDataframe.loc[category, projectCode] = donorPercTumorType

# Order dataframe columns (tumor types) first according to "High" and then to "None" category
transposedDf = categoryPercDataframe.transpose()
sortedDf = transposedDf.sort_values(['High', 'Moderate', 'Low', 'None'], ascending=[True, True, True, False])
categoryPercSortedDataframe = sortedDf.transpose()

print "categoryPercSortedDataframe: ", categoryPercSortedDataframe

### 6. Make list per category containing the percentage of samples belonging to a given category in each tumor type
# list 0 insertions [%ProjectCode1, %ProjectCode2, ... ]
# list 1-10 insertions [%ProjectCode1, %ProjectCode2, ... ]
# ...

percHighList, percLowList, percModerateList, percNoneList = categoryPercSortedDataframe.values.tolist()

print "percHighList,percLowList,percModerateList,percNoneList: ", percHighList, percLowList, percModerateList, percNoneList

### 7. Make ordered list containing the number of donors per tumor type
tumorTypeList = categoryPercSortedDataframe.columns.values.tolist()

print "tumorTypeList: ", tumorTypeList

donorCountsSortedSeries = donorsPerTumorTypeSerie.reindex(tumorTypeList)

print "donorCountsSortedSeries: ", donorCountsSortedSeries

###  8. Make bar plot
ypos = np.arange(1, len(percHighList) + 1)    # the y locations for the groups

height = 0.75      # the width of the bars: can also be len(x) sequence
#fig = plt.figure(figsize=(5,3))
fig = plt.figure(figsize=(5,7))
ax = fig.add_subplot(111)
ax.yaxis.set_label_position("right")
plt.ylabel('Number of samples', fontsize=10, labelpad=35)
plt.xlabel('% Samples', fontsize=10)

p1 = ax.barh(ypos, percHighList, color='#ff0000', alpha=0.90, edgecolor='#000000', height=height, align='center')

p2 = ax.barh(ypos, percModerateList, color='#fdb913', alpha=0.90, edgecolor='#000000', height=height, align='center',
             left=[i for i in percHighList])

p3 = ax.barh(ypos, percLowList, color='#fff68f', alpha=0.90, edgecolor='#000000', height=height, align='center',
             left=[i+j for i,j in zip(percHighList, percModerateList)])

p4 = ax.barh(ypos, percNoneList, color='#DCDCDC', alpha=0.90, edgecolor='#000000', height=height, align='center',
             left=[i+j+x for i,j,x in zip(percHighList, percModerateList, percLowList)])


# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=1)
ax.set_axisbelow(True)

## Customize axis
plot_xmargin = 20
plot_ymargin = 0

x0, x1, y0, y1 = plt.axis()
plt.axis((x0,
          x1 - plot_xmargin,
          y0,
          y1 - plot_ymargin))

## Customize ticks
plt.xticks(np.arange(0, 100.001, 10), fontsize=8)
tumorTypeList = categoryPercSortedDataframe.columns.values.tolist()
plt.yticks(ypos, tumorTypeList, fontsize=8)

# Rotate them
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)

## Add the number of samples per tumor type on the top of each bar
donorCountsList = donorCountsSortedSeries.values.tolist()
autolabel(p4, ax, donorCountsList) ## autolabel function

## Make legend
#circle1 = mpatches.Circle((0, 0), 5, color='#DCDCDC', alpha=0.90)
#circle2 = mpatches.Circle((0, 0), 5, color='#fff68f', alpha=0.90)
#circle3 = mpatches.Circle((0, 0), 5, color='#fdb913', alpha=0.90)
#circle4 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.90)

#l = plt.figlegend((circle1, circle2, circle3, circle4), ('0', '[1-10]', '(10-100]', '>100'), loc='upper center', ncol=4, labelspacing=0.75, title="Number of retrotransposition events", fontsize=10, fancybox=True)

## Save figure
fileName = outDir + "/PCAWG_retrotranspositionRates_barPlot.pdf"

plt.savefig(fileName)

## End ##
print
print "***** Finished! *****"
print
