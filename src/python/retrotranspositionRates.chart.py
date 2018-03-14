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
    def __init__(self, excludedDonorsList):
        """
        """
        self.excludedDonorsList = excludedDonorsList
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
            tumorTypeList = line[1].split(",")
            VCFfile = line[2]
            
            if (donorId not in self.excludedDonorsList) and ("UNK" not in tumorTypeList):

                # Create VCF object
                VCFObj = formats.VCF()

                info("Reading " + VCFfile + "...")

                # Input VCF available
                if os.path.isfile(VCFfile):

                    # Read VCF and add information to VCF object
                    VCFObj.read_VCF(VCFfile)

                    for tumorType in tumorTypeList:
                    
                        # Initialize the donor list for a given project if needed 
                        if tumorType not in self.VCFdict:

                            self.VCFdict[tumorType] = []   

                        # Add donor VCF to cohort
                        self.VCFdict[tumorType].append(VCFObj)
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
parser.add_argument('inputPath', help='Tabular text file containing one row per donor with the following consecutive fields: donorId tumorType vcf_path')
parser.add_argument('metadata', help='Tabular text file containing donor metadata info')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
metadata = args.metadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "metadata: ", metadata
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ##

### 0. Read metadata file and create list of donors to be excluded
donorMetadataFile = open(metadata, 'r')
excludedDonorsList = []

## For donor's metadata:
for line in donorMetadataFile:
    if not line.startswith("#"):

        line = line.rstrip('\n')
        line = line.split("\t")
        donorId = line[1]
        traficExclusionStatus = line[3]
        tumorTypeExclusionStatus = line[11]

        if (traficExclusionStatus == "Excluded") or (tumorTypeExclusionStatus == "Excluded"):
            excludedDonorsList.append(donorId)

# print len(excludedDonorsList) # 91

### 1. Initialize cohort object
cohortObj = cohort(excludedDonorsList)

### 2. Read VCF files, create VCF objects and organize them
cohortObj.read_VCFs(inputPath)

### 3. Make a dictionary containing per tumor type the number of samples with the following retrotransposition activity categories:
# None: 0 insertions
# Low: 1-10 insertions
# Moderate: 10-100 insertions
# High: > 100 insertions

categoryCountsDict = {}

for tumorType in cohortObj.VCFdict:

    ## Initialize category counts    
    categoryCountsDict[tumorType] = {}
    categoryCountsDict[tumorType]["None"] = 0
    categoryCountsDict[tumorType]["Low"] = 0
    categoryCountsDict[tumorType]["Moderate"] = 0
    categoryCountsDict[tumorType]["High"] = 0

    ## Count total number of donors per tumor type and number of donor per category
    ## For each donor    
    for VCFObj in cohortObj.VCFdict[tumorType]:

        nbInsertions = 0

        # For each MEI in a given donor
        for MEIObj in VCFObj.lineList:
        
            # Only count MEI that passes the filters:
            if (MEIObj.filter == "PASS"):
                nbInsertions += 1
        
        ## a) None (0 insertions):
        if (nbInsertions == 0):
            categoryCountsDict[tumorType]["None"] += 1

        ## b) Low ([1-10] insertions):
        elif (nbInsertions > 0) and (nbInsertions <= 10):
            categoryCountsDict[tumorType]["Low"] += 1

        ## c) Moderate ((10-100] insertions):
        elif (nbInsertions > 10) and (nbInsertions <= 100):
            categoryCountsDict[tumorType]["Moderate"] += 1

        ## d) High (>100 insertions):
        elif (nbInsertions > 100):
            categoryCountsDict[tumorType]["High"] += 1                

        ## d) Error
        else:
            print "ERROR..."

### 4. Make dataframe with the number of samples per category for each tumor type
# Project codes: columns
# Categories number of samples:
#                   tumorType1 tumorType2 tumorType3....
# 0_insertions      X1           Y1           Z1     
# 1-10_insertions   X2           Y2           Z2
# ...
categoryCountsDataframe =  pd.DataFrame(categoryCountsDict)

### 5. Make dataframe with the percentage of samples per category for each tumor type
# Project codes: columns
# Categories % samples:
#                   tumorType1 tumorType2 tumorType3....
# 0_insertions      X1%          Y1%          Z1%     
# 1-10_insertions   X2%          Y2%          Z2%

donorsPerTumorTypeSerie = categoryCountsDataframe.sum(axis=0)

categories = categoryCountsDataframe.index
projecCodes = categoryCountsDataframe.columns

categoryPercDataframe = pd.DataFrame(index=categories, columns=projecCodes)

# Iterate over row index labels (activity categories)
for category in categories:

    # Iterate over column index labels (project codes)
    for tumorType in projecCodes:
    
        categoryCountTumorType = categoryCountsDataframe.loc[category, tumorType]
        nbDonorsInTumortype = donorsPerTumorTypeSerie.loc[tumorType]

        # Compute the percentage
        donorPercTumorType = float(categoryCountTumorType)/float(nbDonorsInTumortype) * 100

        # Add source element contribution to dataframe
        categoryPercDataframe.loc[category, tumorType] = donorPercTumorType

# Order dataframe columns according to z-scores
tumorTypeList = ['Eso-AdenoCA', 'Head-SCC', 'Lung-SCC', 'ColoRect-AdenoCA', 'Stomach-AdenoCA', 'Uterus-AdenoCA', 'Ovary-AdenoCA', 'Cervix-SCC', 'Myeloid-MDS', 'Biliary-AdenoCA', 'Bladder-TCC', 'Bone-Epith', 'Lymph-BNHL', 'Panc-AdenoCA', 'Lung-AdenoCA', 'Lymph-CLL', 'Prost-AdenoCA', 'Bone-Leiomyo', 'CNS-PiloAstro', 'CNS-Oligo', 'Breast-AdenoCA', 'Myeloid-MPN', 'CNS-GBM', 'Thy-AdenoCA', 'Bone-Osteosarc', 'Panc-Endocrine', 'Kidney-ChRCC', 'Myeloid-AML', 'Skin-Melanoma', 'CNS-Medullo', 'Kidney-RCC', 'Liver-HCC']

## Exclude "Myeloid-MDS" as single sample:
tumorTypeList.remove('Myeloid-MDS')

tumorTypeList = list(reversed(tumorTypeList))

transposedDf = categoryPercDataframe.transpose()
sortedDf = transposedDf.loc[tumorTypeList]
categoryPercSortedDataframe = sortedDf.transpose()

### 6. Make list per category containing the percentage of samples belonging to a given category in each tumor type
# list 0 insertions [%tumorType1, %tumorType2, ... ]
# list 1-10 insertions [%tumorType1, %tumorType2, ... ]
# ...

percHighList, percLowList, percModerateList, percNoneList = categoryPercSortedDataframe.values.tolist()

### 7. Make ordered list containing the number of donors per tumor type
donorCountsSortedSeries = donorsPerTumorTypeSerie.reindex(tumorTypeList)

###  8. Make bar plot
ypos = np.arange(1, len(percHighList) + 1)    # the y locations for the groups

height = 0.75      # the width of the bars: can also be len(x) sequence
#fig = plt.figure(figsize=(5,3))
fig = plt.figure(figsize=(3,7))
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
