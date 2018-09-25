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

def classify(nbMEI):
    """
    """
    ## a) None (0 insertions):
    if (nbMEI == 0):
        category = "None"

    ## b) Low ([1-10] insertions):
    elif (nbMEI > 0) and (nbMEI <= 10):
        category = "Low"

    ## c) Moderate ((10-100] insertions):
    elif (nbMEI > 10) and (nbMEI <= 100):
        category = "Moderate"        

    ## d) High (>100 insertions):
    else:
        category = "High"        

    return category

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
parser.add_argument('rtCounts', help='Tabular text file containing one row per sample with retrotransposition counts')
parser.add_argument('histologyOrder', help='File containing histology ordering. One row per histology')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
rtCounts = args.rtCounts
histologyOrder = args.histologyOrder
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "rtCounts: ", rtCounts
print "histologyOrder: ", histologyOrder
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print


## Start ##

### 0. Read histology and create list with histology ordering
histologyFile = open(histologyOrder, 'r')
histologyList = []

for line in histologyFile:
    line = line.rstrip('\n')
    line = line.split("\t")
    histology = line[0]
    histologyList.append(histology)

## Move Myeloid-MPN to the end as no z-score nor p-value due to no insertion found
histologyList.append(histologyList.pop(histologyList.index('Myeloid-MPN')))

### 1. Read counts file
rtCountsDf = pd.read_csv(rtCounts, header=0, index_col=0, sep='\t')

### 2. Classify each sample in one of these categories:
# None: 0 insertions
# Low: 1-10 insertions
# Moderate: 10-100 insertions
# High: > 100 insertions
rtCountsDf["category"] = rtCountsDf["nbTotal"].apply(classify)

# Save dataframe into tsv
outFilePath = outDir + '/supplementaryTable2.selectedHistologies.categories.df'
rtCountsDf.to_csv(outFilePath, sep='\t') 

### 3. Make dataframe with the number of samples per category for each tumor type
# Project codes: columns
# Categories number of samples:
#           tumorType1   tumorType2   tumorType3....
# None      X1           Y1           Z1     
# Low       X2           Y2           Z2
# ...
catCountsDf = rtCountsDf.groupby(["tumorType", "category"]).size().reset_index(name="nbSamples")

## Rearrange dataframe
catCountsDf = catCountsDf.set_index('category')

catCountsDf = catCountsDf.pivot_table(values='nbSamples', index=catCountsDf.index, columns='tumorType', aggfunc='first')

## Replace NaNs by 0s
catCountsDf = catCountsDf.fillna(0)

### 4. Make dataframe with the percentage of samples per category for each tumor type
# Project codes: columns
# Categories % samples:
#                   tumorType1 tumorType2 tumorType3....
# 0_insertions      X1%          Y1%          Z1%     
# 1-10_insertions   X2%          Y2%          Z2%

nbSamplesPerTumorTypeSerie = catCountsDf.sum(axis=0)
categories = catCountsDf.index
tumorTypes = catCountsDf.columns

catPercDf = pd.DataFrame(index=categories, columns=tumorTypes)

# Iterate over row index labels (activity categories)
for category in categories:

    # Iterate over column index labels (tumorTypes)
    for tumorType in tumorTypes:
    
        catCountTumorType = catCountsDf.loc[category, tumorType]
        nbSamplesTumorType = nbSamplesPerTumorTypeSerie.loc[tumorType]

        # Compute the percentage
        catPercTumorType = float(catCountTumorType)/float(nbSamplesTumorType) * 100

        # Add source element contribution to dataframe
        catPercDf.loc[category, tumorType] = catPercTumorType

### 5. Rearrange dataframe
histologyList = list(reversed(histologyList))

transposedDf = catPercDf.transpose()
sortedDf = transposedDf.loc[histologyList]
catPercSortedDf = sortedDf.transpose()

### 6. Make list per category containing the percentage of samples belonging to a given category in each tumor type
# list 0 insertions [%tumorType1, %tumorType2, ... ]
# list 1-10 insertions [%tumorType1, %tumorType2, ... ]
# ...
percHighList, percLowList, percModerateList, percNoneList = catPercSortedDf.values.tolist()

###  7. Make bar plot
ypos = np.arange(1, len(percHighList) + 1)    # the y locations for the groups
height = 0.75      # the width of the bars: can also be len(x) sequence

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
plt.yticks(ypos, histologyList, fontsize=8)

# Rotate them
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)

## Add the number of samples per tumor type on the top of each bar
sampleCountsSortedSeries = nbSamplesPerTumorTypeSerie.reindex(histologyList)
sampleCountsList = sampleCountsSortedSeries.values.tolist()
autolabel(p4, ax, sampleCountsList) ## autolabel function

## Make legend
circle1 = mpatches.Circle((0, 0), 5, color='#DCDCDC', alpha=0.90)
circle2 = mpatches.Circle((0, 0), 5, color='#fff68f', alpha=0.90)
circle3 = mpatches.Circle((0, 0), 5, color='#fdb913', alpha=0.90)
circle4 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.90)

l = plt.figlegend((circle1, circle2, circle3, circle4), ('0', '[1-10]', '(10-100]', '>100'), loc='upper center', ncol=4, labelspacing=0.75, title="# MEI", fontsize=10, fancybox=True)

## Save figure
fileName = outDir + "/Pictures/rtRates_histologies.pdf"
plt.savefig(fileName)

fileName = outDir + "/Pictures/rtRates_histologies.svg"
plt.savefig(fileName)

## End ##
print
print "***** Finished! *****"
print
