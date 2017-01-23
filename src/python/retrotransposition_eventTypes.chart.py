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
    # Get x-axis height to calculate label position from.
    (x_left, x_right) = ax.get_xlim()    
    x_length = x_right - x_left

    index = 0
    for rect in rects:
        value = valuesList[index]
        ax.text(1.04*x_length, rect.get_y(), 
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
parser.add_argument('inputPath', help='Tabular text file containing one row per sample with the following consecutive fields: projectCode sampleId vcf_path')
parser.add_argument('pseudoPath', help='[PROVISIONAL] Tabular text file containing one row per projectCode with the following consecutive fields: projectCode Nb.pseudogenes')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.')

args = parser.parse_args()
inputPath = args.inputPath
pseudoPath = args.pseudoPath
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "inputPath: ", inputPath
print "pseudoPath: ", pseudoPath
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ##

############# BAR CHART ###############
### 1. Initialize cohort object
cohortObj = cohort()

### 2. Read VCF files, create VCF objects and organize them
cohortObj.read_VCFs(inputPath)

### 3. Make a dictionary containing per tumor type the total number of events for each retrotrotransposition insertion type:
# * L1-solo
# * L1-td
# * Alu
# * SVA
# * ERVK
# * processed-pseudogene

eventCountsDict = {}

for projectCode in cohortObj.VCFdict:

    ## Initialize category counts    
    eventCountsDict[projectCode] = {}
    eventCountsDict[projectCode]["L1-solo"] = 0
    eventCountsDict[projectCode]["L1-transduction"] = 0
    eventCountsDict[projectCode]["Alu"] = 0
    eventCountsDict[projectCode]["SVA"] = 0
    eventCountsDict[projectCode]["ERVK"] = 0
    eventCountsDict[projectCode]["processed-pseudogene"] = 0

    ## Count total number of donors per tumor type and number of donor per category
    for VCFObj in cohortObj.VCFdict[projectCode]:

        for MEIObj in VCFObj.lineList:
            
            MEIClass = MEIObj.infoDict["CLASS"]
            MEIType = MEIObj.infoDict["TYPE"]

            ## a) L1-solo:
            if (MEIClass == "L1") and (MEIType == "TD0"):
                eventCountsDict[projectCode]["L1-solo"] += 1
                
            ## b) L1-transduction
            elif (MEIType == "TD1") or (MEIType == "TD2"):
                eventCountsDict[projectCode]["L1-transduction"] += 1

            ## c) Alu
            # Note: I added a provisional SCORE filtering
            # to ask for at least one breakpoint reconstructed 
            elif (MEIClass == "Alu") and (int(MEIObj.infoDict["SCORE"]) > 2):
                eventCountsDict[projectCode]["Alu"] += 1
               
            ## d) SVA
            # Note: I added a provisional SCORE filtering
            # to ask for at least one breakpoint reconstructed
            elif (MEIClass == "SVA") and (int(MEIObj.infoDict["SCORE"]) > 2):
                eventCountsDict[projectCode]["SVA"] += 1
                 
            ## e) ERVK
            # Note: I added a provisional SCORE filtering
            # to ask for at least one breakpoint reconstructed
            elif (MEIClass == "ERVK") and (int(MEIObj.infoDict["SCORE"]) > 2):
                eventCountsDict[projectCode]["ERVK"] += 1
                
            ## f) Processed pseudogene            
            elif (MEIClass == "PSD"):
                eventCountsDict[projectCode]["processed-pseudogene"] += 1
                            
            ## g) Unexpected value
            #else:
            #    print MEIObj.infoDict["CLASS"], "[ERROR] Unexpected MEI Class value"
    

print "eventCountsDict: ", eventCountsDict

### 4. Make dataframe with the number of events per type for each tumor type
# Project codes: columns
# Event type number of events:
#                   ProjectCode1 ProjectCode2 ProjectCode3....
# L1-solo      X1           Y1           Z1     
# L1-transduction   X2           Y2           Z2
# ...

eventCountsDataframe =  pd.DataFrame(eventCountsDict)

print "eventCountsDataframe: ", eventCountsDataframe


### PROVISIONAL -- ADD PSEUDOGENE INFORMATION
pseudoCounts = open(pseudoPath, 'r')

# Read file line by line
for line in pseudoCounts:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
        
        projectCode = str(fieldsList[0])
        nbPseudogenes = fieldsList[1]
   
        # Add nb. pseudogenes to the counts dataframe
        eventCountsDataframe.set_value("processed-pseudogene", projectCode, nbPseudogenes)

print "eventCountsDataframePseudo: ", eventCountsDataframe


### 5. Make dataframe with the percentage of events per type for each tumor type
# Project codes: columns
# Categories % samples:
#                   ProjectCode1 ProjectCode2 ProjectCode3....
# L1-solo           X1%          Y1%          Z1%     
# L1-transduction   X2%          Y2%          Z2%

nbEventsPerTumorTypeSerie = eventCountsDataframe.sum(axis=0)

print "nbEventsPerTumorTypeSerie:" , nbEventsPerTumorTypeSerie

eventTypes = eventCountsDataframe.index
projecCodes = eventCountsDataframe.columns

eventPercDataframe = pd.DataFrame(index=eventTypes, columns=projecCodes)

# Iterate over row index labels (activity categories)
for eventType in eventTypes:

    # Iterate over column index labels (project codes)
    for projectCode in projecCodes:
    
        eventCountProjectCode = eventCountsDataframe.loc[eventType, projectCode]
        nbEventsInTumortype = nbEventsPerTumorTypeSerie.loc[projectCode]

        # Compute the percentage
        eventPercTumorType = float(eventCountProjectCode)/float(nbEventsInTumortype) * 100

        # Add source element contribution to dataframe
        eventPercDataframe.set_value(eventType, projectCode, eventPercTumorType)
        
print "eventPercDataframe: ", eventPercDataframe


## Order dataframe columns (tumor types) in the same order as for the chart generated in "retrotranspositionRates.chart.py"
tumorTypeList = ['ESAD', 'HNSC', 'COAD', 'LUSC', 'STAD', 'UCEC', 'PRAD', 'MELA', 'BOCA', 'PACA', 'BRCA', 'LIRI', 'READ', 'CESC', 'OV', 'SARC', 'LIHC', 'GBM', 'THCA', 'BLCA', 'GACA', 'PAEN', 'KICH', 'BTCA', 'ORCA', 'SKCM', 'LINC', 'KIRP', 'LGG', 'LUAD', 'KIRC', 'DLBC', 'EOPC', 'LAML', 'RECA', 'CMDI', 'LICA', 'MALY', 'PBCA', 'CLLE']

# I need to reverse it to have the bars correctly placed...
tumorTypeList.reverse()

eventPercSortedDataframe = eventPercDataframe.reindex_axis(tumorTypeList, axis=1)

print "eventPercSortedDataframe: ", eventPercSortedDataframe

### 6. Make list per event type containing the percentage of events in each tumor category
# list 0 insertions [%ProjectCode1, %ProjectCode2, ... ]
# list 1-10 insertions [%ProjectCode1, %ProjectCode2, ... ]
# ...

AluList, ERVKList, L1SoloList, L1TDList, SVAList, PSDList = eventPercSortedDataframe.values.tolist()

### 7. Make ordered list containing the total number of insertions per tumor type
nbEventsPerTumorTypeSortedSerie = nbEventsPerTumorTypeSerie.reindex(tumorTypeList)

###  8. Make  bar plot
# Note: I will not represent ERVK as we only have one insertion in the full PCAWG cohort...

ypos = np.arange(1, len(AluList) + 1)    # the y locations for the groups

height = 0.75      # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(7, 12))
# fig.suptitle('Number of samples', fontsize=12)
ax = fig.add_subplot(111)
ax.yaxis.set_label_position("right")
plt.ylabel('Total number of MEI', fontsize=10, labelpad=40)
plt.xlabel('% MEI', fontsize=10)
p1 = ax.barh(ypos, L1SoloList, color='#aed3e3', alpha=0.90, edgecolor='#000000', height=height, align='center')

p2 = ax.barh(ypos, L1TDList, color='#ed1f24', alpha=0.90, edgecolor='#000000', height=height, align='center',
             left=[i for i in L1SoloList])

p3 = ax.barh(ypos, AluList, color='#59bd7d', alpha=0.90, edgecolor='#000000', height=height, align='center',
             left=[i+j for i,j in zip(L1SoloList, L1TDList)])

p4 = ax.barh(ypos, SVAList, color='#faa41a', alpha=0.90, edgecolor='#000000', height=height, align='center',
             left=[i+j+x for i,j,x in zip(L1SoloList, L1TDList, AluList)])

p5 = ax.barh(ypos, PSDList, color='#8B4513', alpha=0.80, edgecolor='#000000', height=height, align='center',
             left=[i+j+x+z for i,j,x,z in zip(L1SoloList, L1TDList, AluList, SVAList)])

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

## Customize axis
plot_xmargin = 20
plot_ymargin = 4

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
nbEventsPerTumorList = nbEventsPerTumorTypeSortedSerie.values.tolist()

print "nbEventsPerTumorTypeSortedSerie: ", nbEventsPerTumorTypeSortedSerie

autolabel(p5, ax, nbEventsPerTumorList) ## autolabel function

## Make legend
circle1 = mpatches.Circle((0, 0), 5, color='#aed3e3', alpha=0.90)
circle2 = mpatches.Circle((0, 0), 5, color='#ed1f24', alpha=0.90)
circle3 = mpatches.Circle((0, 0), 5, color='#59bd7d', alpha=0.90)
circle4 = mpatches.Circle((0, 0), 5, color='#faa41a', alpha=0.90)
circle5 = mpatches.Circle((0, 0), 5, color='#8B4513', alpha=0.90)

l = plt.figlegend((circle1, circle2, circle3, circle4, circle5), ("L1-solo", "L1-transduction", "Alu", "SVA", "processed-pseudogene"), loc = 'upper center', ncol=5, labelspacing=0.75, fontsize=8, fancybox=True)

## Save figure
fileName = outDir + "/PCAWG_retrotransposition_events_tumorTypes.pdf"

plt.savefig(fileName)

############# PIE CHART ###############

## 1. Gather data

regionsList = []

## For project code
for projectCode in cohortObj.VCFdict:

    ## For donor
    for VCFObj in cohortObj.VCFdict[projectCode]:

        ## For MEI
        for MEIObj in VCFObj.lineList:

            if (MEIObj.infoDict['REGION']=="splicing") or (MEIObj.infoDict['REGION']=="upstream,downstream") or (MEIObj.infoDict['REGION']=="upstream") or (MEIObj.infoDict['REGION']=="downstream"):
                region = "Other"

            elif (MEIObj.infoDict['REGION']=="UTR5") or (MEIObj.infoDict['REGION']=="UTR3") or (MEIObj.infoDict['REGION']=="UTR5,UTR3") or (MEIObj.infoDict['REGION']=="UTR3,UTR5"):
                region = "UTR"

            elif (MEIObj.infoDict['REGION']=="ncRNA_exonic") or (MEIObj.infoDict['REGION']=="ncRNA_intronic") or (MEIObj.infoDict['REGION']=="ncRNA_splicing"):
                region = "ncRNA"

            else:
                region = MEIObj.infoDict['REGION']

            regionsList.append(region)

regionTuples = [(x, int(regionsList.count(x))) for x in set(regionsList)]
regionList =  [list(t) for t in zip(*regionTuples)]

labels = regionList[0]
sizes = regionList[1]

## 2. Make pie chart
fig = plt.figure(figsize=(6,6))
fig.suptitle('Somatic MEI functional spectrum', fontsize=16)
colors = ['#008000', '#A67D3D', '#87CEFA', '#ff0000', '#FFD700', '#FFA500']
patches, texts, perc = plt.pie(sizes, colors=colors, startangle=90, autopct='%1.1f%%', pctdistance=1.2, labeldistance=1)
plt.legend(patches, labels, loc="best", fontsize=11)

##### Save figure
fileName = outDir + "/PCAWG_somatic_funcSpectrum_piechart.pdf"
plt.savefig(fileName)

## End ##
print
print "***** Finished! *****"
print
