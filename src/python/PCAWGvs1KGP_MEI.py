#!/usr/bin/env python
#coding: utf-8

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
    print timeInfo, string,

#### CLASSES ####
class germlineMEIDb():
    """
            .....................

            Methods:
            -

    """
    def __init__(self):
        """

        """
        self.germlineMEIDict = {}

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
                pos = line[1]
                MEIClass = line[3]
                family = line[4]
                TSD = line[5]
                orientation = line[6]
                length = line [7]
                germlineMEIObj = germlineMEI(chrom, pos, MEIClass, family, TSD, orientation, length)

                ### Create a nested dictionary containing for each insertion class and chromosome a dictionary of
                ## germline MEI objects organized by coordinates
                # key1(MEIclass) -> value(dict2) -> key2(chrId) -> value(dict3) -> key3(coordinate) -> value(germlineMEIObj)

                ## Create nested dictionaries as needed
                # A) First MEI of a given class
                if germlineMEIObj.MEIClass not in self.germlineMEIDict:
                    self.germlineMEIDict[germlineMEIObj.MEIClass] = {}
                    self.germlineMEIDict[germlineMEIObj.MEIClass][germlineMEIObj.chrom] = {}

                # B) First MEI of a given class in the chromosome
                elif germlineMEIObj.chrom not in self.germlineMEIDict[germlineMEIObj.MEIClass]:
                    self.germlineMEIDict[germlineMEIObj.MEIClass][germlineMEIObj.chrom] = {}

                ## Add germline MEI object to the dictionary
                self.germlineMEIDict[germlineMEIObj.MEIClass][germlineMEIObj.chrom][germlineMEIObj.pos] = germlineMEIObj

        # Close germline MEI database file
        inputBedFile.close()

class germlineMEI():
    """
            .....................

            Methods:
            -

    """
    def __init__(self, chrom, pos, MEIClass, family, TSD, orientation, length):
        """

        """
        self.chrom = chrom
        self.pos = pos
        self.MEIClass = MEIClass
        self.family = family
        self.TSD = TSD
        self.orientation = orientation
        self.length = length


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import formats
import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import matplotlib.patches as mpatches
import scipy.stats as stats
from operator import itemgetter, attrgetter, methodcaller


## Get user's input ##
parser = argparse.ArgumentParser(description= "Compare PCAWG with 1KGP MEI and obtain different statistics")
parser.add_argument('VCF', help='PCAWG MEI')
parser.add_argument('oneKGPdb', help='BED file containing 1KGP MEI')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
oneKGPdb = args.oneKGPdb
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "oneKGPdb: ", oneKGPdb
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

###########
## Start ## 
###########

############################################
## 1. Create VCF object and read input VCF #
############################################
header("Read input VCF")
VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)

###########################################
## 2. Create germline MEI database object #
###########################################
header("Read germline MEI database")
germlineMEIDbObj = germlineMEIDb()
germlineMEIDbObj.read_bed(oneKGPdb)

##################################################################
## 3. Compare TraFiC MEI characteristics with matching 1KGP MEI #
##################################################################
header("Compare PCAWG with 1KGP MEI")

# Gather all this information into arrays for plotting.
bkpDistances = []
bkpDistL1 = []
bkpDistAlu = []
bkpDistSVA = []
strandComparisons = []
TSDSeqComparisons = []
TSDLenComparisons = []
MEILenL1 = []
MEILenL1Distances = []
MEILenAlu = []
MEILenAluDistances = []
MEILenSVA = []
MEILenSVADistances = []

MEIObjList = VCFObj.lineList

for MEIObj in MEIObjList:

    ## Select only those MEIs that passes the filters
    if (MEIObj.filter == "PASS"):

        info("Checking if " + MEIObj.chrom + ":" + str(MEIObj.pos) + " MEI is in 1KGP germline database..." )

        # A) MEI already reported in 1KGP
        if ( 'GERMDB' in MEIObj.infoDict ):

            print 'yes'

            ## Make numpy array with the chromosomal positions for the germline MEI on the same chromosome in the 1KGP database
            MEIDbPosArr = np.array(sorted(germlineMEIDbObj.germlineMEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom].keys(), key=int), dtype='int')

            ## Identify those MEI in the germline database overlapping the first insertion breakpoint
            # Determine bkpA beg and end range to search for overlap with 1KGP MEI database
            begBkpA = MEIObj.pos - int(MEIObj.infoDict["CIPOS"]) - 10
            endBkpA = MEIObj.pos + int(MEIObj.infoDict["CIPOS"]) + 10

            # Select those MEI in the database within bkpA beg and end range
            filteredArrBkpA = MEIDbPosArr[(MEIDbPosArr >= begBkpA) & (MEIDbPosArr <= endBkpA)]

            ## A) Target site info available -> Identify those MEI in the germline database overlapping the second insertion breakpoint
            if 'TSLEN' in MEIObj.infoDict:

                # a) Consistent TSD
                if (MEIObj.infoDict["TSLEN"] != "inconsistent"):

                    # Determine second breakpoint coordinates and bkpB beg and end range to search for overlap with 1KGP MEI database
                    bkpB = MEIObj.pos + abs(int(MEIObj.infoDict["TSLEN"]))
                    begBkpB = bkpB - int(MEIObj.infoDict["CIPOS"]) - 10
                    endBkpB = bkpB + int(MEIObj.infoDict["CIPOS"]) + 10

                    # Select those MEI in the database within bkpA beg and end range
                    filteredArrBkpB = MEIDbPosArr[(MEIDbPosArr >= begBkpB) & (MEIDbPosArr <= endBkpB)]

                # b) Inconsistent TSD
                else:

                    # Empty array
                    filteredArrBkpB = np.array([ ])

            ## B) Target site info not available
            else:
                # Empty array
                filteredArrBkpB = np.array([ ])

            ## Make a single array containing a not redundant list of 1KGP MEI coordinates overlapping the first and/or second insertion breakpoints.
            filteredArr = np.array(np.unique(np.concatenate((filteredArrBkpA, filteredArrBkpB), axis=0)), dtype='int')

            ### Compare insertion characteristics (orientation, TSD, length...)
            ## A) Single 1KGP MEI overlapping insertion breakpoint
            if (len(filteredArr) == 1):
                MEI1KGP = germlineMEIDbObj.germlineMEIDict[MEIObj.infoDict["CLASS"]][MEIObj.chrom][str(filteredArr[0])]
    
                # a) Breakpoint position
                if ('TSLEN' in MEIObj.infoDict):
                    bkpB = abs(MEIObj.pos + abs(int(MEIObj.infoDict["TSLEN"])))
                    diffA = abs(MEIObj.pos - int(MEI1KGP.pos))
                    diffB = abs(bkpB - int(MEI1KGP.pos))
                    difference = min(diffA, diffB)
                    bkpDistances.append(difference)
    
                    if (MEIObj.infoDict['CLASS'] == "L1"):
                        bkpDistL1.append(difference)
    
                    elif (MEIObj.infoDict['CLASS'] == "Alu"):
                        bkpDistAlu.append(difference)
    
                    elif (MEIObj.infoDict['CLASS'] == "SVA"):
                        bkpDistSVA.append(difference)
    
                # b) Orientation (+ or -)
                if (MEI1KGP.orientation != "UNK"):
                    comparisonResult = 'consistent' if MEIObj.infoDict["STRAND"] == MEI1KGP.orientation else 'inconsistent'
                    strandComparisons.append(comparisonResult)
    
                # c) Target site duplication length
                if (MEI1KGP.TSD != "UNK") and ('TSSEQ' in MEIObj.infoDict):
    
                    # Compare sequence
                    comparisonResult = 'consistent' if MEIObj.infoDict["TSSEQ"] == MEI1KGP.TSD else 'inconsistent'
                    TSDSeqComparisons.append(comparisonResult)
    
                    # Compare length
                    comparisonResult = 'consistent' if int(MEIObj.infoDict["TSLEN"]) == len(MEI1KGP.TSD) else 'inconsistent'
                    TSDLenComparisons.append(comparisonResult)
    
                # d) transposon length
                if (MEI1KGP.length != "UNK") and ('LEN' in MEIObj.infoDict):
    
                    length = int(MEIObj.infoDict['LEN'])
                    difference = abs(int(MEI1KGP.length) - length)
                    lenTuple = ( length, int(MEI1KGP.length) )
    
                    if (MEIObj.infoDict['CLASS'] == "L1"):
                        MEILenL1Distances.append(difference)
                        MEILenL1.append(lenTuple)
    
                    elif (MEIObj.infoDict['CLASS'] == "Alu"):
                        MEILenAluDistances.append(difference)
                        MEILenAlu.append(lenTuple)
    
                    elif (MEIObj.infoDict['CLASS'] == "SVA"):
                        MEILenSVADistances.append(difference)
                        MEILenSVA.append(lenTuple)

            ## B) Multiple 1KGP MEI overlapping insertion breakpoint -> Not consider bkp
            else:
                print "skip-ambiguous-MEI"

        # B) Novel MEI
        else:
            print 'no'

#######################################
## 3. Plot information gathered in 2) #
#######################################

## 3.1 Bkp distances, strand and TSD consistency/inconsistency bar plot
########################################################################
header("Make barplot")

## 3.1.1 Gather information
## Bkp distances
# Count number of consistent and inconsistent bkp instances
# Consistent bkp when exact breakpoint identified
bkpCategories = [ "consistent" if (x==0) else "inconsistent" for x in bkpDistances ]
bkpDict = dict([(x, bkpCategories.count(x)) for x in set(bkpCategories)])

# Convert into percentages
total = (bkpDict['consistent'] + bkpDict['inconsistent'])
bkpDict['consistent'] = float(bkpDict['consistent']) / total * 100
bkpDict['inconsistent'] = float(bkpDict['inconsistent']) / total * 100
#print "bkp-distances: ", sorted(bkpDistances), np.mean(bkpDistances), np.std(bkpDistances), bkpDict

## TSD seq
# Count number of consistent and inconsistent TSD sequence instances
TSDSeqDict = dict([(x, TSDSeqComparisons.count(x)) for x in set(TSDSeqComparisons)])

# Convert into percentages
total = (TSDSeqDict['consistent'] + TSDSeqDict['inconsistent'])
TSDSeqDict['consistent'] = float(TSDSeqDict['consistent']) / total * 100
TSDSeqDict['inconsistent'] = float(TSDSeqDict['inconsistent']) / total * 100

## TSD length
# Count number of consistent and inconsistent TSD length instances
TSDLenDict = dict([(x, TSDLenComparisons.count(x)) for x in set(TSDLenComparisons)])

# Convert into percentages
total = (TSDLenDict['consistent'] + TSDLenDict['inconsistent'])
TSDLenDict['consistent'] = float(TSDLenDict['consistent']) / total * 100
TSDLenDict['inconsistent'] = float(TSDLenDict['inconsistent']) / total * 100

## Strand
# Count number of consistent and inconsistent strand instances
strandDict = dict([(x, strandComparisons.count(x)) for x in set(strandComparisons)])

# Convert into percentages
total = (strandDict['consistent'] + strandDict['inconsistent'])
strandDict['consistent'] = float(strandDict['consistent']) / total * 100
strandDict['inconsistent'] = float(strandDict['inconsistent']) / total * 100

## 3.1.2 Put percentages into two lists
consistentList = [ bkpDict['consistent'], TSDLenDict['consistent'], strandDict['consistent'] ]
inconsistentList = [ bkpDict['inconsistent'], TSDLenDict['inconsistent'], strandDict['inconsistent'] ]

## 3.1.3 Make figure
# Make bar plot
xpos = np.arange(3)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(5,6))
fig.suptitle('MEI features consistency')
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
plt.xticks(xpos, ('BKP', 'TSD-LEN', 'STRAND'))

## Make legend
circle1 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.75)
circle2 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.75)
plt.figlegend((circle1, circle2), ('Consistent', 'Inconsistent'), loc = 'lower center', ncol=2, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/PCAWG_1KGP_barPlot.pdf"
plt.savefig(fileName)

## 3.2 Transposon length distances correlation plot
#####################################################
header("Compute correlations and make scatterplots")

fig = plt.figure(figsize=(10,14))
fig.suptitle('MEI length correlation', fontsize=18)

## L1
tmpList = map(list, zip(*MEILenL1))
L1lenPCAWG = tmpList[0]
L1len1KGP = tmpList[1]

# Compute correlation
corr = stats.pearsonr(L1lenPCAWG, L1len1KGP)
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

# Make scatterplot
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("LINE-1", fontsize=14)
plt.scatter(L1lenPCAWG, L1len1KGP, color='#008000', alpha=.4)
plt.xlim((0, (max(L1lenPCAWG) + 500)))
plt.ylim((0, (max(L1len1KGP) + 600)))
plt.xlabel('PCAWG', fontsize=12)
plt.ylabel('1KGP', fontsize=12)
ax1.text(0.5, 0.1, text, transform = ax1.transAxes)

## Alu
tmpList = map(list, zip(*MEILenAlu))
AlulenPCAWG = tmpList[0]
Alulen1KGP = tmpList[1]

# Compute correlation
corr = stats.pearsonr(AlulenPCAWG, Alulen1KGP )
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

# Make scatterplot
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("ALU", fontsize=14)
plt.scatter(AlulenPCAWG, Alulen1KGP, color='#008000', alpha=.4)
plt.xlim((0, (max(AlulenPCAWG) + 10)))
plt.ylim((0, (max(Alulen1KGP) + 15)))
plt.xlabel('PCAWG', fontsize=12)
plt.ylabel('1KGP', fontsize=12)
ax2.text(0.5, 0.1, text, transform = ax2.transAxes)

## SVA
tmpList = map(list, zip(*MEILenSVA))
SVAlenPCAWG = tmpList[0]
SVAlen1KGP = tmpList[1]

# Compute correlation
corr = stats.pearsonr(SVAlenPCAWG, SVAlen1KGP )
coefficient = format(corr[0], '.3f')
pvalue = corr[1]
text = 'pearson_corr: ' + str(coefficient)

# Make scatterplot
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=14)
plt.scatter(SVAlenPCAWG, SVAlen1KGP, color='#008000', alpha=.4)
plt.xlim((0, (max(SVAlenPCAWG) + 100)))
plt.ylim((0, (max(SVAlen1KGP) + 150)))
plt.xlabel('PCAWG', fontsize=12)
plt.ylabel('1KGP', fontsize=12)
ax3.text(0.5, 0.1, text, transform = ax3.transAxes)

## Save figure
fileName = outDir + "/PCAWG_1KGP_correlation.pdf"
plt.savefig(fileName)

## 3.3 Venn diagram
######################
header("Make venn diagrams")

# Generate needed data
classesList1KGP = [line.rstrip('\n').split('\t')[3] for line in open(oneKGPdb)]

# Make figure
fig = plt.figure(figsize=(5,6))
#fig.suptitle('Overlap between MEI calls')

## Alu
# Generate needed data
AluPCAWG =  sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'Alu') and (MEIObj.filter == "PASS") and ('GERMDB' not in MEIObj.infoDict))
AluBoth = sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'Alu') and (MEIObj.filter == "PASS") and ('GERMDB' in MEIObj.infoDict))
totalAlu1KGP = classesList1KGP.count("Alu")
Alu1KGP = totalAlu1KGP - AluBoth

# Make venn diagram
ax1 = fig.add_subplot(2, 2, 2)
ax1.set_title("ALU", fontsize=10)
venn2(subsets=(AluPCAWG, Alu1KGP, AluBoth), set_labels = ('', ''))


## L1
# Generate needed data
L1PCAWG =  sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'L1') and (MEIObj.filter == "PASS") and ('GERMDB' not in MEIObj.infoDict))
L1Both = sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'L1') and (MEIObj.filter == "PASS") and ('GERMDB' in MEIObj.infoDict))
totalL11KGP = classesList1KGP.count("L1")
L11KGP = totalL11KGP - L1Both

# Make venn diagram
ax2 = fig.add_subplot(2, 2, 1)
ax2.set_title("LINE-1", fontsize=10)
venn2(subsets=(L1PCAWG, L11KGP, L1Both), set_labels = ('', ''))

## SVA
# Generate needed data
SvaPCAWG =  sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'SVA') and (MEIObj.filter == "PASS") and ('GERMDB' not in MEIObj.infoDict))
SvaBoth = sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'SVA') and (MEIObj.filter == "PASS") and ('GERMDB' in MEIObj.infoDict))
totalSva1KGP = classesList1KGP.count("SVA")
Sva1KGP = totalSva1KGP - SvaBoth

# Make venn diagram
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("SVA", fontsize=10)
venn2(subsets=(SvaPCAWG, Sva1KGP, SvaBoth), set_labels = ('', ''))

## ERVK
# Generate needed data
ERVKPCAWG =  sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'ERVK') and (MEIObj.filter == "PASS") and ('GERMDB' not in MEIObj.infoDict))
ERVKBoth = sum(1 for MEIObj in MEIObjList if (MEIObj.infoDict['CLASS'] == 'ERVK') and (MEIObj.filter == "PASS") and ('GERMDB' in MEIObj.infoDict))
totalERVK1KGP = classesList1KGP.count("ERVK")
ERVK1KGP = totalERVK1KGP - ERVKBoth

# Make venn diagram
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("ERVK", fontsize=10)
venn2(subsets=(ERVKPCAWG, ERVK1KGP, ERVKBoth), set_labels = ('', ''))

## Make legend
circle1 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.5)
circle2 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.5)
plt.figlegend((circle1, circle2), ('PCAWG', '1KGP'), loc = 'upper center', ncol=2, labelspacing=0., fontsize=8)

## Save figure
fileName = outDir + "/PCAWG_1KGP_venn.pdf"
plt.savefig(fileName)

#########################################################
## 4. Compute bkp distance mean and standard desviation #
#########################################################

header("Compute bkp distance mean and standard desviation")

#### Alu
print "ALU-mean: ", np.mean(bkpDistAlu)
print "ALU-std: ", np.std(bkpDistAlu)

#### L1
print "L1-mean: ", np.mean(bkpDistL1)
print "L1-std: ",np.std(bkpDistL1)

#### SVA
print "SVA-mean: ", np.mean(bkpDistSVA)
print "SVA-std: ", np.std(bkpDistSVA)


## End ##
print
print "***** Finished! *****"
print
