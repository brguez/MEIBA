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

#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import time
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('lenghts', help='L1-mediated deletion size ditribution')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
lenghts = args.lenghts
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "lenghts: ", lenghts
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"


## Start ## 

## Generate array with deletion lengths
########################################
lenghtsFile = open(lenghts, 'r')

lenghtList = []

# For each line
for line in lenghtsFile: 

    line = line.rstrip('\r\n')
    line = line.split('\t')

    if line[0] != "UNK":
        lenghtList.append(int(line[0]))


print len(lenghtList), lenghtList

## Make plots
#############
sns.set(style="whitegrid")

### Striplot
fig = plt.figure(figsize=(2,4))
#fig.suptitle('Variant allele frequencies', fontsize=20)

ax = sns.stripplot(lenghtList, size=8, jitter=True, orient="v", edgecolor="black", color='#ff0000', alpha=1, linewidth=1)

ax.set_yscale('log')

## Save figure
outFile = outDir + '/' + "L1DEL_size_distr.striplot.pdf"
fig.savefig(outFile)

### Histogram
fig = plt.figure(figsize=(4,4))
plt.hist(lenghtList, range=[0, 8000], bins=15, color='#ff0000', alpha=1)
plt.yticks(np.arange(0, 21, 2))
plt.xticks(np.arange(0, 8001, 1000))

## Save figure
outFile = outDir + '/' + "L1DEL_size_distr.hist.pdf"
fig.savefig(outFile)

#sns.distplot(lenghtList, bins=range(1, , 10), kde=False)
#sns.distplot(a, bins=range(1, 110, 10), ax=ax, kde=False)
#plt.xscale('log')




## Customize ticks
# X axis
#xPosList = [ 500, 1000, 10000, 100000, 1000000, 10000000, max(lenghtList) ]
#ax.set_xticks(xPosList)
#ax.set_xticklabels(xPosList)
#locs, labels = plt.xticks()
#plt.setp(labels, rotation=90)

#plt.hist(lenghtList, bins=40, color='#008000', alpha=0.75)
#plt.xlabel("VAF", fontsize=14)
#plt.ylabel("# MEI", fontsize=14)
#plt.xlim(0, 1)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
#ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)
#ax.set_axisbelow(True)

## Customize ticks
#plt.xticks(np.arange(0, 1.01, 0.1))
#locs, labels = plt.xticks()
#plt.setp(labels, rotation=30)


#### END
header("Finished")


