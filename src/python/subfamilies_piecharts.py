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

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.1f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct


#### MAIN ####

## Import modules ##
import argparse
import sys
import os.path
import time
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl

## Conf:
mpl.rcParams['patch.force_edgecolor'] = True

## Get user's input ##
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('counts', help='tsv with event counts per subfamily')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
counts = args.counts
outDir = args.outDir
scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "counts: ", counts
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 


## 1. Load input table:
#######################
header("1. Load input table")

countsDf = pd.read_csv(counts, header=0, index_col=0, sep='\t')

## 2. Make pie charts for L1:
##############################
#     subfamily  counts group
# L1         Ta    8608      Ta
# L1       Ta-0     273      Ta
# L1       Ta-1    1296      Ta
# L1     pre-Ta     367  pre-Ta

header("2. Make pie charts for L1")
countsL1Df = countsDf.loc[["L1"]]

### Subfamily ###
labels = countsL1Df['subfamily'].tolist()
counts = countsL1Df['counts'].tolist()

## make pie chart
fig1, ax1 = plt.subplots()
ax1.pie(counts, labels=labels, autopct=make_autopct(counts), startangle=90, pctdistance=0.85)

## draw circle
centre_circle = plt.Circle((0,0),0.60,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

## Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')  
plt.tight_layout()

## Save figure
fileName = outDir + "/L1_subfamilies.pdf"
plt.savefig(fileName)


## 3. Make pie chart for Alu:
##############################
#     subfamily  counts group
# L1         Ta    8608      Ta
# L1       Ta-0     273      Ta
# L1       Ta-1    1296      Ta
# L1     pre-Ta     367  pre-Ta

header("3. Make pie chart for Alu")

### Subfamily
countsAluDf = countsDf.loc[["Alu"]]
labels = countsAluDf['subfamily'].tolist()
counts = countsAluDf['counts'].tolist()

## make pie chart
fig1, ax1 = plt.subplots()
ax1.pie(counts, labels=labels, autopct=make_autopct(counts), startangle=90, pctdistance=0.85)

## draw circle
centre_circle = plt.Circle((0,0),0.60,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

## Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')  
plt.tight_layout()

## Save figure
fileName = outDir + "/Alu_subfamilies.pdf"
plt.savefig(fileName)

## 4. Make pie chart for SVA:
##############################
#     subfamily  counts group
# SVA     SVA_A       1   SVA
# SVA     SVA_D       5   SVA
# SVA     SVA_E       6   SVA
# SVA     SVA_F      11   SVA

header("4. Make pie chart for SVA")

### Subfamily
countsSVADf = countsDf.loc[["SVA"]]
labels = countsSVADf['subfamily'].tolist()
counts = countsSVADf['counts'].tolist()

## make pie chart
fig1, ax1 = plt.subplots()
ax1.pie(counts, labels=labels, autopct=make_autopct(counts), startangle=90, pctdistance=0.85)

## draw circle
centre_circle = plt.Circle((0,0),0.60,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

## Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')  
plt.tight_layout()

## Save figure
fileName = outDir + "/SVA_subfamilies.pdf"
plt.savefig(fileName)

####
header("Finished")

