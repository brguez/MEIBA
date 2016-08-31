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
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('VCF', help='...')
parser.add_argument('germlineMEIBed', help='...')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
VCF = args.VCF
germlineMEIBed = args.germlineMEIBed
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "VCF: ", VCF
print "germlineMEIBed: ", germlineMEIBed
print "outDir: ", outDir
print 
print "***** Executing ", scriptName, ".... *****"
print 


## Start ## 

## 1. Create VCF object and read input VCF
VCFObj = formats.VCF()
VCFObj.read_VCF(VCF)

## 2. Create germline MEI database object
germlineMEIDbObj = germlineMEIDb()
germlineMEIDbObj.read_bed(germlineMEIBed)

## 3. Compare TraFiC MEI characteristics with matching 1000G MEI
# Gather all this information into arrays for plotting. 
bkpDistances = [] 
strandComparisons = []
TSDComparisons = [] 
TSDLenDistances = []
MEILenL1 = []
MEILenL1Distances = []
MEILenAlu = []
MEILenAluDistances = []
MEILenSVA = []
MEILenSVADistances = []

VCFlineObjList = VCFObj.lineList

for VCFlineObj in VCFlineObjList:

	info("Checking if " + VCFlineObj.chrom + ":" + str(VCFlineObj.pos) + " MEI is in 1000 genomes germline database..." )

	# A) MEI already reported in 1000 Genomes
	if ( 'GERMDB' in VCFlineObj.infoDict ):

		print 'yes'

		## Make numpy array with the chromosomal positions for the germline MEI on the same chromosome in the 1000 genomes database
		MEIDbPosArr = np.array(sorted(germlineMEIDbObj.germlineMEIDict[VCFlineObj.infoDict["CLASS"]][VCFlineObj.chrom].keys(), key=int), dtype='int')

		## Identify those MEI in the germline database overlapping the first insertion breakpoint
		# Determine bkpA beg and end range to search for overlap with 1000 genomes MEI database
		begBkpA = VCFlineObj.pos - int(VCFlineObj.infoDict["CIPOS"]) - 10
		endBkpA = VCFlineObj.pos + int(VCFlineObj.infoDict["CIPOS"]) + 10
	
		# Select those MEI in the database within bkpA beg and end range
		filteredArrBkpA = MEIDbPosArr[(MEIDbPosArr >= begBkpA) & (MEIDbPosArr <= endBkpA)] 

		## A) Target site info available -> Identify those MEI in the germline database overlapping the second insertion breakpoint
		if 'TSLEN' in VCFlineObj.infoDict:
	
			# a) Consistent TSD
			if (VCFlineObj.infoDict["TSLEN"] != "inconsistent"):	
				# Determine second breakpoint coordinates and bkpB beg and end range to search for overlap with 1000 genomes MEI database
				bkpB = VCFlineObj.pos + abs(int(VCFlineObj.infoDict["TSLEN"]))
				begBkpB = bkpB - int(VCFlineObj.infoDict["CIPOS"]) - 10
				endBkpB = bkpB + int(VCFlineObj.infoDict["CIPOS"]) + 10
	
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
	
		## Make a single array containing a not redundant list of 1000 genomes MEI coordinates overlapping the first and/or second insertion breakpoints.
		filteredArr = np.array(np.unique(np.concatenate((filteredArrBkpA, filteredArrBkpB), axis=0)), dtype='int')

		### Compare insertion characteristics (orientation, TSD, length...)
		## A) Single 1000G MEI overlapping insertion breakpoint 
		if (len(filteredArr) == 1):
			MEI1000G = germlineMEIDbObj.germlineMEIDict[VCFlineObj.infoDict["CLASS"]][VCFlineObj.chrom][str(filteredArr[0])]
			
			# a) Breakpoint position
			if ('TSLEN' in VCFlineObj.infoDict):	
				bkpB = abs(VCFlineObj.pos + abs(int(VCFlineObj.infoDict["TSLEN"])))
				diffA = abs(VCFlineObj.pos - int(MEI1000G.pos))
				diffB = abs(bkpB - int(MEI1000G.pos)) 
				difference = min(diffA, diffB)
				bkpDistances.append(difference)

			# b) Orientation (+ or -)
			if (MEI1000G.orientation != "UNK"):
				comparisonResult = 'consistent' if VCFlineObj.infoDict["STRAND"] == MEI1000G.orientation else 'inconsistent'
				strandComparisons.append(comparisonResult)
			
			# c) Target site duplication length
			if (MEI1000G.TSD != "UNK") and ('TSSEQ' in VCFlineObj.infoDict):
				
				# Compare sequence
				comparisonResult = 'consistent' if VCFlineObj.infoDict["TSSEQ"] == MEI1000G.TSD else 'inconsistent'
				TSDComparisons.append(comparisonResult)

				# Compare length
				length = abs(int(VCFlineObj.infoDict['TSLEN']))				
				difference = abs(len(MEI1000G.TSD) - length)	
				TSDLenDistances.append(difference)
			
			# d) transposon length 
			if (MEI1000G.length != "UNK") and ('LEN' in VCFlineObj.infoDict):

				## note: substract 31, since it is the length of the polyA sequence added to the  Alu 
				# and SVA reference sequences. For L1, the length is 33... This is temporal. I need to modify 'findInsertionBkp' script
				# to substract polyA length internally. 
				length = int(VCFlineObj.infoDict['LEN']) - 32		
				difference = abs(int(MEI1000G.length) - length) 
				lenTuple = ( length, MEI1000G.length )

				if (VCFlineObj.infoDict['CLASS'] == "L1"):			
					MEILenL1Distances.append(difference)				
					MEILenL1.append(lenTuple)

				elif (VCFlineObj.infoDict['CLASS'] == "Alu"):
					MEILenAluDistances.append(difference)
					MEILenAlu.append(lenTuple)
				
				elif (VCFlineObj.infoDict['CLASS'] == "SVA"):
					MEILenSVADistances.append(difference)		
					MEILenSVA.append(lenTuple)

		## B) Multiple 1000G MEI overlapping insertion breakpoint -> Not consider bkp
		else:
			print "skip-ambiguous-MEI"

	# B) Novel MEI	
	else:
		print 'no'		

## 3. Plot information gathered in 2)
## 3.1 Bkp distances, strand and TSD consistency/inconsistency bar plot
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

## TSD 
# Count number of consistent and inconsistent TSD length instances
TSDDict = dict([(x, TSDComparisons.count(x)) for x in set(TSDComparisons)])

# Convert into percentages
total = (TSDDict['consistent'] + TSDDict['inconsistent'])
TSDDict['consistent'] = float(TSDDict['consistent']) / total * 100
TSDDict['inconsistent'] = float(TSDDict['inconsistent']) / total * 100

## Strand 
# Count number of consistent and inconsistent strand instances
strandDict = dict([(x, strandComparisons.count(x)) for x in set(strandComparisons)])

# Convert into percentages
total = (strandDict['consistent'] + strandDict['inconsistent'])
strandDict['consistent'] = float(strandDict['consistent']) / total * 100
strandDict['inconsistent'] = float(strandDict['inconsistent']) / total * 100

## 3.1.2 Put percentages into two lists
consistentList = [ bkpDict['consistent'], TSDDict['consistent'], strandDict['consistent'] ]
inconsistentList = [ bkpDict['inconsistent'], TSDDict['inconsistent'], strandDict['inconsistent'] ]
#print "TSD-len-distances: ", sorted(TSDLenDistances), np.mean(TSDLenDistances), np.std(TSDLenDistances), TSDDict

print consistentList, inconsistentList

## 3.1.3 Make figure
# Make bar plot 
xpos = np.arange(3)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fig = plt.figure(figsize=(5,6))
fig.suptitle('PCAWG and 1K-GENOMES consistency')
plt.ylabel('%', fontsize=12)
ax = fig.add_subplot(111)
p1 = ax.bar(xpos, consistentList, color='#008000', alpha=0.5, edgecolor='#000000', width=width, align='center')
p2 = ax.bar(xpos, inconsistentList, color='#ff0000', alpha=0.5, edgecolor='#000000', width=width, align='center',
             bottom=consistentList)

## Customize ticks
plt.yticks(np.arange(0, 101, 10))
plt.xticks(xpos, ('BKP', 'TSD', 'STRAND'))

## Make legend 
circle1 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.5)
circle2 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.5)
plt.figlegend((circle1, circle2), ('Consistent', 'Inconsistent'), loc = 'lower center', ncol=2, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/PCAWG_1KGENOMES_barPlot.pdf"
plt.savefig(fileName)

## 3.2 Transposon length distances correlation plot
fig = plt.figure(figsize=(5,6))  

## L1
# print "L1-len-distances: ", sorted(MEILenL1Distances), np.mean(MEILenL1Distances), np.std(MEILenL1Distances)
tmpList = map(list, zip(*MEILenL1))
L1lenPCAWG = tmpList[0]
L1len1000G = tmpList[1]

## Alu
# print "Alu-len-distances: ", sorted(MEILenAluDistances), np.mean(MEILenAluDistances), np.std(MEILenAluDistances)
tmpList = map(list, zip(*MEILenAlu))
AlulenPCAWG = tmpList[0]
Alulen1000G = tmpList[1]

plt.scatter(AlulenPCAWG, Alulen1000G)
plt.show()
print "Alu-corr: ", stats.pearsonr(AlulenPCAWG, Alulen1000G)

## SVA 
# print "SVA-len-distances: ", sorted(MEILenSVADistances), np.mean(MEILenSVADistances), np.std(MEILenSVADistances)
tmpList = map(list, zip(*MEILenSVA))
SVAlenPCAWG = tmpList[0]
SVAlen1000G = tmpList[1]

## 3.3 Venn diagram
info("Plot PCAWG vs 1KGENOMES MEI venn diagrams\n")

# Generate needed data
classesList1000G = [line.rstrip('\n').split('\t')[3] for line in open(germlineMEIBed)]

# Make figure
fig = plt.figure(figsize=(5,6))                
fig.suptitle('Overlap between MEI germline calls')

## L1
# Generate needed data
L1PCAWG =  sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'L1') and ('GERMDB' not in VCFlineObj.infoDict))
L1Both = sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'L1') and ('GERMDB' in VCFlineObj.infoDict))
totalL11000G = classesList1000G.count("L1")
L11000G = totalL11000G - L1Both

# Plot venn diagram
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("LINE-1", fontsize=10)
venn2(subsets=(L1PCAWG, L11000G, L1Both), set_labels = ('', ''))

## Alu
# Generate needed data
AluPCAWG =  sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'Alu') and ('GERMDB' not in VCFlineObj.infoDict))
AluBoth = sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'Alu') and ('GERMDB' in VCFlineObj.infoDict))
totalAlu1000G = classesList1000G.count("Alu")
Alu1000G = totalAlu1000G - AluBoth

# Plot venn diagram
ax1 = fig.add_subplot(2, 2, 2)
ax1.set_title("ALU", fontsize=10)
venn2(subsets=(AluPCAWG, Alu1000G, AluBoth), set_labels = ('', ''))

## SVA
# Generate needed data
SvaPCAWG =  sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'SVA') and ('GERMDB' not in VCFlineObj.infoDict))
SvaBoth = sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'SVA') and ('GERMDB' in VCFlineObj.infoDict))
totalSva1000G = classesList1000G.count("SVA")
Sva1000G = totalSva1000G - SvaBoth

# Plot venn diagram
ax1 = fig.add_subplot(2, 2, 3)
ax1.set_title("SVA", fontsize=10)
venn2(subsets=(SvaPCAWG, Sva1000G, SvaBoth), set_labels = ('', ''))

## ERVK
# Generate needed data
ERVKPCAWG =  sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'ERVK') and ('GERMDB' not in VCFlineObj.infoDict))
ERVKBoth = sum(1 for VCFlineObj in VCFlineObjList if (VCFlineObj.infoDict['CLASS'] == 'ERVK') and ('GERMDB' in VCFlineObj.infoDict))
totalERVK1000G = classesList1000G.count("ERVK")
ERVK1000G = totalERVK1000G - ERVKBoth

# Plot venn diagram
ax1 = fig.add_subplot(2, 2, 4)
ax1.set_title("ERVK", fontsize=10)
venn2(subsets=(ERVKPCAWG, ERVK1000G, ERVKBoth), set_labels = ('', ''))

## Make legend 
circle1 = mpatches.Circle((0, 0), 5, color='#ff0000', alpha=0.5)
circle2 = mpatches.Circle((0, 0), 5, color='#008000', alpha=0.5)
plt.figlegend((circle1, circle2), ('PCAWG', '1K-GENOMES'), loc = 'lower center', ncol=2, labelspacing=0., fontsize=8, fancybox=True )

## Save figure
fileName = outDir + "/PCAWG_1KGENOMES_venn.pdf"
plt.savefig(fileName)

## End ##
print 
print "***** Finished! *****"
print 


