#!/usr/bin/env python
#coding: utf-8

### CLASSES ###

class GTF:
    """

    """
    
    def __init__(self, line):
        """

        """
        line = line.rstrip('\n\r')
        line = line.split("\t")

        self.chrom = line[0]
        self.source = line[1] 
        self.feature = line[2] 
        self.beg = line[3]
        self.end = line[4]
        self.score = line[5] 
        self.strand = line[6] 
        self.frame = line[7]
        self.attribute = line[8]
        self.attributeDict = self.parse_attribute()

    def parse_attribute(self):
        """

        """
        attributeList = self.attribute.split(" ")
        attributeDict = {}

        tagList = []

        for i in xrange(0,len(attributeList),2):
            key = attributeList[i]
            value = attributeList[i+1] #"ENSG00000139618.10";
            value = value.replace("\"", "") # remove quotes
            value = value.replace(";", "") # remove comma

            if (key == "tag"):            
                tagList.append(value)
            else:
                attributeDict[key] = value
        
        attributeDict["tags"] = tagList

        return attributeDict

class gene(GTF):
    """

    """

    def __init__(self, line):
        """

        """
        GTF.__init__(self, line)
        self.representativeTranscript = ""
        self.transcriptList = []

    def addTranscript(self, transcriptObj):
    
        self.transcriptList.append(transcriptObj)

class transcript(GTF):
    """

    """

    def __init__(self, line):
        """

        """
        GTF.__init__(self, line)
        exonList = []
         


## Load modules/libraries
import sys
import argparse
import os
import errno

## Get user's input
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('annotGTF', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
annotGTFpath = args.annotGTF
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "annotation: ", annotGTFpath
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print

##### 1) Generate 3 independent GTF for genes, transcripts and exons, respectively.  
annotation = open(annotGTFpath, 'r')
allPseudogenesDict = {}

## Create output files
genesPath = outDir + "/genes.gtf"
transcriptsPath = outDir + "/transcripts.gtf"
exonsPath = outDir + "/exons.gtf"

genes = open(genesPath, 'w')
transcripts = open(transcriptsPath, 'w')
exons = open(exonsPath, 'w')

# Read file line by line
for line in annotation:
    
    ## Discard header
    if not line.startswith("#"):
        
        feature = line.split("\t")[2]

        # a) Gene
        if (feature == "gene"):
            genes.write(line)

        # b) Transcript
        elif (feature == "transcript"):
            transcripts.write(line)

        # c) Exon
        elif (feature == "exon"):
            exons.write(line)
    
##### 2) Read genes GTF and generate a gene object for each protein coding gene
## Save gene objects in a dictionary using the gene_name as key and gene_object as value
## Gene selection criteria:
# - Protein coding biotype

genes = open(genesPath, 'r')

genesDict = {} 

# Read file line by line
for line in genes:
    
    geneObj = gene(line)        
    geneName = geneObj.attributeDict["gene_name"]

    ## Select only protein coding genes
    if (geneObj.attributeDict["gene_type"] == "protein_coding"):
        
        genesDict[geneName] = geneObj
        
    
print "Nb.genes: ",len(genesDict)


##### 2) Add transcripts to the gene objects
## Transcript selection criteria:
# - 'protein_coding' biotype
# - 'appris_principal' tag. Transcript expected to code for the main functional isoform based solely on the core modules in the APPRIS database.
# - 'basic' tag. Identifies a subset of representative transcripts for each gene; prioritises full-length protein coding transcripts over partial or non-protein coding transcripts within the same gene, and intends to highlight those transcripts that will be useful to the majority of users.

transcripts = open(transcriptsPath, 'r')

# Read file line by line
for line in transcripts:

    transcriptObj = transcript(line)        


    # Select transcripts from protein coding genes and tagged with 'appris_principal' and 'basic' tags
    if (transcriptObj.attributeDict["gene_type"] == "protein_coding") and ('appris_principal' in  transcriptObj.attributeDict['tags']) and ('basic' in  transcriptObj.attributeDict['tags']):
        
        geneName = transcriptObj.attributeDict["gene_name"]
        
        # Expected gene name
        if geneName in genesDict:

            genesDict[geneName].addTranscript(transcriptObj)
            
        # Unexpected gene name
        else:
            print "WARNING! Expected gene name"
        
   
##### 3) Per gene object, select a set of representative transcripts. 
# The representative transcript will be the ones with lowest level. If still several possibilities, pick the first one in the list arbitrarialy. 



        
print "***** Finished! *****"
print

