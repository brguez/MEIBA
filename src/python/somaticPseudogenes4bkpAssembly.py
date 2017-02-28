#!/usr/bin/env python
#coding: utf-8

## Load modules/libraries
import sys
import argparse
import os
import errno

## Get user's input
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('insertions', help='')
parser.add_argument('metadata', help='')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
insertionsPath = args.insertions
metadataPath = args.metadata
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "insertions: ", insertionsPath
print "metadata: ", metadataPath
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print




### 1) Read insertions file and organize data into a dictionary with the following format:
# Key1 [tumor_wgs_aliquot_id] -> value [list of pseudogene dictionaries: pseudogeneDict1, pseudogeneDict2, ...]
# Note: The number of elements in the list will be equal to the total number of pseudogene insertions in a given tumor sample. 
# Note: The pseudogene dictionary contains all the info regarding a given pseudogene insertion

insertions = open(insertionsPath, 'r')
allPseudogenesDict = {}

# Read file line by line
for line in insertions:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
        
        # Generic fields
        chromPlus = fieldsList[0]
        begPlus = fieldsList[1] 
        endPlus = fieldsList[2]  
        nbReadsPlus = fieldsList[3]
        classPlus = "NA"
        readListPlus = fieldsList[5]
        chromMinus = fieldsList[6]
        begMinus = fieldsList[7] 
        endMinus = fieldsList[8]   
        nbReadsMinus = fieldsList[9] 
        classMinus = "NA"
        readListMinus = fieldsList[11]
        insertionType = "PSD"
        
        # L1 transductions specific fields
        chromSource = "NA" 
        begSource = "NA"
        endSource = "NA"
        strandSource = "NA"
        tdBeg = "NA"
        tdEnd = "NA"
        tdRnaLen = "NA"
        tdLen = "NA"

        # Processed pseudogene specific fields (PSD)
        psdGene = fieldsList[23]
        chromExonA = fieldsList[24]
        begExonA = fieldsList[25]
        endExonA = fieldsList[26]
        chromExonB = fieldsList[27]
        begExonB = fieldsList[28]
        endExonB = fieldsList[29]

        # Tumor aliquot Id
        tumor_wgs_aliquot_id = fieldsList[14]

        # Gather all the relevant pseudogene info in a dictionary
        pseudogeneDict = {}
        pseudogeneDict["chromPlus"] = chromPlus 
        pseudogeneDict["begPlus"] = begPlus  
        pseudogeneDict["endPlus"] = endPlus   
        pseudogeneDict["nbReadsPlus"] = nbReadsPlus 
        pseudogeneDict["classPlus"] = classPlus 
        pseudogeneDict["readListPlus"] = readListPlus 
        pseudogeneDict["chromMinus"] = chromMinus 
        pseudogeneDict["begMinus"] = begMinus  
        pseudogeneDict["endMinus"] = endMinus    
        pseudogeneDict["nbReadsMinus"] = nbReadsMinus  
        pseudogeneDict["classMinus"] = classMinus 
        pseudogeneDict["readListMinus"] = readListMinus 
        pseudogeneDict["insertionType"] = insertionType 
        pseudogeneDict["chromSource"] = chromSource  
        pseudogeneDict["begSource"] = begSource 
        pseudogeneDict["endSource"] = endSource 
        pseudogeneDict["strandSource"] = strandSource 
        pseudogeneDict["tdBeg"] = tdBeg 
        pseudogeneDict["tdEnd"] = tdEnd 
        pseudogeneDict["tdRnaLen"] = tdRnaLen 
        pseudogeneDict["tdLen"] = tdLen 
        pseudogeneDict["psdGene"] = psdGene
        pseudogeneDict["chromExonA"] = chromExonA
        pseudogeneDict["begExonA"] = begExonA
        pseudogeneDict["endExonA"] = endExonA
        pseudogeneDict["chromExonB"] = chromExonB
        pseudogeneDict["begExonB"] = begExonB
        pseudogeneDict["endExonB"] = endExonB

        if tumor_wgs_aliquot_id not in allPseudogenesDict:
            allPseudogenesDict[tumor_wgs_aliquot_id] = []
            
        allPseudogenesDict[tumor_wgs_aliquot_id].append(pseudogeneDict)
  
      
### 2) For each PCAWG tumor sample generate a properly formated file containing its somatic processed pseudogenes and
# organize these files into the following directory architecture:
# outDir -> dcc_project_code > submitted_donor_id -> tumor_wgs_aliquot_id -> **somatic pseudogenes tsv**

metadata = open(metadataPath, 'r')

# Read file line by line
for line in metadata:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        submitted_donor_id = fieldsList[0]	
        wgs_exclusion_white_gray = fieldsList[1]  	
        wgs_exclusion_trafic = fieldsList[2]   	  	
        dcc_project_code = fieldsList[5]
 
        tumor_wgs_aliquot_ids = fieldsList[11]
        
        # Discard excluded donors        
        if (wgs_exclusion_trafic != "Excluded"):

            # For each tumor sample from a given donor (>1 for multitumor donors)
            tumor_wgs_aliquot_ids_list = tumor_wgs_aliquot_ids.split(",")
        
            for tumor_wgs_aliquot_id in tumor_wgs_aliquot_ids_list:

                ## Create donor directory:
                outDonorDir = outDir + "/" + dcc_project_code + "/" + submitted_donor_id + "/" + tumor_wgs_aliquot_id
    
                try:
                    os.makedirs(outDonorDir)

                # Do not create directory if already exists
                except OSError as exc:
                    if exc.errno == errno.EEXIST and os.path.isdir(outDonorDir):  
                        pass
                    else:
                        raise
        
                ## Create text file with processed pseudogenes:
                pseudogenesPath = outDonorDir + "/" + tumor_wgs_aliquot_id + ".somatic.pseudogenes.tsv"
                outFile = open(pseudogenesPath, "w" )

                # Print header into the output file
                header = "#chromPlus" + "\t" + "begPlus" + "\t" + "endPlus" + "\t" + "nbReadsPlus" + "\t" + "classPlus" + "\t" + "readListPlus" + "\t" + "chromMinus" + "\t" + "begMinus" + "\t" + "endMinus" + "\t" + "nbReadsMinus" + "\t" + "classMinus" + "\t" + "readListMinus" + "\t" + "insertionType" + "\t" + "chromSource" + "\t" + "begSource" + "\t" + "endSource" + "\t" + "strandSource" + "\t" + "tdBeg" + "\t" + "tdEnd" + "\t" + "tdRnaLen" + "\t" + "tdLen" + "\t" + "psdGene" + "\t" + "chromExonA" + "\t" + "begExonA" + "\t" + "endExonA" + "\t" + "chromExonB" + "\t" + "begExonB" + "\t" + "endExonB" + "\n"

                outFile.write(header)

                # Print insertions for those tumor samples with at least one pseudogene insertions
                if (tumor_wgs_aliquot_id in allPseudogenesDict):
                
                    # For each processed pseudogene:
                    for pseudogeneDict in allPseudogenesDict[tumor_wgs_aliquot_id]:
                        row = pseudogeneDict["chromPlus"] + "\t" + pseudogeneDict["begPlus"] + "\t" + pseudogeneDict["endPlus"] + "\t" + pseudogeneDict["nbReadsPlus"] + "\t" + pseudogeneDict["classPlus"] + "\t" + pseudogeneDict["readListPlus"] + "\t" + pseudogeneDict["chromMinus"] + "\t" + pseudogeneDict["begMinus"] + "\t" + pseudogeneDict["endMinus"] + "\t" + pseudogeneDict["nbReadsMinus"] + "\t" + pseudogeneDict["classMinus"] + "\t" + pseudogeneDict["readListMinus"] + "\t" + pseudogeneDict["insertionType"] + "\t" + pseudogeneDict["chromSource"] + "\t" + pseudogeneDict["begSource"] + "\t" + pseudogeneDict["endSource"] + "\t" + pseudogeneDict["strandSource"] + "\t" + pseudogeneDict["tdBeg"] + "\t" + pseudogeneDict["tdEnd"] + "\t" + pseudogeneDict["tdRnaLen"] + "\t" + pseudogeneDict["tdLen"] + "\t" + pseudogeneDict["psdGene"] + "\t" + pseudogeneDict["chromExonA"] + "\t" + pseudogeneDict["begExonA"] + "\t" + pseudogeneDict["endExonA"] + "\t" + pseudogeneDict["chromExonB"] + "\t" + pseudogeneDict["begExonB"] + "\t" + pseudogeneDict["endExonB"] + "\n"

                        # Print pseudogene information into the output file
                        outFile.write(row)   

print "***** Finished! *****"
print

