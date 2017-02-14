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
# Key1 [tumor_wgs_aliquot_id] -> value [list of insertion dictionaries: insertionDict1, insertionDict2, ...]
# Note: The number of elements in the list will be equal to the total number of insertions in a given tumor sample. 
# Note: The insertion dictionary contains all the info regarding a given insertion

insertions = open(insertionsPath, 'r')
allInsertionsDict = {}

# Read file line by line
for line in insertions:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
        
        chromPlus = fieldsList[0]
        begPlus = fieldsList[1] 
        endPlus = fieldsList[2]  
        nbReadsPlus = fieldsList[3]
        classPlus = fieldsList[4]
        readListPlus = fieldsList[5]
        chromMinus = fieldsList[6]
        begMinus = fieldsList[7] 
        endMinus = fieldsList[8]   
        nbReadsMinus = fieldsList[9] 
        classMinus = fieldsList[10]
        readListMinus = fieldsList[11]

        tumor_wgs_aliquot_id = fieldsList[15]
            
        insertionType = fieldsList[21] 
        
        # A) Solo insertion
        if (insertionType == "td0"):

            insertionType = "TD0" 
            chromSource = "NA" 
            begSource = "NA"
            endSource = "NA"
            strandSource = "NA"
            tdBeg = "NA"
            tdEnd = "NA"
            tdRnaLen = "NA"
            tdLen = "NA"

        # B) Partnered transduction
        elif (insertionType == "td1"):
                
            insertionType = "TD1" 
            chromSource = fieldsList[22]
            begSource = fieldsList[23]
            endSource = fieldsList[24]
            strandSource = fieldsList[25]

            # a) Source element in forward
            if (strandSource == "plus"):

                strandSource = "+"
                tdBeg = endSource
                tdEnd = fieldsList[27]
                tdRnaLen = fieldsList[28]
                tdLen = fieldsList[29]        

            # b) Source element in reverse
            elif (strandSource == "minus"):

                strandSource = "-"
                tdBeg = fieldsList[26]
                tdEnd = begSource
                tdRnaLen = fieldsList[28]
                tdLen = fieldsList[29]

            # c) Unknown orientation  
            else:
            
                tdBeg = fieldsList[26]
                tdEnd = fieldsList[27]
                tdRnaLen = "NA"
                tdLen = "NA"

        # C) Orphan transduction:
        else:                

            insertionType = "TD2"                    
            chromSource = fieldsList[22]
            begSource = fieldsList[23]
            endSource = fieldsList[24]
            strandSource = fieldsList[25]
            tdBeg = fieldsList[26]
            tdEnd = fieldsList[27]

            # a) Source element in forward
            if (strandSource == "plus"):

                strandSource = "+"
                tdRnaLen = fieldsList[28]
                tdLen = fieldsList[29]

            # b) Source element in reverse
            elif (strandSource == "minus"):
        
                strandSource = "-"
                tdRnaLen = fieldsList[28]
                tdLen = fieldsList[29]
        
            # c) Unknown orientation
            else:
        
                tdRnaLen = "NA"
                tdLen = "NA" 

        insertionDict = {}

        insertionDict["chromPlus"] = chromPlus 
        insertionDict["begPlus"] = begPlus  
        insertionDict["endPlus"] = endPlus   
        insertionDict["nbReadsPlus"] = nbReadsPlus 
        insertionDict["classPlus"] = classPlus 
        insertionDict["readListPlus"] = readListPlus 
        insertionDict["chromMinus"] = chromMinus 
        insertionDict["begMinus"] = begMinus  
        insertionDict["endMinus"] = endMinus    
        insertionDict["nbReadsMinus"] = nbReadsMinus  
        insertionDict["classMinus"] = classMinus 
        insertionDict["readListMinus"] = readListMinus 
        insertionDict["insertionType"] = insertionType 
        insertionDict["chromSource"] = chromSource  
        insertionDict["begSource"] = begSource 
        insertionDict["endSource"] = endSource 
        insertionDict["strandSource"] = strandSource 
        insertionDict["tdBeg"] = tdBeg 
        insertionDict["tdEnd"] = tdEnd 
        insertionDict["tdRnaLen"] = tdRnaLen 
        insertionDict["tdLen"] = tdLen
        insertionDict["psdGene"] = "NA"
        insertionDict["chromExonA"] = "NA"
        insertionDict["begExonA"] = "NA"
        insertionDict["endExonA"] = "NA"
        insertionDict["chromExonB"] = "NA"
        insertionDict["begExonB"] = "NA"
        insertionDict["endExonB"] = "NA"

        if tumor_wgs_aliquot_id not in allInsertionsDict:
            allInsertionsDict[tumor_wgs_aliquot_id] = []

        allInsertionsDict[tumor_wgs_aliquot_id].append(insertionDict)
        

### 2) For each PCAWG tumor sample generate a properly formated file containing its somatic MEI and
# organize these files into the following directory architecture:
# outDir -> dcc_project_code > submitted_donor_id -> tumor_wgs_aliquot_id -> **somaticMEI tsv**

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
        
        # Discard TraFiC excluded donors         
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
        
                ## Create text file with insertions:
                insertionsPath = outDonorDir + "/" + tumor_wgs_aliquot_id + ".somatic.td0_td1_td2.tsv"
                outFile = open(insertionsPath, "w" )

                # Print header into the output file
                header = "#chromPlus" + "\t" + "begPlus" + "\t" + "endPlus" + "\t" + "nbReadsPlus" + "\t" + "classPlus" + "\t" + "readListPlus" + "\t" + "chromMinus" + "\t" + "begMinus" + "\t" + "endMinus" + "\t" + "nbReadsMinus" + "\t" + "classMinus" + "\t" + "readListMinus" + "\t" + "insertionType" + "\t" + "chromSource" + "\t" + "begSource" + "\t" + "endSource" + "\t" + "strandSource" + "\t" + "tdBeg" + "\t" + "tdEnd" + "\t" + "tdRnaLen" + "\t" + "tdLen" + "\t" + "psdGene" + "\t" + "chromExonA" + "\t" + "begExonA" + "\t" + "endExonA" + "\t" + "chromExonB" + "\t" + "begExonB" + "\t" + "endExonB" + "\n"
                outFile.write(header)

                # Print insertions for those tumor samples with at least one insertion
                if (tumor_wgs_aliquot_id in allInsertionsDict):
                
                    # For each insertion:
                    for insertionDict in allInsertionsDict[tumor_wgs_aliquot_id]:
                        row = insertionDict["chromPlus"] + "\t" + insertionDict["begPlus"] + "\t" + insertionDict["endPlus"] + "\t" + insertionDict["nbReadsPlus"] + "\t" + insertionDict["classPlus"] + "\t" + insertionDict["readListPlus"] + "\t" + insertionDict["chromMinus"] + "\t" + insertionDict["begMinus"] + "\t" + insertionDict["endMinus"] + "\t" + insertionDict["nbReadsMinus"] + "\t" + insertionDict["classMinus"] + "\t" + insertionDict["readListMinus"] + "\t" + insertionDict["insertionType"] + "\t" + insertionDict["chromSource"] + "\t" + insertionDict["begSource"] + "\t" + insertionDict["endSource"] + "\t" + insertionDict["strandSource"] + "\t" + insertionDict["tdBeg"] + "\t" + insertionDict["tdEnd"] + "\t" + insertionDict["tdRnaLen"] + "\t" + insertionDict["tdLen"] + "\t" + insertionDict["psdGene"] + "\t" + insertionDict["chromExonA"] + "\t" + insertionDict["begExonA"] + "\t" + insertionDict["endExonA"] + "\t" + insertionDict["chromExonB"] + "\t" + insertionDict["begExonB"] + "\t" + insertionDict["endExonB"] + "\n"

                        # Print insertion information into the output file
                        outFile.write(row)   

print "***** Finished! *****"
print

