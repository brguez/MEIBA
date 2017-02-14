#!/usr/bin/env python
#coding: utf-8

## Load modules/libraries
import sys
import argparse
import os

## Get user's input
parser = argparse.ArgumentParser(description= """""")
parser.add_argument('mainMetadata', help='Main file containing metadata information')
parser.add_argument('multitumor', help='File with representative tumor samples for multitumor donors')
parser.add_argument('histology', help='File with donor histology information')
parser.add_argument('ancestry', help='File with donor ancestry information')
parser.add_argument('traficBlackList', help='File with trafic blacklisted samples')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
mainMetadataPath = args.mainMetadata
multitumorPath = args.multitumor
histologyPath = args.histology
ancestryPath = args.ancestry
traficBlackListPath = args.traficBlackList
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "mainMetadataPath: ", mainMetadataPath
print "multitumorPath: ", multitumorPath
print "histologyPath: ", histologyPath
print "ancestryPath: ", ancestryPath
print "traficBlackList: ", traficBlackListPath
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print


### 1) Read file with representative tumor_aliquot_ids for multitumor donors
# Make dictionary with donorUniqueId (projectCode::donorId) as key and representative tumor_aliquot_ids as values

multitumor = open(multitumorPath, 'r')
multitumorDict = {}

# Read file line by line
for line in multitumor:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
        donorUniqueId = fieldsList[0]
        tumorAliquotId = fieldsList[1]
        multitumorDict[donorUniqueId] = tumorAliquotId


### 2) Read file with specimenIds and histology information. 
# Store specimen and histology information into a nested dictionary
# Key1 (donorUniqueId) -> value (dict)
# Key2:
#   - "specimenId" -> value
#   - "abbreviation" -> value
#   - "tier1" -> value
#   - "tier2" -> value

histology = open(histologyPath, 'r')
histologyDict = {}

# Read file line by line
for line in histology:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        donorUniqueId = fieldsList[0] 
        histologyDict[donorUniqueId] = {}
        histologyDict[donorUniqueId]["specimenId"] = fieldsList[1]
        histologyDict[donorUniqueId]["abbreviation"] = fieldsList[2]
        histologyDict[donorUniqueId]["tier1"] = fieldsList[3]
        histologyDict[donorUniqueId]["tier2"] = fieldsList[4]
 

### 3) Read file with representative donor's ancestry 
# Make dictionary with normalAliquotId as keys and donor's ancestry as values

ancestry = open(ancestryPath, 'r')
ancestryDict = {}

# Read file line by line
for line in ancestry:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
    
        normalAliquotId = fieldsList[0]
        ancestry = fieldsList[6]
        ancestryDict[normalAliquotId] = ancestry


### 4) Read TraFiC blacklist of excluded samples
traficBlackList = open(traficBlackListPath, 'r')
blackList = []

# Read file line by line
for line in traficBlackList:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
            
        fieldsList = line.split("\t")
        donorUniqueId = fieldsList[0] 
        blackList.append(donorUniqueId)


### 5) Read main metadata file and make output file containing all the metadata information
# Open output file
fileName = "PCAWG_donors_metadata.tsv"
outFilePath = outDir + "/" + fileName
outFile = open( outFilePath, "w" )

# Write header:
row = "#submitted_donor_id" + "\t" + "wgs_exclusion_white_gray" + "\t" + "wgs_exclusion_trafic" + "\t" + "icgc_specimen_id" + "\t" + "ancestry_primary" + "\t" + "dcc_project_code" + "\t" + "histology_abbreviation"	 + "\t" + "histology_tier1" + "\t" +	 "histology_tier2" + "\t" + "normal_wgs_aliquot_id" + "\t" +	"tumor_wgs_specimen_count" + "\t" +	"tumor_wgs_aliquot_id" + "\t" + "tumor_wgs_representative_aliquot_id" + "\t" + "normal_wgs_has_matched_rna_seq"	 + "\t" + "tumor_wgs_has_matched_rna_seq" + "\n"

outFile.write(row)


mainMetadata = open(mainMetadataPath, 'r')

# Read file line by line
for line in mainMetadata:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        donorId = fieldsList[0]	
        donorStatus = fieldsList[1]	    
        projectCode = fieldsList[2]
        normalAliquotId	= fieldsList[3]
        tumorSpecimenCount = fieldsList[4]     	
        tumorAliquotId = fieldsList[5]  	
        normalMatchedRnaSeq = fieldsList[6]
        tumorMatchedRnaSeq = fieldsList[7] 

        donorUniqueId = projectCode + "::" + donorId
        specimenId = histologyDict[donorUniqueId]["specimenId"] 
        abbreviation = histologyDict[donorUniqueId]["abbreviation"] 
        tier1 = histologyDict[donorUniqueId]["tier1"] 
        tier2 = histologyDict[donorUniqueId]["tier2"]    

        # a) Donor with a single tumor sample
        if (tumorSpecimenCount == "1"):
            tumorRepresentativeAliquotId = "NA"
        
        # b) Multi tumor sample donor 
        else:
            tumorRepresentativeAliquotId = multitumorDict[donorUniqueId]
    
        # a) Donor without ancestry information
        # Note: there are 16 donors we do not have ancestry information...         
        if normalAliquotId not in ancestryDict:
            ancestry = "UNK"  
    
        # b) Donor with ancestry information available         
        else:     
            ancestry = ancestryDict[normalAliquotId]

        # a) Donor in TraFiC blacklist
        if donorUniqueId in blackList:
            donorStatus2 = "Excluded"
            
        # b) Donor not in blacklist
        else:
            donorStatus2 = "Whitelist"


        # Write metadata row into the output file
        row = donorId + "\t" + donorStatus + "\t" + donorStatus2 + "\t" + specimenId + "\t" + ancestry + "\t" + projectCode + "\t" + abbreviation + "\t" + tier1 + "\t" +	tier2 + "\t" + normalAliquotId + "\t" +	tumorSpecimenCount + "\t" +	tumorAliquotId + "\t" + tumorRepresentativeAliquotId + "\t" + normalMatchedRnaSeq + "\t" + tumorMatchedRnaSeq + "\n" 

        outFile.write(row)


print "***** Finished! *****"
print

