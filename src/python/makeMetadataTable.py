#!/usr/bin/env python
#coding: utf-8

## Load modules/libraries
import sys
import argparse
import os

## Get user's input
parser = argparse.ArgumentParser(description= "Takes a set of tables containing different pieces of PCWAG donors metadata and produces a master metadata table with all the relevant info")
parser.add_argument('mainMetadata', help='Main file containing metadata information')
parser.add_argument('multitumor', help='File with representative tumor samples for multitumor donors')
parser.add_argument('histology', help='File with donor histology information')
parser.add_argument('excludedTumorTypes', help='File with the list of excluded tumor subtypes')
parser.add_argument('ancestry', help='File with donor ancestry information')
parser.add_argument('traficBlackList', help='File with trafic blacklisted samples')
parser.add_argument('clinical', help='File with donor clinical information (gender, age at diagnosis...)')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
mainMetadataPath = args.mainMetadata
multitumorPath = args.multitumor
histologyPath = args.histology
excludedTumorTypesPath = args.excludedTumorTypes
ancestryPath = args.ancestry
traficBlackListPath = args.traficBlackList
clinicalPath = args.clinical
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "mainMetadataPath: ", mainMetadataPath
print "multitumorPath: ", multitumorPath
print "histologyPath: ", histologyPath
print "excludedTumorTypesPath: ", excludedTumorTypesPath
print "ancestryPath: ", ancestryPath
print "traficBlackList: ", traficBlackListPath
print "clinical: ", clinicalPath
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


### 2) Read file with histology information. 
# Store histology information into a nested dictionary
# Key1 (donorUniqueId) -> value (dict)
# Key2:
#   - "abbreviation" -> value
#   - "tier1" -> value
#   - "tier2" -> value

# #donor_unique_id    	histology_abbreviation	histology_tier1	histology_tier2

histology = open(histologyPath, 'r')
histologyDict = {}

# Read file line by line
for line in histology:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        donorUniqueId = fieldsList[0] 
        histologyAbbrv = fieldsList[1] 
        histologyTier1 = fieldsList[2]
        histologyTier2 = fieldsList[3]

        histologyDict[donorUniqueId] = {}
        histologyDict[donorUniqueId]["abbreviation"] = histologyAbbrv
        histologyDict[donorUniqueId]["tier1"] = histologyTier1
        histologyDict[donorUniqueId]["tier2"] = histologyTier2

### 3) Read file with with the list of excluded tumor subtypes
excludedTumorTypes = open(excludedTumorTypesPath, 'r')
excludedTumorTypesList = []

# Read file line by line
for line in excludedTumorTypes:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")
        excludedTumorType = fieldsList[0]
        excludedTumorTypesList.append(excludedTumorType)


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
        icgc_donor_id = fieldsList[0] 
        blackList.append(icgc_donor_id)

### 5) Read file with donor clinical information. 
# Store clinical information into a nested dictionary
# Key1 (donorUniqueId) -> value (dict)
# Key2:
#   - "sex" -> value
#   - "diagnosisAge" -> value
#   - "survivalTime" -> value
#   - "intervalLastFollowup" -> value

clinical = open(clinicalPath, 'r')
clinicalDict = {}

# Read file line by line
for line in clinical:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        donorUniqueId = fieldsList[0] 
        clinicalDict[donorUniqueId] = {}
        clinicalDict[donorUniqueId]["sex"] = fieldsList[1]
        clinicalDict[donorUniqueId]["diagnosisAge"] = fieldsList[2]
        clinicalDict[donorUniqueId]["survivalTime"] = fieldsList[3]
        clinicalDict[donorUniqueId]["intervalLastFollowup"] = fieldsList[4]

### 6) Read main metadata file and make output file containing all the metadata information
# Open output file
fileName = "PCAWG_donors_metadata.tsv"
outFilePath = outDir + "/" + fileName
outFile = open( outFilePath, "w" )

# Write header:
row = "#submitted_donor_id" + "\t" + "icgc_donor_id" + "\t" + "wgs_exclusion_white_gray" + "\t" + "wgs_exclusion_trafic" + "\t" + "ancestry_primary" + "\t" + "donor_sex" + "\t" + "donor_age_at_diagnosis" + "\t" + "donor_survival_time" + "\t" + "donor_interval_of_last_followup" + "\t" + "dcc_project_code" + "\t" + "histology_count" + "\t" + "histology_exclusion_status" + "\t" + "histology_abbreviation"	 + "\t" + "histology_tier1" + "\t" +	 "histology_tier2" + "\t" + "normal_wgs_icgc_specimen_id" + "\t" + "normal_wgs_icgc_sample_id" + "\t" + "normal_wgs_aliquot_id" + "\t" + "tumor_wgs_specimen_count" + "\t" + "tumor_wgs_icgc_specimen_id" + "\t" + "tumor_wgs_icgc_sample_id" + "\t" + "tumor_wgs_aliquot_id" + "\t" + "tumor_wgs_representative_aliquot_id" + "\t" + "normal_wgs_has_matched_rna_seq" + "\t" + "tumor_wgs_has_matched_rna_seq" + "\t" + 	"normal_rna_seq_icgc_specimen_id" + "\t" + "normal_rna_seq_icgc_sample_id" + "\t" + "tumor_rna_seq_specimen_count" + "\t" + "tumor_rna_seq_icgc_specimen_id" + "\t" + "tumor_rna_seq_icgc_sample_id" + "\t" + "tumor_rna_seq_aliquot_id" + "\n" 

outFile.write(row)

mainMetadata = open(mainMetadataPath, 'r')

# Read file line by line
for line in mainMetadata:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):
        
        fieldsList = line.split("\t")

        submitted_donor_id = fieldsList[0]	
        icgc_donor_id = fieldsList[1]	
        wgs_exclusion_white_gray = fieldsList[2]	
        dcc_project_code = fieldsList[3]	
        normal_wgs_icgc_specimen_id = fieldsList[4]	
        normal_wgs_icgc_sample_id = fieldsList[5]	
        normal_wgs_aliquot_id = fieldsList[6]	
        tumor_wgs_specimen_count = fieldsList[7]	
        tumor_wgs_icgc_specimen_id = fieldsList[8]	
        tumor_wgs_icgc_sample_id = fieldsList[9]	
        tumor_wgs_aliquot_id = fieldsList[10]	
        normal_wgs_has_matched_rna_seq = fieldsList[11]	
        tumor_wgs_has_matched_rna_seq = fieldsList[12]	
        normal_rna_seq_icgc_specimen_id = fieldsList[13]	
        normal_rna_seq_icgc_sample_id = fieldsList[14]	
        tumor_rna_seq_specimen_count = fieldsList[15]	
        tumor_rna_seq_icgc_specimen_id = fieldsList[16]	
        tumor_rna_seq_icgc_sample_id = fieldsList[17]	
        tumor_rna_seq_aliquot_id = fieldsList[18]

        donorUniqueId = dcc_project_code + "::" + submitted_donor_id
        
        ## Histology        
        histology_abbreviation = histologyDict[donorUniqueId]["abbreviation"] 
        histology_count = len(histology_abbreviation.split(','))
        histology_tier1 = histologyDict[donorUniqueId]["tier1"] 
        histology_tier2 = histologyDict[donorUniqueId]["tier2"]   
        histology_exclusion_status = 'Excluded' if histology_abbreviation in excludedTumorTypesList else 'included'

        ## Clinical
        donor_sex = clinicalDict[donorUniqueId]["sex"] 
        donor_age_at_diagnosis = clinicalDict[donorUniqueId]["diagnosisAge"] 
        donor_survival_time = clinicalDict[donorUniqueId]["survivalTime"] 
        donor_interval_of_last_followup = clinicalDict[donorUniqueId]["intervalLastFollowup"] 

        # a) Donor with a single tumor sample
        if (tumor_wgs_specimen_count == "1"):
            tumor_wgs_representative_aliquot_id = "NA"
        
        # b) Multi tumor sample donor 
        else:
            tumor_wgs_representative_aliquot_id = multitumorDict[donorUniqueId]
    
        ## Ancestry
        # a) Donor without ancestry information
        # Note: there are 16 donors we do not have ancestry information...         
        if normal_wgs_aliquot_id not in ancestryDict:
            ancestry_primary = "UNK"  
    
        # b) Donor with ancestry information available         
        else:     
            ancestry_primary = ancestryDict[normal_wgs_aliquot_id]

        ## Blacklist
        # a) Donor in TraFiC blacklist
        if icgc_donor_id in blackList:
            wgs_exclusion_trafic = "Excluded"
            
        # b) Donor not in blacklist
        else:
            wgs_exclusion_trafic = "Whitelist"


        ### Write metadata row into the output file
        row = submitted_donor_id + "\t" + icgc_donor_id + "\t" + wgs_exclusion_white_gray + "\t" + wgs_exclusion_trafic + "\t" + ancestry_primary + "\t" + donor_sex + "\t" + donor_age_at_diagnosis + "\t" + donor_survival_time + "\t" + donor_interval_of_last_followup + "\t" + dcc_project_code + "\t" + str(histology_count) + "\t" + histology_exclusion_status + "\t" + histology_abbreviation	 + "\t" + histology_tier1 + "\t" +	 histology_tier2 + "\t" + normal_wgs_icgc_specimen_id + "\t" + normal_wgs_icgc_sample_id + "\t" + normal_wgs_aliquot_id + "\t" + tumor_wgs_specimen_count	 + "\t" + tumor_wgs_icgc_specimen_id + "\t" + tumor_wgs_icgc_sample_id + "\t" + tumor_wgs_aliquot_id + "\t" + tumor_wgs_representative_aliquot_id + "\t" + normal_wgs_has_matched_rna_seq + "\t" + tumor_wgs_has_matched_rna_seq + "\t" + normal_rna_seq_icgc_specimen_id + "\t" + normal_rna_seq_icgc_sample_id + "\t" + tumor_rna_seq_specimen_count + "\t" + tumor_rna_seq_icgc_specimen_id + "\t" + tumor_rna_seq_icgc_sample_id + "\t" + tumor_rna_seq_aliquot_id + "\n" 

        outFile.write(row)


print "***** Finished! *****"
print

