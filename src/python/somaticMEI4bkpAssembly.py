#!/usr/bin/env python
#coding: utf-8


### Functions ###
def overlap(begA, endA, begB, endB):

    """
    Check if both ranges overlap. 2 criteria for defining overlap: 
    ## A) Begin of the range A within the range B         
    #       *beg* <---------range_A---------->                         
    # <---------range_B----------> 
                
    #    *beg* <-------range_A----->
    # <-------------range_B------------------>

    ## B) Begin of the range B within the range A     
    # <---------range_A----------> 
    #               *beg* <---------range_B---------->
          
    # <-------------range_A----------------->
    #    *beg* <-------range_B------>

    """    
       
    # a) Begin of the range A within the range B   
    if ((begA >= begB) and (begA <= endB)):
        overlap = True
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True

    # c) Ranges do not overlapping
    else:
        overlap = False

    return overlap


## Load modules/libraries
import sys
import argparse
import os
import errno

## Get user's input
parser = argparse.ArgumentParser(description= "Produce a correctly formated file for executing breakpoint assembly pipeline on the somatic insertions from a given donor")
parser.add_argument('traficOut', help='Main TraFiC output file')
parser.add_argument('clusterPaths', help='Path to td0, td1 and td2 TraFiC cluster files')
parser.add_argument('sourceMetadataPath', help='Path to source elements metadata file')
parser.add_argument('outFileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory.' )

args = parser.parse_args()
traficOut = args.traficOut
clusterPaths = args.clusterPaths
sourceMetadataPath = args.sourceMetadataPath
outFileName = args.outFileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])


## Display configuration to standard output
print
print "***** ", scriptName, " configuration *****"
print "traficOut: ", traficOut
print "clusters: ", clusterPaths
print "sourceMetadataPath: ", sourceMetadataPath
print "outFileName: ", outFileName
print "outDir: ", outDir
print

print "***** Executing ", scriptName, " *****"
print
print "..."
print


### 1) Read main TraFiC output file and gather transductions info in a dictionary
traficOut = open(traficOut, 'r')
tdInfoDict = {}

# Read file line by line
for line in traficOut:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):

        # Extract needed info
        fieldsList = line.split("\t")
        
        chrom = fieldsList[1]
        beg = fieldsList[2]
        end = fieldsList[3]
        insertionClass = fieldsList[4]
        insertionId = insertionClass + '_' + chrom + "_" + beg + "_" + end
        tdInfoDict[insertionId] = {}
        
        insertionType = fieldsList[6]
        chromSource = fieldsList[7]
        begSource = fieldsList[8] 
        endSource = fieldsList[9]
        strandSource = fieldsList[10]
        tdBeg = fieldsList[11]
        tdEnd = fieldsList[12]
        tdRnaLen = fieldsList[13]
        tdLen = fieldsList[14]

        ## Add info to the nested dictionary:
        tdInfoDict[insertionId]["insertionType"] = insertionType
        tdInfoDict[insertionId]["chromSource"] = chromSource
        tdInfoDict[insertionId]["begSource"] = begSource
        tdInfoDict[insertionId]["endSource"] = endSource
        tdInfoDict[insertionId]["strandSource"] = strandSource
        tdInfoDict[insertionId]["tdBeg"] = tdBeg
        tdInfoDict[insertionId]["tdEnd"] = tdEnd
        tdInfoDict[insertionId]["tdRnaLen"] = tdRnaLen
        tdInfoDict[insertionId]["tdLen"] = tdLen


### 2) Read source elements metadata dictionary and generate a dictionary with the following structure:
# dict1 -> Key(chrom) -> tuple(cytobandId, sourceCoord)

sourceMetadata = open(sourceMetadataPath, 'r')
sourceMetadataDict = {}

# Read file line by line
for line in sourceMetadata:
    line = line.rstrip('\r\n')

    ## Discard header
    if not line.startswith("#"):

        # Extract needed info
        fieldsList = line.split("\t")
        cytobandId =	 fieldsList[0]
        chrom = fieldsList[1]
        beg = fieldsList[2]
        end = fieldsList[3]
        sourceCoord = chrom + ':' + beg + '-' + end
        
        if chrom not in sourceMetadataDict:
            sourceMetadataDict[chrom] = []

        sourceTuple = (cytobandId, sourceCoord)
        sourceMetadataDict[chrom].append(sourceTuple)


### 3) Read TraFiC clusters output files and generate a single properly formated file for breakpoint assembly. 
## Create output file with insertions:
outPath = outDir + "/" + outFileName + ".TraFiC.tsv"
outFile = open(outPath, "w" )

# Print header into the output file
header = "#chromPlus" + "\t" + "begPlus" + "\t" + "endPlus" + "\t" + "nbReadsPlus" + "\t" + "classPlus" + "\t" + "readListPlus" + "\t" + "chromMinus" + "\t" + "begMinus" + "\t" + "endMinus" + "\t" + "nbReadsMinus" + "\t" + "classMinus" + "\t" + "readListMinus" + "\t" + "insertionType" + "\t" + "cytobandId" + "\t" + "sourceType" + "\t" + "chromSource" + "\t" + "begSource" + "\t" + "endSource" + "\t" + "strandSource" + "\t" + "tdBeg" + "\t" + "tdEnd" + "\t" + "tdRnaLen" + "\t" + "tdLen" + "\t" + "psdGene" + "\t" + "chromExonA" + "\t" + "begExonA" + "\t" + "endExonA" + "\t" + "chromExonB" + "\t" + "begExonB" + "\t" + "endExonB" + "\t" + "grType" + "\n"

outFile.write(header) 

clusterPaths = open(clusterPaths, 'r')
insertionList = []

# Read file line by line
for line in clusterPaths:
    line = line.rstrip('\r\n')
 
    fieldsList = line.split("\t")
    insertionType = fieldsList[0]
    clusterFile = fieldsList[1]
    clusterFile = open(clusterFile, 'r')
        
    for cluster in clusterFile:

        ## Discard header
        if not cluster.startswith("#"):

            # Extract needed info
            fieldsList = cluster.split("\t")
        
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

            insertionId = classPlus + '_' + chromPlus + "_" + str(int(endPlus)+100) + "_" + begMinus           

            ## Select only insertions into TraFiC final output file and discard repeated insertions
            if (insertionId in tdInfoDict) and (insertionId not in insertionList):
                                                          
                insertionList.append(insertionId)  
                insertionType = tdInfoDict[insertionId]["insertionType"]
                
                # A) Solo insertion
                if (insertionType == "td0"):
                    
                    insertionType = "TD0"
                    cytobandId = "NA" 
                    sourceType = "NA"
                    chromSource = "NA" 
                    begSource = "NA"
                    endSource = "NA"
                    strandSource = "NA"
                    tdBeg = "NA"
                    tdEnd = "NA"
                    tdRnaLen = "NA"
                    tdLen = "NA"
                
                # B) Partnered transduction
                elif (insertionType == "td1"):
                                                                                                   
                    insertionType = "TD1" 
                    chromSource = tdInfoDict[insertionId]["chromSource"] 
                    begSource = tdInfoDict[insertionId]["begSource"] 
                    endSource = tdInfoDict[insertionId]["endSource"] 
                    strandSource = tdInfoDict[insertionId]["strandSource"] 
                             
                    #### Determine if germline source element gave rise to the transduction:
                    ## By default consider the source element as somatic
                    cytobandId = "NA"
                    sourceType = "SOMATIC"
        
                    ### Assess if it overlaps with a germline source element -> then germline 
                    # There are germline source elements in the chromosome
                    if (chromSource in sourceMetadataDict):
        
                        # For germline source element 
                        for sourceTuple in sourceMetadataDict[chromSource]:
        
                            ## Define transduction source element coordinates range
                            begA = int(begSource) - 1000
                            endA = int(endSource) + 1000
        
                            ## Define germline source element range
                            coords = sourceTuple[1].split(":")[1]
                            begB, endB = coords.split("-")
                            
                            begB = int(endB)-6500 if (begB == "UNK") else int(begB)
                            endB = int(begB)+6500 if (endB == "UNK") else int(endB)
        
                            ## Search for overlap                    
                            if overlap(int(begA), int(endA), int(begB), int(endB)):
                                cytobandId = sourceTuple[0]
                                sourceType = "GERMLINE"
        
                                # Update source element beg and end coordinates
                                sourceCoord = sourceTuple[1]
                                chrom, coords = sourceCoord.split(":")
                                begSource, endSource = coords.split("-")
                   
                                break
                                                                    
                    ## Determine transduction beginning and end coordinates depending on source element orientation
                    # a) Source element in forward
                    if (strandSource == "plus"):
                     
                        strandSource = "+"
                        tdBeg = endSource
                        tdEnd = tdInfoDict[insertionId]["tdEnd"]
                        tdRnaLen = str(int(tdEnd) - int(tdBeg))
                        tdLen = tdRnaLen     
        
                    # b) Source element in reverse
                    elif (strandSource == "minus"):

                        strandSource = "-"
                        tdBeg = tdInfoDict[insertionId]["tdBeg"]
                        tdEnd = begSource
                        tdRnaLen = str(int(tdEnd) - int(tdBeg))
                        tdLen = tdRnaLen
        
                    # c) Unknown orientation  
                    else:

                        tdBeg = tdInfoDict[insertionId]["tdBeg"]
                        tdEnd = tdInfoDict[insertionId]["tdEnd"]
                        tdRnaLen = "NA"
                        tdLen = "NA"
                        
                # C) Orphan transduction:
                else:                
        
                    insertionType = "TD2"                   
     
                    chromSource = tdInfoDict[insertionId]["chromSource"]
                    begSource = tdInfoDict[insertionId]["begSource"]
                    endSource = tdInfoDict[insertionId]["endSource"]
                    strandSource = tdInfoDict[insertionId]["strandSource"]
                    tdBeg = tdInfoDict[insertionId]["tdBeg"]
                    tdEnd = tdInfoDict[insertionId]["tdEnd"]
        
                    #### Determine if germline source element gave rise to the transduction:
                    ## By default consider the source element as somatic
                    cytobandId = "NA"
                    sourceType = "SOMATIC"
        
                    ### Assess if it overlaps with a germline source element -> then germline 
                    # Ther are germline source elements in the chromosome
                    if (chromSource in sourceMetadataDict):
        
                        # For germline source element 
                        for sourceTuple in sourceMetadataDict[chromSource]:
        
                            ## Define transduction source element coordinates range
                            begA = int(begSource) - 1000
                            endA = int(endSource) + 1000
        
                            ## Define germline source element range
                            coords = sourceTuple[1].split(":")[1]
                            begB, endB = coords.split("-")
                            
                            begB = int(endB)-6500 if (begB == "UNK") else int(begB)
                            endB = int(begB)+6500 if (endB == "UNK") else int(endB)
        
                            ## Search for overlap                    
                            if overlap(int(begA), int(endA), int(begB), int(endB)):
                                cytobandId = sourceTuple[0]
                                sourceType = "GERMLINE"
        
                                # Update source element beg and end coordinates
                                sourceCoord = sourceTuple[1]
                                chrom, coords = sourceCoord.split(":")
                                begSource, endSource = coords.split("-")
                   
                                break 
        
                    # a) Source element in forward
                    if (strandSource == "plus"):
        
                        strandSource = "+"
                        tdRnaLen = str(abs(int(tdInfoDict[insertionId]["tdRnaLen"])))
                        tdLen = str(abs(int(tdInfoDict[insertionId]["tdLen"]))) 

                    # b) Source element in reverse
                    elif (strandSource == "minus"):
                
                        strandSource = "-"
                        tdRnaLen = str(abs(int(tdInfoDict[insertionId]["tdRnaLen"])))
                        tdLen = str(abs(int(tdInfoDict[insertionId]["tdLen"]))) 
                
                    # c) Unknown orientation
                    else:
                
                        tdRnaLen = "NA"
                        tdLen = str(abs(int(tdInfoDict[insertionId]["tdLen"]))) 
                                        
                ## Set as NA pseudogene specific fields:
                psdGene = "NA"
                chromExonA = "NA"
                begExonA = "NA"
                endExonA = "NA"
                chromExonB = "NA"
                begExonB = "NA"
                endExonB = "NA"

                ## Set as NA rearrangement specific fields:
                grType = "NA"                

                ## Write properly formated line into the output file:            
                row = chromPlus + "\t" + begPlus + "\t" + endPlus + "\t" + nbReadsPlus + "\t" + classPlus + "\t" + readListPlus + "\t" + chromMinus + "\t" + begMinus + "\t" + endMinus + "\t" + nbReadsMinus + "\t" + classMinus + "\t" + readListMinus + "\t" + insertionType + "\t" + cytobandId + "\t" + sourceType + "\t" + chromSource + "\t" + begSource + "\t" + endSource + "\t" + strandSource + "\t" + tdBeg + "\t" + tdEnd + "\t" + tdRnaLen + "\t" + tdLen + "\t" + psdGene + "\t" + chromExonA + "\t" + begExonA + "\t" + endExonA + "\t" + chromExonB + "\t" + begExonB + "\t" + endExonB + "\t" + grType + "\n"        

                outFile.write(row) 


print "***** Finished! *****"
print

