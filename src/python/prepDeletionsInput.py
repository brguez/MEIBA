#!/usr/bin/env python
#coding: utf-8

# Description:
#   Parse TraFiC output, extract TD* event meta data from "+"- and "-" cluster

#   CSV files.
# Input:
#   1) eventsFile:  TAB-separated file containing at least 7 fields:
#                   cancer_id	donor_id	tumor_sample	chr	start	end	length
#   2) samplesPath: Directory containing cluster meta info for samples.
#                   Expected structure:
#                     data/L1-DEL/samples/
#                     └── ESAD-UK
#                         ├── OCCAMS-AH-011
#                         │   ├── dd7d623b-b9af-4147-9aa6-e09793691f10.clusters_mas.filteredNormalcihg19.txt
#                         │   ├── dd7d623b-b9af-4147-9aa6-e09793691f10.clusters_menos.filteredNormalcihg19.txt
#                         │   ├── dd7d623b-b9af-4147-9aa6-e09793691f10.indepclusters.deftd2.DEL.mas.ok.txt
#                         │   ├── dd7d623b-b9af-4147-9aa6-e09793691f10.indepclusters.deftd2.DEL.menos.ok.txt
#                         [...]
#  3) outDir: Directory to write output to.
#             The following directory structure will be created:
#               outputDir/
#               ├── ESAD-UK
#               │   ├── OCCAMS-AH-011
#               │   │   └── dd7d623b-b9af-4147-9aa6-e09793691f10
#               │   │       └── L1-del.txt
#               [...]


from __future__ import print_function
import argparse
import os
import sys

# input files are expexcted to follow this naming convention:
input_fn_norm_p = "%s.clusters_mas.filteredNormalcihg19.txt"   # insert: sample_id
input_fn_norm_m = "%s.clusters_menos.filteredNormalcihg19.txt" # insert: sample_id
input_fn_del_p  = "%s.indepclusters.deftd2.DEL.mas.ok.txt"     # insert: sample_id
input_fn_del_m  = "%s.indepclusters.deftd2.DEL.menos.ok.txt"   # insert: sample_id

# output lines follow this scheme:
# see description in main pipeline script "TEIBA.sh" for documentation of fields
out_line_fmt = "%s" + ("\t%s" * 12) + ("\tNA" * 17)  + "\t" + "DEL"

### Functions ###
def makeSourceElementDict(sourceMetadataPath):
    """
    Read source elements metadata file and generate a dictionary with the following structure:
        dict1 -> Key(chrom) -> tuple(cytobandId, sourceCoord)
    """
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

    return sourceMetadataDict


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


def isGermlineSource(tdCoords, sourceMetadataDict):
    """
    Determine if germline source element gave rise to the transduction
    """

    chromSource, coords = tdCoords.split(":")
    begSource, endSource = coords.split("-")

    ## By default consider the source element as somatic
    cytobandId = "NA"
    sourceType = "SOMATIC"
        
    ### Assess if it overlaps with a germline source element -> then germline 
    # There are germline source elements in the chromosome
    if (chromSource in sourceMetadataDict):
        
        # For germline source element 
        for sourceTuple in sourceMetadataDict[chromSource]:
       
            ## Define transduction source element coordinates range
            begA = int(begSource) - 5000
            endA = int(endSource) + 5000
        
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

    return(cytobandId, sourceType, begSource, endSource)


def computeTdRnaLen(begSource, endSource, strandSource, tdBeg, tdEnd):
    """
    Determine transduction length at RNA level (to be tested)
    """
    ## Determine transduction RNA length depending on source element orientation
    # a) Source element in forward
    if (strandSource == "plus"):
                     
        strandSource = "+"
        tdBeg = endSource
        tdEnd = tdEnd
        tdRnaLen = str(int(tdEnd) - int(tdBeg))
        
    # b) Source element in reverse
    elif (strandSource == "minus"):

        strandSource = "-"
        tdBeg = tdBeg
        tdEnd = begSource
        tdRnaLen = str(int(tdEnd) - int(tdBeg))
        
    # c) Unknown orientation  
    else:

        tdBeg = tdBeg
        tdEnd = tdEnd
        tdRnaLen = "NA"

    return (strandSource, tdBeg, tdEnd, tdRnaLen)


def log(msg, event_type="INFO"):
    outfile = sys.stdout
    if event_type == "ERROR":
        outfile = sys.stderr
    print("[%s] (prepDeletionsInput.py) %s" % (event_type, msg), file=outfile)

def parse_args():
    parser = argparse.ArgumentParser(description="""Generates TEIBA input files for deletions from TraFiC output.""")
    parser.add_argument('eventsFile', help='TraFiC output file containing the coordinates of deletion events (one file for all samples).')
    parser.add_argument('samplesPath', help='Directory that contains sample folders.')
    parser.add_argument('sourceMetadataPath', help='L1 source elements metadata')
    parser.add_argument('outputDir', help='Output directory (will be created if not exists)')

    args = parser.parse_args()
    if not os.path.isfile(args.eventsFile):
        log("Input file does not exist: '%s'" % args.eventsFile, "ERROR")
        sys.exit(1)
    if not os.path.isdir(args.samplesPath):
        log("Input directory does not exist: '%s'" % args.samplesPath, "ERROR")
        sys.exit(1)
    if not os.path.isdir(args.outputDir):
        try:
            os.makedirs(args.outputDir)
        except:
            log("Output directory could not be created: '%s'" % args.outputDir, "ERROR")
            sys.exit(1)

    return args

if __name__ == "__main__":
    args = parse_args()

    scriptName = os.path.basename(sys.argv[0])

    ## Display configuration to standard output ##
    log("Events file: '%s'" % args.eventsFile, "INPUT")
    log("Sample path: '%s'" % args.samplesPath, "INPUT")
    log("Source metadata: '%s'" % args.sourceMetadataPath, "INPUT")
    log("Output directory: '%s'" % args.outputDir, "INPUT")

    # 0) Generate dictionary with source elements metadata
    sourceMetadataDict = makeSourceElementDict(args.sourceMetadataPath)

    # 1) collect event meta data for all samples
    num_samples = 0
    num_events = 0
    samples_events = {}
    for line in open(args.eventsFile, 'rt'):
        if line.startswith('#'):
            continue
        row = line.strip().split('\t')
        if len(row) >= 7:
          study, donor, sample, chrom, start, end = row[:6]
          id_sample = (study, donor, sample)
          if id_sample not in samples_events:
              samples_events[id_sample] = []
              num_samples += 1
          samples_events[id_sample] += [(chrom, start, end)]
          num_events += 1

    log("Found %d deletion events for %d samples" % (num_events, num_samples))

    # 2) parse sample files, outputting TEIBA input info for previously registered events
    num_events_missed = 0
    for study, donor, sample in sorted(samples_events):
        #log("%s_%s_%s" % (study, donor, sample))
        # create output folder
        outdir = os.path.join(args.outputDir, study, donor, sample)
        try:
            os.makedirs(outdir)
        except:
            log("Could not create output directory: '%s'" % outdir, "ERROR")

        # make sure sample files exist
        fn_norm_p = os.path.join(args.samplesPath, study, donor, input_fn_norm_p % sample)
        fn_norm_m = os.path.join(args.samplesPath, study, donor, input_fn_norm_m % sample)
        fn_del_p = os.path.join(args.samplesPath, study, donor, input_fn_del_p % sample)
        fn_del_m = os.path.join(args.samplesPath, study, donor, input_fn_del_m % sample)
        if not os.path.isfile(fn_norm_p):
            log("Missing input file: '%s'" % fn_norm_p, "WARN")
            continue
        if not os.path.isfile(fn_norm_m):
            log("Missing input file: '%s'" % fn_norm_m, "WARN")
            continue
        if not os.path.isfile(fn_del_p):
            log("Missing input file: '%s'" % fn_del_p, "WARN")
            continue
        if not os.path.isfile(fn_del_m):
            log("Missing input file: '%s'" % fn_del_m, "WARN")
            continue

        ## Parse sample input files
        ############################       
        #   columns: chrom, beg, end, num, class, reads, ...
        #print("## %s, %s, %s ##" % (study, donor, sample))

        ### A) Read positive normal clusters
        clust_norm_p_end = {}
        for line in open(fn_norm_p, 'rt'):
            #print("+#%s#" % line)
            row = line.strip().split()  
            chrom, beg, end, num, klass, reads = row[:6]
            clust_norm_p_end[(chrom, end)] = (chrom, beg, end, num, klass, reads)
        
        ### B) Read minus normal clusters
        clust_norm_m_beg = {}
        for line in open(fn_norm_m, 'rt'):
            #print("-#%s#" % line)
            row = line.strip().split()

            chrom, beg, end, num, klass, reads = row[:6]
            clust_norm_m_beg[(chrom, beg)] = (chrom, beg, end, num, klass, reads)
        
        ### C) Read positive transduction clusters
        clust_del_p_end = {}
        for line in open(fn_del_p, 'rt'):
            #print("+#%s#" % line)
            row = line.strip().split()

            if len(row) == 14:    
                chrom, beg, end, num, klass, reads, chromSource, begSource, endSource, strandSource, tdChrom, tdBeg, tdEnd = row[:13]
                clust_del_p_end[(chrom, end)] = (chrom, beg, end, num, klass, reads, chromSource, begSource, endSource, strandSource, tdChrom, tdBeg, tdEnd)
        
        ### D) Read positive transduction clusters
        clust_del_m_beg = {}
        for line in open(fn_del_m, 'rt'):
            #print("-#%s#" % line)
            row = line.strip().split()

            if len(row) == 14:    
                chrom, beg, end, num, klass, reads, chromSource, begSource, endSource, strandSource, tdChrom, tdBeg, tdEnd = row[:13]
                clust_del_m_beg[(chrom, beg)] = (chrom, beg, end, num, klass, reads, chromSource, begSource, endSource, strandSource, tdChrom, tdBeg, tdEnd)


        # combine cluster info for registered deletions
        #################################################
        id_sample = (study, donor, sample)
        with open(os.path.join(outdir, "L1-del.txt"), 'wt') as outfile:
            for chrom, beg, end in samples_events[id_sample]:
                insertionType = "NA"

                #### determine transduction type (TD0: solo, TD1: partnered, TD2: orphan)
                ## A) Solo insertion. Normal positive and negative cluster
                if ((chrom,beg) in clust_norm_p_end) and ((chrom,end) in clust_norm_m_beg):
                    insertionType = "TD0"                    
                    chromPlus, begPlus, endPlus, nbReadsPlus, classPlus, readListPlus = clust_norm_p_end[chrom,beg]
                    chromMinus, begMinus, endMinus, nbReadsMinus, classMinus, readListMinus = clust_norm_m_beg[chrom,end]
                    
                    cytobandId, sourceType, chromSource, begSource, endSource, strandSource, tdBeg, tdEnd, tdRnaLen, tdLen, psdGene, chromExonA, begExonA, endExonA, chromExonB, begExonB, endExonB = ["NA"] *  17

                ## B) Orphan transduction. Transduction-like positive and negative clusters
                elif ((chrom,beg) in clust_del_p_end) and ((chrom,end) in clust_del_m_beg):
                    insertionType = "TD2"
                    chromPlus, begPlus, endPlus, nbReadsPlus, classPlus, readListPlus, chromSource, begSource, endSource, strandSource, tdChrom, tdBegPlus, tdEndPlus = clust_del_p_end[chrom,beg]
                    chromMinus, begMinus, endMinus, nbReadsMinus, classMinus, readListMinus, chromSource, begSource, endSource, strandSource, tdChrom, tdBegMinus, tdEndMinus = clust_del_m_beg[chrom,end]

                    tdBeg = min([tdBegPlus, tdBegMinus])
                    tdEnd = max([tdEndPlus, tdEndMinus])

                    tdCoords = tdChrom + ":" + tdBeg + "-" + tdEnd
                    cytobandId, sourceType, begSource, endSource = isGermlineSource(tdCoords, sourceMetadataDict)
                    strandSource, tdBeg, tdEnd, tdRnaLen = computeTdRnaLen(begSource, endSource, strandSource, tdBeg, tdEnd)
                    tdLen = str(int(tdEnd) - int(tdBeg)) 


                ## C) Partnered transduction. Normal positive cluster and transduction-like negative cluster
                elif ((chrom,beg) in clust_norm_p_end) and ((chrom,end) in clust_del_m_beg):
                    insertionType = "TD1"
                    chromPlus, begPlus, endPlus, nbReadsPlus, classPlus, readListPlus = clust_norm_p_end[chrom,beg]
                    chromMinus, begMinus, endMinus, nbReadsMinus, classMinus, readListMinus, chromSource, begSource, endSource, strandSource, tdChrom, tdBeg, tdEnd = clust_del_m_beg[chrom,end]

                    tdCoords = tdChrom + ":" + tdBeg + "-" + tdEnd
                    cytobandId, sourceType, begSource, endSource = isGermlineSource(tdCoords, sourceMetadataDict)
                    strandSource, tdBeg, tdEnd, tdRnaLen = computeTdRnaLen(begSource, endSource, strandSource, tdBeg, tdEnd)
                    tdLen = tdRnaLen

                ## D) Partnered transduction. Transduction-like positive cluster and normal negative cluster
                elif ((chrom,beg) in clust_del_p_end) and ((chrom,end) in clust_norm_m_beg):
                    insertionType = "TD1"
                    chromPlus, begPlus, endPlus, nbReadsPlus, classPlus, readListPlus, chromSource, begSource, endSource, strandSource, tdChrom, tdBeg, tdEnd = clust_del_p_end[chrom,beg]
                    chromMinus, begMinus, endMinus, nbReadsMinus, classMinus, readListMinus = clust_norm_m_beg[chrom,end]

                    tdCoords = tdChrom + ":" + tdBeg + "-" + tdEnd
                    cytobandId, sourceType, begSource, endSource = isGermlineSource(tdCoords, sourceMetadataDict)
                    strandSource, tdBeg, tdEnd, tdRnaLen = computeTdRnaLen(begSource, endSource, strandSource, tdBeg, tdEnd)
                    tdLen = tdRnaLen

                #print("# %s, %s, %s, %s, %s, %s #" % (study, donor, sample, chrom, beg, end))
                if insertionType == "NA":
                    log("Missing cluster info for deletion", "WARN")
                    num_events_missed += 1
                    continue

                outline = chromPlus + "\t" + begPlus + "\t" + endPlus + "\t" + nbReadsPlus + "\t" + classPlus + "\t" + readListPlus + "\t" + chromMinus + "\t" + begMinus + "\t" + endMinus + "\t" + nbReadsMinus + "\t" + classMinus + "\t" + readListMinus + "\t" + insertionType + "\t" + cytobandId + "\t" + sourceType + "\t" + chromSource + "\t" + begSource + "\t" + endSource + "\t" + strandSource + "\t" + tdBeg + "\t" + tdEnd + "\t" + tdRnaLen + "\t" + tdLen + "\t" + psdGene + "\t" + chromExonA + "\t" + begExonA + "\t" + endExonA + "\t" + chromExonB + "\t" + begExonB + "\t" + endExonB + "\t" + "DEL"         
                print(outline, file=outfile)

    log("Skipped events: %d (of %d)" % (num_events_missed, num_events))



