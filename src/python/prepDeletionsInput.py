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
out_line_fmt = "%s" + ("\t%s" * 12) + ("\tNA" * 15)


def log(msg, event_type="INFO"):
    outfile = sys.stdout
    if event_type == "ERROR":
        outfile = sys.stderr
    print("[%s] (prepDeletionsInput.py) %s" % (event_type, msg), file=outfile)

def parse_args():
    parser = argparse.ArgumentParser(description="""Generates TEIBA input files for deletions from TraFiC output.""")
    parser.add_argument('eventsFile', help='TraFiC output file containing the coordinates of deletion events (one file for all samples).')
    parser.add_argument('samplesPath', help='Directory that contains sample folders.')
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
            sys.exit(1)

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

        # parse sample input files
        #   columns: chrom, beg, end, num, class, reads, ...
        clust_norm_p_end = {}
        #print("## %s, %s, %s ##" % (study, donor, sample))
        for line in open(fn_norm_p, 'rt'):
            #print("+#%s#" % line)
            row = line.strip().split()
            chrom, beg, end, num, klass, reads = row[:6]
            clust_norm_p_end[(chrom, end)] = (chrom, beg, end, num, klass, reads)
        clust_norm_m_beg = {}
        for line in open(fn_norm_m, 'rt'):
            #print("-#%s#" % line)
            row = line.strip().split()
            chrom, beg, end, num, klass, reads = row[:6]
            clust_norm_m_beg[(chrom, beg)] = (chrom, beg, end, num, klass, reads)
        clust_del_p_end = {}
        for line in open(fn_del_p, 'rt'):
            #print("+#%s#" % line)
            row = line.strip().split()
            chrom, beg, end, num, klass, reads = row[:6]
            clust_del_p_end[(chrom, end)] = (chrom, beg, end, num, klass, reads)
        clust_del_m_beg = {}
        for line in open(fn_del_m, 'rt'):
            #print("-#%s#" % line)
            row = line.strip().split()
            chrom, beg, end, num, klass, reads = row[:6]
            clust_del_m_beg[(chrom, beg)] = (chrom, beg, end, num, klass, reads)


        # combine cluster info for registered deletions
        id_sample = (study, donor, sample)
        with open(os.path.join(outdir, "L1-del.txt"), 'wt') as outfile:
            for chrom, beg, end in samples_events[id_sample]:
                td_type = "NA"
                # determine transduction type (TD0: solo, TD1: partnered, TD2: orphan)
                if ((chrom,beg) in clust_norm_p_end) and ((chrom,end) in clust_norm_m_beg):
                    td_type = "TD0"
                    outline = out_line_fmt % sum((clust_norm_p_end[chrom,beg], clust_norm_m_beg[chrom,end], (td_type,)), ())
                elif ((chrom,beg) in clust_del_p_end) and ((chrom,end) in clust_del_m_beg):
                    td_type = "TD2"
                    outline = out_line_fmt % sum((clust_del_p_end[chrom,beg], clust_del_m_beg[chrom,end], (td_type,)), ())
                elif ((chrom,beg) in clust_norm_p_end) and ((chrom,end) in clust_del_m_beg):
                    td_type = "TD1"
                    outline = out_line_fmt % sum((clust_norm_p_end[chrom,beg], clust_del_m_beg[chrom,end], (td_type,)), ())
                elif ((chrom,beg) in clust_del_p_end) and ((chrom,end) in clust_norm_m_beg):
                    td_type = "TD1"
                    outline = out_line_fmt % sum((clust_del_p_end[chrom,beg], clust_norm_m_beg[chrom,end], (td_type,)), ())
                #print("# %s, %s, %s, %s, %s, %s #" % (study, donor, sample, chrom, beg, end))
                if td_type == "NA":
                    log("Missing cluster info for deletion", "WARN")
                    num_events_missed += 1
                    continue


                print(outline, file=outfile)

    log("Skipped events: %d (of %d)" % (num_events_missed, num_events))
