#!/usr/bin/env awk

# ******************************************************************************

#         addOffset2template.awk
# 
#         Copyright (c) 2016 Bernardo Rodríguez-Martín

#         Mobile Genomes & Disease Lab.
#         Universidad de Vigo (Spain)

#         Licensed under the GNU General Public License 3.0 license.
# ******************************************************************************

# Description
################


# Usage
#######

# awk -v OFS='\t' -f  addOffset2template.awk L1:8_37122847_37122890:-.tmp.psl

## Input:
# 286    0    0    0    0    0    0    0    +    NODE_2_length_614_cov_6.804560    634    0    286    L1    6021    5735    6021    1    286,    0,    5735,
# 330    0    0    0    0    0    0    0    +    NODE_2_length_614_cov_6.804560    634    304    634    8:37121847-37123890    2044    1043    1373    1    330,    304,    1043,

{
    # Get input
    matches=$1;
    misMatches=$2;
    repMatches=$3;
    nCount=$4;
    qNumInsert=$5;
    qBaseInsert=$6;
    tNumInsert=$7;
    tBaseInsert=$8;
    strand=$9;
    qName=$10;
    qSize=$11;
    qStart=$12;
    qEnd=$13;
    tName=$14;
    tSize=$15;
    tStart=$16;
    tEnd=$17;
    blockCount=$18;
    blockSizes=$19;
    qStarts=$20;
    tStarts=$21;

    # get chromosome id and coordinates (if present) from template id
    split(tName, name, ":");
    tName=name[1];

    if (length(name)==2)
    {
        split(name[2], coord, "-");

        if (coord[1]~/[0-9]+/ && coord[2]~/[0-9]+/) {

            offset = coord[1] - 1;

            # Add offset to template start and end (chromosomic coordinates)
            tStart=tStart + offset;
            tEnd=tEnd + offset;

            # Add offset to the list of template starts
            nb=split(tStarts, startsList, ",");

            tStarts="";

            for (counter=1; counter<=(nb-1); counter++)
            {
                startPos = startsList[counter];

                if (tStarts=="")
                {
                    startPos=startPos + offset;
                    tStarts=startPos",";
                }
                else
                {
                    startPos=startPos + offset;
                    tStarts=tStarts""startPos",";
                }
            }
        }
    }

    # Print a psl row with updated template positions
row=matches"\t"misMatches"\t"repMatches"\t"nCount"\t"qNumInsert"\t"qBaseInsert"\t"tNumInsert"\t"tBaseInsert"\t"strand"\t"qName"\t"qSize"\t"qStart"\t"qEnd"\t"tName"\t"tSize"\t"tStart"\t"tEnd"\t"blockCount"\t"blockSizes"\t"qStarts"\t"tStarts;
print row;
}
