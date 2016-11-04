#!/usr/bin/env awk

# ******************************************************************************
        
#         vcf2bed.MEI.awk
#         
#         Copyright (c) 2016 Bernardo Rodríguez-Martín
        
#         Mobile Genomes & Disease Lab.
#         Universidad de Vigo (Spain)

#         Licenced under the GNU General Public License 3.0 license.
# ******************************************************************************

# Description
################ 
        
# Takes as input a vcf containing TraFiC MEI and convert insertion calls into bed format

# Usage
#######

# awk -f vcf2bed.MEI.awk input.vcf

#### Input (VCF with MEI calls):
##fileformat=VCFv4.2
##fileDate=20160620 
##source=TraFiCvX 
##reference=hs37d5
##contig=<ID=1,assembly=GRCh37,length=249250621,species=human>
##contig=<ID=2,assembly=GRCh37,length=243199373,species=human>
##contig=<ID=3,assembly=GRCh37,length=198022430,species=human>
# ........
# 10 7440548 SVTYPE=<MEI>;CLASS=L1;TYPE=TD0;SCORE=2A;CIPOS=0;STRAND=-;STRUCT=DEL;LEN=1641;CONTIGA=GAATCAATATCGTGTAAATGGCCATACTGCCCAAGGTAATTTACACATTTAATGCCATCCCCATCAAGTTACAAATGACTTTCTTCACAGAATTGGAAAAAACTACTTTAAAGTTCATATGGAACCAAAAAAGAGCCCGCATTGCCAAGTCAATCCTAAGCCAAAAGAACAAAGCTGGAGGCATCACACTACCTGACTTCAAACTATACTACAAGGCTATAGTAACCAAAACAGCATGGTACTGGTACCAAAACAGAGATATAGATCAATGGAACAGAACAGAGCCCTCAGAAATAATGCCGCATATCTACAACTATCTGATCTTTGACAAACCTGAGAAAAACAAGCAATGGGGAAAGGATTCCCTATTTAATAAATGGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTTTTTTTTATTTTTGTGGGTACATGGTGTAAATGGGGCCTTATTTAAAATGCAGCAGTTTTAGAATTAGCAGAGATGTAGTCACACTGGGTTCAAACCTATTGGACAAACCCTACATACAGGATTTCTCTGTCCCTTATAAATATTAAAAGGAATAAATTAGAGCATCCCCCCAGGTTTGTAAAAAGAGCTAAGTAAACCTTATATTCCCCAAAATATTTGTGTAATGCATTACATACAAATGCAGATAAGTCATTTTAAGTAGATCATGTGTTGTCCATACCGCCATGAGACTCAAAAAGTAAAGAAACTTACAAGGAAGGCTTCTGTAAGCCAGACATGCATGGATAATGCCTAATTAGGGCATAATTTTCTGGGTGGATTTTAAGAATTGATCACTTCAAGCCAAACAAAC
# 1 246817033 SVTYPE=<MEI>;CLASS=L1;TYPE=TD0;SCORE=3;CIPOS=30

#### Output (Text file containing MEI calls in proper format for variant annotation with Annovar):
# 10       7440548 7440548
#  1        246817003       246817063

function abs(v) {
    return v < 0 ? -v : v
}

! /^#/ {
    # Get input 
    chrom = $1;
    pos = $2;
    info = $8;
    
    # Parse info field and select class and CIPOS fields
    split(info, infoList, ";");

    for (field in infoList)
    {

        split(infoList[field], fieldList, "=");
        tag = fieldList[1];
        value = fieldList[2];

        # A) Class
        if (tag == "CLASS")
        {
            class = value
        }

        # B) Score
        else if (tag == "SCORE")
        {
            score = value
        }        
    
        # C) CIPOS (confindence interval around position)
        else if (tag == "CIPOS")
        {
            CIPOS = value;    
        }

        # C) TSLEN (target site length)
        else if (tag == "TSLEN")
        {
            targetSiteLen = abs(value); # Absolute value (deletion length -> positive integer)    
        }
    }    

    # A) Target site identified
    if ( score == "5")
    {
        # Compute beg and end
        beg = pos - 1;    # Substract 1 to convert from 1-based (VCF) to 0-based (bed) coordinate system
        end = pos + targetSiteLen;                
    }
    # B) Target site not identified
    else
    {
        # Compute beg and end
        beg = pos - CIPOS - 1;    # Substract 1 to convert from 1-based (VCF) to 0-based (bed) coordinate system
        end = pos + CIPOS;
    }

    # Print bed row
    row=chrom"\t"beg"\t"end"\t"class;        
    print row;
}
