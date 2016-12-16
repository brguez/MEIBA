#!/usr/bin/env awk

# ******************************************************************************
        
#         fixFormatGermlineMEI.awk
#         
#         Copyright (c) 2016 Bernardo Rodríguez-Martín
        
#         Mobile Genomes & Disease Lab.
#         Universidad de Vigo (Spain)

#         Licenced under the GNU General Public License 3.0 license.
# ******************************************************************************

# Description
################ 
        
# Fix format of the germlin insertions file for making it compatible with current TEIBA version input format

# Usage
#######

#### Input:

### TraFic MEI, tsv with the following format:
# 1. chrom_plus
# 2. beg_plus
# 3. end_plus
# 4. nbReads_plus
# 5. class_plus
# 6. readList_plus
# 7. chrom_minus
# 8. beg_minus
# 9. end_minus
# 10. nbReads_minus
# 11. class_minus
# 12. readList_minus

### Output
## Appends to each input row the following 9 rows
# 13. insertion_type, "TD0"
# 14. chrom_source_element, "NA"  
# 15. beg_source_element, "NA" 
# 16. end_source_element, "NA" 
# 17. orientation_source_element, "NA" 
# 18. transduction_beg, "NA"  
# 19. transduction_end, "NA" 
# 20. transduction_rna_length, "NA" 
# 21. transduction_length, "NA" 

# Notes:
# - TD0: solo-insertion, TD1: partnered-transduction and TD2: orphan-transduction
# - For partnered-transductions the ${transduct_rna_len} is equal to ${transduct_len}. 
# - NA: not applicable
 

{
    # Get input
    chromPlus=$1;
    begPlus=$2;
    endPlus=$3;
    nbReadsPlus=$4;
    classPlus=$5;
    readListPlus=$6;
    chromMinus=$7;
    begMinus=$8;
    endMinus=$9;
    nbReadsMinus=$10;
    classMinus=$11;
    readListMinus=$12;

    # New fields:
    type="TD0";
    chromSrc="NA";
    begSrc="NA"; 
    endSrc="NA"; 
    orientationSrc="NA"; 
    tdBeg="NA";  
    tdEnd="NA";
    tdRNALen="NA";
    tdLen="NA";

    # Produce output
        row=chromPlus"\t"begPlus"\t"endPlus"\t"nbReadsPlus"\t"classPlus"\t"readListPlus"\t"chromMinus"\t"begMinus"\t"endMinus"\t"nbReadsMinus"\t"classMinus"\t"readListMinus"\t"type"\t"chromSrc"\t"begSrc"\t"endSrc"\t"orientationSrc"\t"tdBeg"\t"tdEnd"\t"tdRNALen"\t"tdLen
    print row;
}

