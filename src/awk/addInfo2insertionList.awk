#!/usr/bin/env awk

# ******************************************************************************

#         addInfo2insertionList.awk
# 
#         Copyright (c) 2016 Bernardo Rodríguez-Martín

#         Mobile Genomes & Disease Lab.
#         Universidad de Vigo (Spain)

#         Licenced under the GNU General Public License 3.0 license.
# ******************************************************************************

# Description
################

# Add insertion information to a list of insertion ids.

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
# 13. insertion_type
## transduction-specific fields:
# 14. chrom_source_element, "NA" for TD0 
# 15. beg_source_element, "NA" for TD0
# 16. end_source_element, "NA" for TD0
# 17. orientation_source_element, "NA" for TD0
# 18. transduction_beg, "NA" for TD0 and TD1
# 19. transduction_end, "NA" for TD0
# 20. transduction_rna_length, "NA" for TD0
# 21. transduction_length, "NA" for TD0
## pseudogene-specific fields ("NA" for TD* events)
# 22. psd_gene
# 23. chr_exon_a
# 24. beg_exon_a
# 25. end_exon_a
# 26. chr_exon_b
# 27. beg_exon_b
# 28. end_exon_b

## insertionsList.txt:
# L1:TD0:1_29590838_29590787
# L1:TD1:4_847239_847339
# L1:TD2:8_543209_543509
# NA:PSD:2_41370097_41370069
# ...

#### Output
# L1:TD0:1_29590838_29590787   ${readList_plus}    ${readList_minus}    NA  NA
# L1:TD1:4_847239_847339   ${readList_plus}    ${readList_minus}    ${chrom_source}_${beg_source}-${end_source}_${orientation_source}   ${transduct_beg}:${transduct_end}:${transduct_rna_len}:${transduct_len}   NA:NA_NA-NA:NA_NA-NA
# L1:TD2:8_543209_543509   ${readList_plus}    ${readList_minus}    ${chrom_source}_${beg_source}-${end_source}_${orientation_source}   NA:${transduct_end}:${transduct_rna_len}:${transduct_len}   NA:NA_NA-NA:NA_NA-NA
# NA:PSD:2_41370097_41370069   ${readList_plus}    ${readList_minus}    NA:NA-NA_NA   NA:NA:NA:NA   ${psd_gene}:${chr_exon_a}_${beg_exon_a}-${end_exon_a}:${chr_exon_b}_${beg_exon_b}-${end_exon_b}
# ...

# Notes:
# - TD0: solo-insertion, TD1: partnered-transduction and TD2: orphan-transduction
# - For partnered-transductions the ${transduct_rna_len} is equal to ${transduct_len}.
# - NA: not applicable


BEGIN{
    while (getline < fileRef >0)
    {
        ##### Get input
        ### Insertion coordinates
        chrPlus = $1;
        endPlus = $3 + 100; # Done because this coordinate is the beginning of the read. So, we need to sum the readlength
                    # (this is provisional, input file should already have this corrected)
        begMinus = $8;

        ### Family:
        if ($5 == "Other")
        {
            familyPlus = "SVA";
        }
        else
        {
            familyPlus = $5;
        }

        ### Supporting read pair identifiers list:
        readPairsPlus = $6;
        readPairsMinus = $12;

        ### Insertion type
        tdType = $13

        ### Transduction info
        # A) Solo insertion
        if  (tdType == "TD0")
        {
            sourceElementInfo = "NA"
            transductionInfo = "NA"
            pseudogeneInfo = "NA"
        }
        # B) Partnered or orphan transduction
        else if (tdType == "TD1" || tdType == "TD2")
        {
            chromSource = $14
            begSource = $15
            endSource = $16
            orientationSource = $17
            transductBeg = $18
            transductEnd = $19
            transductRnaLen = $20
            transductLen = $21
            sourceElementInfo = chromSource"_"begSource"_"endSource"_"orientationSource
            transductionInfo = chromSource"_"transductBeg"_"transductEnd"_"transductRnaLen"_"transductLen
            pseudogeneInfo = "NA"
        }
        # C) pseudogene insertion
        else if (tdType == "PSD")
        {
            psdGene = $22
            chrExonA = $23
            begExonA = $24
            endExonA = $25
            chrExonB = $26
            begExonB = $27
            endExonB = $28

            sourceElementInfo = "NA"
            transductionInfo = "NA"
            pseudogeneInfo = psdGene":"chrExonA"_"begExonA"-"endExonA":"chrExonB"_"begExonB"-"endExonB
        }
        # *) unknown tdType
        else
        {
            sourceElementInfo = "NA"
            transductionInfo = "NA"
            pseudogeneInfo = "NA"
        }    

        ## Keep info in dictionaries
        insertionId = familyPlus":"tdType":"chrPlus"_"endPlus"_"begMinus;

        readPairsPlusDict[insertionId] = readPairsPlus;
        readPairsMinusDict[insertionId] = readPairsMinus;
        sourceElementInfoDict[insertionId] = sourceElementInfo
        transductionInfoDict[insertionId] = transductionInfo
        pseudogeneInfoDict[insertionId] = pseudogeneInfo
    }
}
{
    insertionId = $1;

    # Gather info from dictionaries and print output:
    row = insertionId"\t"readPairsPlusDict[insertionId]"\t"readPairsMinusDict[insertionId]"\t"sourceElementInfoDict[insertionId]"\t"transductionInfoDict[insertionId]"\t"pseudogeneInfoDict[insertionId];
    print row;
}
