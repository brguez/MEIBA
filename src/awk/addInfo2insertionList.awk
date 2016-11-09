#!/usr/bin/env awk

# ******************************************************************************
        
#         addSupReads2insertionList.awk
#         
#         Copyright (c) 2016 Bernardo Rodríguez-Martín
        
#         Mobile Genomes & Disease Lab.
#         Universidad de Vigo (Spain)

#         Licenced under the GNU General Public License 3.0 license.
# ******************************************************************************

# Description
################ 
        
# 

# Usage
#######

# awk 

#### Input:

## TraFic MEI: 
# 1    29590660    29590738    9    Alu    HWI-ST1324:244:C2YAEACXX:1:2103:18894:41914,HWI-ST896:446:C2Y9RACXX:5:2215:8387:90585,HWI-ST896:443:C2P4CACXX:5:1110:1052:31345,HWI-ST898:490:C2RHCACXX:8:2111:13526:98452,HWI-ST898:490:C2RHCACXX:8:2214:3957:43220,HWI-ST896:446:C2Y9RACXX:5:1313:17751:47531,HWI-ST898:490:C2RHCACXX:7:1106:13592:35855,HWI-ST896:446:C2Y9RACXX:5:1102:13259:64695,HWI-ST896:443:C2P4CACXX:5:2102:3584:81108    1    29590787    29591048    13    Alu    HWI-ST898:490:C2RHCACXX:8:2205:7992:57333,HWI-ST896:443:C2P4CACXX:5:2102:19692:51065,HWI-ST896:443:C2P4CACXX:5:2310:11874:44940,HWI-ST898:490:C2RHCACXX:7:2305:4132:90087,HWI-ST896:446:C2Y9RACXX:5:2205:2438:3595,HWI-ST896:446:C2Y9RACXX:5:1112:17167:72893,HWI-ST1324:244:C2YAEACXX:1:1205:5736:57799,HWI-ST896:446:C2Y9RACXX:5:2306:5122:26461,HWI-ST1324:244:C2YAEACXX:1:2301:2327:91679,HWI-ST898:490:C2RHCACXX:8:1308:18526:27394,HWI-ST898:490:C2RHCACXX:7:2304:7219:50597,HWI-ST1324:244:C2YAEACXX:1:2213:8365:31759,HWI-ST898:490:C2RHCACXX:8:1312:5658:33944
# ...

## insertionsList.txt: 
# Alu:1_29590838_29590787
# ...

#### Output 
# Alu:1_29590838_29590787    HWI-ST1324:244:C2YAEACXX:1:2103:18894:41914,HWI-ST896:446:C2Y9RACXX:5:2215:8387:90585,HWI-ST896:443:C2P4CACXX:5:1110:1052:31345,HWI-ST898:490:C2RHCACXX:8:2111:13526:98452,HWI-ST898:490:C2RHCACXX:8:2214:3957:43220,HWI-ST896:446:C2Y9RACXX:5:1313:17751:47531,HWI-ST898:490:C2RHCACXX:7:1106:13592:35855,HWI-ST896:446:C2Y9RACXX:5:1102:13259:64695,HWI-ST896:443:C2P4CACXX:5:2102:3584:81108    HWI-ST898:490:C2RHCACXX:8:2205:7992:57333,HWI-ST896:443:C2P4CACXX:5:2102:19692:51065,HWI-ST896:443:C2P4CACXX:5:2310:11874:44940,HWI-ST898:490:C2RHCACXX:7:2305:4132:90087,HWI-ST896:446:C2Y9RACXX:5:2205:2438:3595,HWI-ST896:446:C2Y9RACXX:5:1112:17167:72893,HWI-ST1324:244:C2YAEACXX:1:1205:5736:57799,HWI-ST896:446:C2Y9RACXX:5:2306:5122:26461,HWI-ST1324:244:C2YAEACXX:1:2301:2327:91679,HWI-ST898:490:C2RHCACXX:8:1308:18526:27394,HWI-ST898:490:C2RHCACXX:7:2304:7219:50597,HWI-ST1324:244:C2YAEACXX:1:2213:8365:31759,HWI-ST898:490:C2RHCACXX:8:1312:5658:33944
# ...

BEGIN{
    while (getline < fileRef >0)
    {
        ## Get input 
        # Insertion coordinates        
        chrPlus = $1;
        endPlus = $3 + 100; # Done because this coordinate is the beginning of the read. So, we need to sum the readlength 
                    # (this is provisional, input file should already have this corrected)
        begMinus = $8;
 
        # Family:
        if ($5 == "Other")
        {
            familyPlus = "SVA";
        }        
        else
        {
            familyPlus = $5;
        }
        
        # Supporting read pair identifiers list:
        readPairsPlus = $6;
        readPairsMinus = $12;

        ## Keep info in dictionaries
        insertionId = familyPlus":"chrPlus"_"endPlus"_"begMinus;        

        readPairsPlusDict[insertionId] = readPairsPlus; 
        readPairsMinusDict[insertionId] = readPairsMinus; 
    }
} 
{
    insertionId = $1;

    # Gather info from dictionaries and print output:
    row = insertionId"\t"readPairsPlusDict[insertionId]"\t"readPairsMinusDict[insertionId];   
    print row;
}

