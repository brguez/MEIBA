#!/usr/bin/env python
#coding: utf-8 

import re


class VCF():
    """ 
	.....................
    
	Methods:
	- 

    """

    def __init__(self):
	""" 
	    ...
            
	    Output:
            - 
	"""
	self.header = "" 
	self.lineList = []  # List of VCFline objects
        	
    #### METHODS ####
    def read_VCF(self, VCFpath):

	VCFfile = open(VCFpath, 'r')
	headerList = []

	for line in VCFfile:
	    
	    # A) Header
	    if line.startswith("#"):
		headerList.append(line)
	    
	    # B) Variants
	    else:
   	        line = line.rstrip('\n')
	        line = line.split('\t')
		VCFlineObj = VCFline(line)
		self.addLine(VCFlineObj)

	self.header = "".join(headerList)

    def addLine(self, VCFlineObj):
	""" 

	"""
	self.lineList.append(VCFlineObj)

    def write_header(self, outFilePath):
	""" 
	    ...
            
	    Input:
	    1) ...
            
	    Output:
	    1) ...
	"""

	outFile = open(outFilePath, 'w')	
	outFile.write(self.header)
	outFile.close()

    def write_variants(self, outFilePath):
	""" 
	    ...
            
	    Input:
	    1) ...
            
	    Output:
	    1) ...
	"""

	outFile = open(outFilePath, 'a')	

	# Iterate and print each VCF line into the output VCF file
	for VCFline in self.lineList:

	   row = VCFline.chrom + "\t" + str(VCFline.pos) + "\t" + VCFline.id + "\t" + VCFline.ref + "\t" + VCFline.alt + "\t" + VCFline.qual + "\t" + VCFline.filter + "\t" + VCFline.info + "\n"
           outFile.write(row)

	## Close output file
	outFile.close()


class VCFline():
    """ 
	.....................
    
	Methods:
	- 

    """

    def __init__(self, VCFlineList):
	""" 
	    ...
            
	    Output:
            - 
	"""

	self.chrom = VCFlineList[0]
	self.pos = int(VCFlineList[1])
	self.id = VCFlineList[2]
	self.ref = VCFlineList[3]  
	self.alt = VCFlineList[4]
	self.qual = VCFlineList[5] 
	self.filter = VCFlineList[6] 
	self.info = VCFlineList[7]
	self.infoDict = self.read_info()

    def read_info(self):
	""" 

	"""	
	infoDict = {}
	infoList = self.info.split(';')
	
	for field in infoList:
	    fieldList = field.split('=')
	    key = fieldList[0]
	    value = fieldList[1]
 	    infoDict[key] = value

        return infoDict

    def make_info(self):
	""" 
	"""
	
	## Create list containing the order of info fields (provisional)
	infoOrder = [ "SVTYPE", "CLASS", "TYPE", "SCORE", "CIPOS", "STRAND", "STRUCT", "LEN", "TSLEN", "TSSEQ", "POLYA", "REGION", "GENE", "REP", "DIV", "CONTIGA", "CONTIGB", "TRDS" ] 

	## Create info string in the correct order from dictionary 
	infoList = []

	for info in infoOrder:
		
	    if (info in self.infoDict.keys()) and (self.infoDict[info] != "unkn"):
  
		infoField = info + "=" +str(self.infoDict[info])  
		infoList.append(infoField)
	    
	info = ';'.join(infoList)
	  
	return(info)


class annovar():
    """ 
	.....................
    
	Methods:
	- 

    """

    def __init__(self):
	""" 
	    ...
            
	    Output:
            - 
	"""
	self.lineList = []  # List of annovarLine objects
        	
    #### METHODS ####
    def read_annovar(self, annovarPath):
	""" 

	"""
	annovarFile = open(annovarPath, 'r')

	for line in annovarFile:
	    
	    line = line.rstrip('\n')
	    annovarLineList = line.split('\t')
            annovarLineObj = annovarLine(annovarLineList)
   	    self.addLine(annovarLineObj)

    def addLine(self, annovarLineObj):
	""" 

	"""
	self.lineList.append(annovarLineObj)

    def info2VCF(self, VCFObj):
	""" 

	"""

	annotDictList = []
	
	for i in range(0, len(self.lineList), 1):
      
	    annovarLineObj = self.lineList[i]
	    VCFlineObj =  VCFObj.lineList[i]

	    if (annovarLineObj.chrom == VCFlineObj.chrom) and (annovarLineObj.beg == VCFlineObj.pos):
		
		annotDict = {}
	    	annotDict["REGION"] = annovarLineObj.region

		if (annotDict["REGION"] != "intergenic"):		
		    annotDict["GENE"] = annovarLineObj.genes 		

		VCFlineObj.infoDict.update(annotDict)
		VCFlineObj.info = VCFlineObj.make_info()

	    else:
		print "ERROR. Unpaired VCF and annovar variants"
		

class annovarLine():
    """ 
	.....................
    
	Methods:
	- 

    """

    def __init__(self, annovarLineList):
	""" 
	    ...
            
	    Output:
            - 
	"""
	self.region = annovarLineList[0]
	self.genes = re.sub('\(.*?\)', '', annovarLineList[1]) # remove not necesary information enclosed by slashes (.*)
	self.chrom = annovarLineList[2]
	self.beg = int(annovarLineList[3])
	self.end = int(annovarLineList[4])
	self.ref = annovarLineList[5]
	self.alt = annovarLineList[6]

