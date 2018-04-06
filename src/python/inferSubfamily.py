#!/usr/bin/env python
#coding: utf-8


#### CLASSES ####
class fasta():
    """
    """

    def __init__(self):
        """
        """
        self.fastaDict = {}

    #### FUNCTIONS ####
    def fasta_reader(self, fastaFile):
        """
        """
        fastaDict = {}

        subHeader("Fasta reader")

        fh = open(fastaFile)
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()

            # drop the info
            header = header.split(" ")[0]

            info("Reading " + header + "...")
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            fastaDict[header] = seq

        self.fastaDict = fastaDict

    def write_fasta(self, outFilePath):
        """
        """
        outFile = open(outFilePath, "w" )

        for header, seq in self.fastaDict.iteritems():
            header = ">" + header

            outFile.write("%s\n" % header)
            outFile.write("%s\n" % seq)

        # Close output fasta file
        outFile.close()

#### FUNCTIONS ####
def header(string):
    """
        Display  header
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print '\n', timeInfo, "****", string, "****"

def subHeader(string):
    """
        Display  subheader
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, "**", string, "**"

def log(label, string):
    """
        Display labelled information
    """
    print "[" + label + "]", string

def info(string):
    """
        Display basic information
    """
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print timeInfo, string


def retrieveAnchorMates(chrom, beg, end, discordantIdList, bamFile):
    '''
    '''
    anchorMateList = []

    ## Open bam file for reading:
    bamFile = pysam.AlignmentFile(bamFile, "rb")

    ## Extract alignments in the interval
    iterator = bamFile.fetch(chrom, beg, end)

    ## Iterate over the alignments
    for alignmentObj in iterator:

        ## Supporting discordant paired-end read
        if alignmentObj.query_name in discordantIdList:

            anchorMateId = alignmentObj.query_name + '/2' if alignmentObj.is_read1 else alignmentObj.query_name + '/1'
            anchorMateList.append(anchorMateId)

    bamFile.close()

    return anchorMateList

def subFamilyL1(vcf)  :
    """
    """

    VCFObj = formats.VCF()
    VCFObj.read_VCF(vcf)

    diagnosticNtDict = {}
    diagnosticNtDict[5535] = "NA"
    diagnosticNtDict[5538] = "NA"
    diagnosticNtDict[5929] = "NA"
    diagnosticNtDict[5930] = "NA"
    diagnosticNtDict[5931] = "NA"


    ### 1. Interrogate diagnostic nucleotide positions
    ## A) At least one variant in the VCF
    if (len(VCFObj.lineList) > 0):
        
        # Iterate over each position in the VCF
        for VCFlineObj in VCFObj.lineList:

            ## If diagnostic nucleotide position:
            if (int(VCFlineObj.pos) in diagnosticNtDict):

                # a) Alternative nucleotide not available
                if (VCFlineObj.alt != "."):

                    diagnosticNtDict[int(VCFlineObj.pos)] = VCFlineObj.alt

                # b) Alternative nucleotide available
                else:
                    diagnosticNtDict[int(VCFlineObj.pos)] = VCFlineObj.ref

    ## B) No variant in the VCF     
    else:   
        subFamily = "NA"

    ### 2. Determine element subfamily based on diagnostic nucleotides status:
    ## A) pre-Ta (diagnostic trinucleotide ACG at 5,929-5,931)
    if (diagnosticNtDict[5929] == "A") and (diagnosticNtDict[5930] == "C") and (diagnosticNtDict[5931] == "G"):

        subFamily = "pre-Ta"

    ## B) Ta (diagnostic trinucleotide ACA at 5,929-5,931)   
    elif (diagnosticNtDict[5929] == "A") and (diagnosticNtDict[5930] == "C") and (diagnosticNtDict[5931] == "A"):

        ## A.a) Ta-0 (G and C at 5,535 and 5,538)
        if (diagnosticNtDict[5535] == "G") and (diagnosticNtDict[5538] == "C"):
            subFamily = "Ta-0"

        ## A.b) Ta-1 (T and G at 5,535 and 5,538)
        elif (diagnosticNtDict[5535] == "T") and (diagnosticNtDict[5538] == "G"):
            subFamily = "Ta-1"
    
        else:
            subFamily = "Ta"

    ## C) Unknown
    else:
        subFamily = "NA"
  

    return(subFamily)

        

#### MAIN ####

## Import modules ##
import argparse
import time
import sys
import os.path
import formats
import os.path
import itertools 
import pysam

## Get user's input ##
parser = argparse.ArgumentParser(description= "")
parser.add_argument('VCF', help='VCF file')
parser.add_argument('bamFile', help='Tumour bam file')
parser.add_argument('fastaFile', help='Fasta file containing insertion supporting reads')
parser.add_argument('refDir', help='Directory containing transposon reference sequences')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

args = parser.parse_args()
inputVCF = args.VCF
bamFile = args.bamFile
fastaFile = args.fastaFile
refDir = args.refDir
fileName = args.fileName
outDir = args.outDir

scriptName = os.path.basename(sys.argv[0])

## Display configuration to standard output ##
print
print "***** ", scriptName, " configuration *****"
print "vcf: ", inputVCF
print "bamFile: ", bamFile
print "fastaFile: ", fastaFile
print "refDir: ", refDir
print "fileName: ", fileName
print "outDir: ", outDir
print
print "***** Executing ", scriptName, ".... *****"
print

## Start ## 

#### 1. Read fasta file containing supporting reads
fastaObj = fasta()
fastaObj.fasta_reader(fastaFile)

#### 2. Read VCF file
VCFObj = formats.VCF()
VCFObj.read_VCF(inputVCF)

## For each MEI
for VCFlineObj in VCFObj.lineList:

    ### Define insertion Id    
    chrom = VCFlineObj.chrom
    pos = str(VCFlineObj.pos)
    insertionType = VCFlineObj.infoDict["TYPE"]
    family = VCFlineObj.infoDict["CLASS"] if "CLASS" in VCFlineObj.infoDict else 'NA'
            
    ## Initialize subFamily as unknown
    percDiv = "NA"
    subFamily = "NA"

    ## Skip orphan transductions
    if (insertionType != "TD2"):
        insertionId = family + "_" + insertionType + "_" + chrom + "_" + pos

        ### Create output directory
        MEIDir = outDir + '/' + insertionId
        command = 'mkdir -p ' + MEIDir 
        os.system(command) 

        ### 2.1 Make list with the discordant read pair ids supporting the MEI
        discordantPlusIdList = VCFlineObj.infoDict["DP"].split(',')
        discordantMinusIdList = VCFlineObj.infoDict["DN"].split(',')
        discordantIdList = discordantPlusIdList + discordantMinusIdList

        ### 2.2 Make list with the mates of the anchor reads
        beg = int(pos) - 500
        end = int(pos) + 500
        anchorMateList = retrieveAnchorMates(chrom, beg, end, discordantIdList, bamFile)

        ### 2.3 Produce fasta with the anchor mate sequences supporting the MEI
        anchorMatesFastaObj = fasta()

        ## For each read pair
        for anchorMateId in anchorMateList:

            anchorMateSeq = fastaObj.fastaDict[anchorMateId]
            anchorMatesFastaObj.fastaDict[anchorMateId] = anchorMateSeq
        
        ## write fasta
        anchorMatesFastaPath = MEIDir + "/anchorMates.fa"
        anchorMatesFastaObj.write_fasta(anchorMatesFastaPath)

        ## A) L1
        if family == "L1":
    
            ### 2.4 Align the reads into the reference transposon sequence
            referenceFasta = refDir + '/' + family + '_consensus.fa'
            anchorMatesBamPath = MEIDir + "/anchorMates.bam"
            command = 'bwa mem ' + referenceFasta + ' ' + anchorMatesFastaPath + ' | samtools view -b - | samtools sort - > ' + anchorMatesBamPath
            print command
            os.system(command) 

            ### 2.5 Index the bam file
            command = 'samtools index ' + anchorMatesBamPath
            print command
            os.system(command) 

            ### 2.6 Call SNP and INDELS
            MEIvcf = MEIDir + "/MEI.vcf" 
            command = 'samtools mpileup -f ' + referenceFasta + ' -g ' + anchorMatesBamPath + ' | bcftools call -c - | bcftools filter -i \'QUAL>20 && DP>1\' > ' + MEIvcf
            print command
            os.system(command) 

            ### 2.7 Infer subfamily based on diagnostic nucleotides
            subFamily = subFamilyL1(MEIvcf)  

        ## B) Alu, SVA or ERVK
        elif (family == "Alu") or (family == "SVA") or (family == "ERVK"): 

            ### 2.4 Assemble the reads corresponding to the inserted element. 
            contigsFastaPath = MEIDir + '/contigs.fa'
            kmerLen='21'

            command = 'velveth ' + MEIDir + ' ' + kmerLen + ' -fasta -short ' + anchorMatesFastaPath
            print command
            os.system(command) 

            command = 'velvetg ' + MEIDir + ' -exp_cov auto -cov_cutoff auto' 
            print command
            os.system(command) 

            ## Do cleaning:
            sequences = MEIDir + '/Sequences '
            roadmaps = MEIDir + '/Roadmaps '
            pregraph = MEIDir + '/PreGraph '
            stats = MEIDir + '/stats.txt '
            lastGraph = MEIDir + '/LastGraph '
            graph2 = MEIDir + '/Graph2 '
            command = 'rm ' + sequences + roadmaps + pregraph + stats + lastGraph + graph2
            print command
            os.system(command) 

            ## 2.5 Run repeats masker if velvet assembled at least one contig
            if (os.path.isfile(contigsFastaPath)) and (os.path.getsize(contigsFastaPath) > 0):
            
                repeatsMaskerOut = MEIDir + '/contigs.fa.out'
                command = 'RepeatMasker -qq -dir ' + MEIDir + ' ' + contigsFastaPath
                print command
                os.system(command)

                ### 2.6 Extract subFamily from repeats masker output
                with open(repeatsMaskerOut) as f:
                    lines = f.readlines()[3:]
                    bestScore = 0

                    ## Analyze one repeatMasker hit per iteration. 
                    # Select hit with highest smith-waterman score as representative
                    for line in lines:
                        line = line.rstrip('\n')
                        line = line.split()
                        swScore = int(line[0])

                        ## Current hit better than previous
                        if (swScore > bestScore):
                        
                            bestScore = swScore
                            percDiv = line[1]
                            subFamily = line[9]

                    ## If sw-score lower than threshold set subFamily as unknown
                    if (bestScore == 0): 
                        percDiv = "NA"
                        subFamily = "NA"
            else:
                print "[WARNING] NO CONTIG ASSEMBLED. SKIP REPEATS MASKER"         

    
    ## Add subFamily and percentage divergence to consensus to the info field:
    VCFlineObj.infoDict["SUBFAMILY"] = subFamily
    VCFlineObj.infoDict["PDIV"] = percDiv
    VCFlineObj.info = VCFlineObj.make_info()

#### 4. Make output VCF
outFilePath = outDir + '/' + fileName + ".subfamily.vcf"

## 4.1 Write header
VCFObj.write_header(outFilePath)

## 4.2 Write variants
VCFObj.write_variants(outFilePath)

## End ##
print
print "***** Finished! *****"
print

