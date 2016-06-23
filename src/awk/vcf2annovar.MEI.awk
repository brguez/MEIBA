#!/usr/bin/env awk

# ******************************************************************************
        
#         vcf2annovar.MEI.awk
#         
#         Copyright (c) 2016 Bernardo Rodríguez-Martín
        
#         Mobile Genomes & Disease Lab.
#         Universidad de Vigo (Spain)

#         Licenced under the GNU General Public License 3.0 license.
# ******************************************************************************

# Description
################ 
		
# Takes as input a vcf containing TraFiC MEI and convert insertion calls into annovar input format

# Usage
#######

# awk -f vcf2annovar.MEI.awk input.vcf

#### Input (VCF with MEI calls):
##fileformat=VCFv4.2
##fileDate=20160620 
##source=TraFiCvX 
##reference=hs37d5
##contig=<ID=1,assembly=GRCh37,length=249250621,species=human>
##contig=<ID=2,assembly=GRCh37,length=243199373,species=human>
##contig=<ID=3,assembly=GRCh37,length=198022430,species=human>
# ........
# 8	135116616	.	A	<MEI>	.	.	SVTYPE=<MEI>;CLASS=L1;TYPE=TD0;SCORE=2A;CIPOS=0;STRAND=+;STRUCT=DEL;LEN=1139;CONTIGA=CACAGCACCCAGCTAGGTTTACTGGCGAGTTCTAACAATTTAAGGAAGAAATTATAACTTTCTTATGCATTTTTTAAAAATTAAAAATATCCAAACCTCCCTTTGTTTTCTGCATTTATCGTCTTTTTTCTAACTATTTGCATAAGATCCAAAATGTGCAAATAAGCCTTTGTTCTCTCTGAACACTATAAAACAAGGTTACTTTAGCTTACACCATGTGCATTTCTTTTGTTTTAAGGAAAAGAAATCTGTATCTTATACAACTCCAGGGTCCCTCAGTACTAAGTTTTCCTTTCCAATTAATAGCATGCCACAAGATGTTTGACTACTTTGCAGCTGATTTTGCATTCACTTACCTCTAGCCAATTGTTTTAATAATGATGCATCCTGAGTAGTGTCTTTCTCCATTTAGGCCTAAAGAAAACCTAGGCATTACCATTCAGGACATAGGCGTGGGCAAGGACTTCATGTCCAAAACACCAAAAGCAATGGCAACAAAAGCCAAAATTGACAAATGGGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTACCATCAGAGTGAACAGGCAACCTACAACATGGGAGAAAATTTTTGCAACCTACTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAACTCAAACAAATTTACAAGAAAAAAACAAACAACCCCATCAAAAAGTGGGCGAAGGACATGAACAGAcacttctcaa
# 8	140170469	.	G	<MEI>	.	.	SVTYPE=<MEI>;CLASS=L1;TYPE=TD0;SCORE=2A;CIPOS=0;STRAND=-;STRUCT=INV;CONTIGA=AGATTTGGTCTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGCTCATTTCTTTTTATTCTTTTTTCTCTAAACTTCCCTTCTCGCTTCATTTCATTCATTTCGTCTTCCATTGCTGATACCCTTTCTTCCAGTTGATCGCATCGGCTCCTGAGGCTTCTGCATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCTCCATCAGCCGCTTTTGAAAGCAAGATGGCATAACAAAGTGCCTCCTTAATGGGTTGTAAAGAAAAGAACCATATATTTCACTTGGTCCCCTTCATCCTTGTTCTTGACAACCAGTTATTCTGAAGTTGGTTTGTCTGAAGCTTACTCAGCACTTAAGTACCTGGATGATAAATGGAATGAAGACAAGACACAGCAGTTAAAGGAACCAGCAGAGATTAAGGAAGTGTTAGGCCAGAGTTCGTCTGGGGAAAAAATGGGCCCCAGGCTTGAGATTCCTCTGAGAAAAGGC

#### Output (Text file containing MEI calls in proper format for variant annotation with Annovar):
# 8	135116616	135116616	0	0
# 8	140170469	140170469	0	0

! /^#/ {
	# Get input 
	chrom = $1;
	beg = $2;
	end = $2;

	# Print row with updated template positions
	row=chrom"\t"beg"\t"end"\t0\t0";		
	print row;
}