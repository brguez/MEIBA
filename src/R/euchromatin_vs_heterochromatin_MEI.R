# Inigo Martincorena - 24.11.2016
# Relative density of transpositions in different colours of chromatin

infile = "/Users/brodriguez/Research/Projects/Pancancer/Somatic/Analysis/TEIBA_0.5.7/GenomicDistribution/pos_class_germline_MEI.txt"

# Libraries
library("GenomicRanges")

# Function to compare TE insertion rates and generate each plot
get_TEins_rates = function(te_gr,filename) {

    bed_file_list = c(
R="/Users/brodriguez/Research/References/Annotations/H.sapiens/hg19/Chromatin/R_intersect.bed",
T="/Users/brodriguez/Research/References/Annotations/H.sapiens/hg19/Chromatin/T_intersect.bed")

    rates = data.frame(chrom=names(bed_file_list), L=NA, n=NA, rateperMb=NA, cilow=NA, cihigh=NA)
    
    for (j in 1:length(bed_file_list)) {
    
        bed = read.table(bed_file_list[j], header=0, sep="\t", stringsAsFactors=F)
        bed[,1] = substr(bed[,1], 4, nchar(bed[,1]))
        bed_gr = reduce(GRanges(bed[,1], IRanges(bed[,2], bed[,3]))) # GRanges object
    
        # Mapping TE insertions to the bed regions
        ol = as.matrix(findOverlaps(te_gr, bed_gr, type="any", select="all"))
        rates[j,2] = sum(end(bed_gr)-start(bed_gr)+1)
        rates[j,3] = length(unique(ol[,1]))
        rates[j,4] = rates$n[j] / rates$L[j] * 1e6
        rates[j,5:6] = poisson.test(rates$n[j])$conf.int / rates$L[j] * 1e6
    }
    pval = poisson.test(x=rates$n, T=rates$L)$p.value

    # Barplot
    dev.new(width=3,height=5)
    bc = barplot(rates$rateperMb, names.arg=names(bed_file_list), col=c("gray","darkolivegreen4"), border=NA, 
        las=1, ylim=c(0,max(rates$cihigh)), ylab="TE insertions per Mb", main=sprintf("%s - pval=%0.4g", substr(filename,1,nchar(filename)-4),pval), cex.main=0.7)
    arrows(bc, rates$cilow, bc, rates$cihigh, lwd=1.5, angle=90, code=3, length=0)
    dev.copy(pdf,file=filename, width=3, height=5)
    dev.off()

    return(rates)
}


# Loading the input data
te_input = read.table(infile, header=0, sep="\t", stringsAsFactors=F)
te = data.frame(chr=sapply(strsplit(te_input[,1],split=":"), function(x) x[1]),
pos=as.numeric(sapply(strsplit(te_input[,1],split=":"), function(x) x[2])), type=te_input[,2], stringsAsFactors=F)
te_gr = GRanges(te[,1], IRanges(te[,2],te[,2]))

# 1. All TEs
r_all = get_TEins_rates(te_gr,filename="chromatin_allTEs.pdf")

# 2. By type
r_Alu = get_TEins_rates(te_gr[which(te$type=="Alu")],filename="chromatin_Alu.pdf")
r_L1 = get_TEins_rates(te_gr[which(te$type=="L1")],filename="chromatin_L1.pdf")
r_SVA = get_TEins_rates(te_gr[which(te$type=="SVA")],filename="chromatin_SVA.pdf")


