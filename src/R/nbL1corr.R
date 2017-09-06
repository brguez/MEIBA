
library(scales)

## Read input file
table <- read.table(file = '/Users/brodriguez/Research/Projects/Pancancer/Germline/Analysis/TEIBA_0.6.5/sourceElements/activity/donorId_nbL1_nbActive.tsv', sep = '\t', header = TRUE)
head(table)

## Remove donors with 0 insertions
# These will mostly  correspond to telomeric, centromeric regions

#### Make the scatterplot:

#Outfile:
outDir <- "/Users/brodriguez/Research/Projects/Pancancer/Germline/Analysis/TEIBA_0.6.5/sourceElements/activity/Pictures/"
outFile <- paste(outDir, "nbL1_nbActiveSrcElements_correlation.pdf", sep="")
pdf(outFile, width=4, height=4)

# Make the base plot
plot(table$nbL1, table$nbActive, xlab="Number L1 insertions / Mb", ylab="Number active source elements", pch=19, cex=.4)
# xlim=c(0,1200), 

##### Assess correlation:
cor.test(table$nbL1, table$nbActive, method="spearman")
dev.off()
