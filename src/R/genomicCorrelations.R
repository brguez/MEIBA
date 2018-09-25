library(RColorBrewer)
library(MASS)

## Read input file
table <- read.table(file = '/Users/brodriguez/Research/Projects/Pancancer/Somatic/Analysis/TEIBA_0.8.3/genomicDistribution/supplementaryTable4.alignable.tsv', sep = '\t', header = TRUE)
table

#### 1. Number of L1 and replication time scatterplot 
######################################################
nbL1RT <- table[,c("nbL1","medianRT")]

dim(nbL1RT)
# [1] 3053    2

head(nbL1RT)
#nbL1 medianRT
#1    0 627.9740
#2    0 106.7028
#3    0 141.8966
#4    0 414.5490
#5    5 798.0827
#6    7 819.6064

# Remove rows with replication time NA values:
nbL1RT <- nbL1RT[complete.cases(nbL1RT), ]

dim(nbL1RT)
# [1] 2734    2

#### Make the scatterplot:
#Outfile:
outDir <- "/Users/brodriguez/Research/Projects/Pancancer/Somatic/Analysis/TEIBA_0.8.3/genomicDistribution/Pictures/"
outFile <- paste(outDir, "nbL1_medianRT_corr4.pdf", sep="")
pdf(outFile, width=4, height=4)

# A color palette from blue to yellow to red
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(nbL1RT$nbL1, nbL1RT$medianRT, n=50)

# Make the base plot
plot(nbL1RT, xlab="L1 insertions / Mb", ylab="Replication time", ylim=c(0,85), pch=19, cex=.4)

# Draw the colored contour lines
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

dev.off()


#### 2. Number of L1 and L1 EN motif
######################################
#Read input data
nbL1Motif <- table[,c("nbL1","nbL1Motif")]

dim(nbL1Motif)
# [1] 3053    2

head(nbL1Motif)
#nbL1 nbL1Motif
#1    0       0.126
#2    4       0.275
#3    1       0.578
#4    3       0.217
#5    8       0.000
#6   10       0.241

##Remove bins with number motif == 0
# These will mostly  correspond to telomeric, centromeric regions
nbL1Motif <- nbL1Motif[(nbL1Motif$nbL1Motif != 0), ]
dim(nbL1Motif)
#[1] 2888    2

#### Make the scatterplot:
#Outfile:
outDir <- "/Users/brodriguez/Research/Projects/Pancancer/Somatic/Analysis/TEIBA_0.8.3/genomicDistribution/Pictures/"
outFile <- paste(outDir, "nbL1_nbL1Motif_corr3.pdf", sep="")
pdf(outFile, width=4, height=4)

# A color palette from blue to yellow to red
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(nbL1Motif$nbL1, nbL1Motif$nbL1Motif, n=50)

# Make the base plot
plot(nbL1Motif, xlab="L1 insertions / Mb", ylim=c(1,16800), ylab="L1 EN Motifs / Mb", pch=19, cex=.4)

# Draw the colored contour lines
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

dev.off()
