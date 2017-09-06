
## Read input file
table <- read.table(file = '/Users/brodriguez/Research/References/Annotations/H.sapiens/hg19/GenomicFeatures/supplementaryTable2.NA.tsv', sep = '\t', header = TRUE)


#### 1. Number of L1 and replication time scatterplot 
######################################################
nbL1RT <- table[,c("nbL1","medianRT")]

dim(nbL1RT)
# [1] 3053    2

head(nbL1RT)
#nbL1 medianRT
#1    0       NA
#2    4       NA
#3    1       NA
#4    3       NA
#5    8       NA
#6   10       NA

# Remove rows with replication time NA values:
nbL1RT <- nbL1RT[complete.cases(nbL1RT), ]

dim(nbL1RT)
# [1] 2734    2

#### Make the scatterplot:
#smoothScatter(nbL1RT)

#Outfile:
outDir <- "/Users/brodriguez/Research/Projects/Pancancer/Somatic/Analysis/TEIBA_0.6.5/genomicDistribution/Pictures/"
outFile <- paste(outDir, "nbL1_medianRT_correlation2.pdf", sep="")
pdf(outFile, width=4, height=4)

# A color palette from blue to yellow to red
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(nbL1RT$nbL1, nbL1RT$medianRT, n=50)

# Make the base plot
plot(nbL1RT, xlab="Number L1 insertions / Mb", ylab="Replication time", ylim=c(0,1200), pch=19, cex=.4)

# Draw the colored contour lines
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

dev.off()

#### 2. Number of L1 and median expression scatterplot 
#######################################################

#Read input data
nbL1Expr <- table[,c("nbL1","medianExpr")]

dim(nbL1Expr)
# [1] 3053    2

head(nbL1Expr)
#nbL1 medianExpr
#1    0       4664
#2    4      39709
#3    1     158014
#4    3      54272
#5    8        753
#6   10        947

# Remove rows with replication time NA values:
nbL1Expr <- nbL1Expr[complete.cases(nbL1Expr), ]
dim(nbL1Expr)

##Remove bins with expression level == 0
# These will mostly  correspond to telomeric, centromeric regions
nbL1Expr <- nbL1Expr[(nbL1Expr$medianExpr != 0), ]
dim(nbL1Expr)
#[1] 2867    2

##Apply log10 to the expression
nbL1Expr$medianLog10Expr <- log10(nbL1Expr$medianExpr)

nbL1Expr

#### Make the scatterplot:
#Outfile:
outDir <- "/Users/brodriguez/Research/Projects/Pancancer/Somatic/Analysis/TEIBA_0.6.5/genomicDistribution/Pictures/"
outFile <- paste(outDir, "nbL1_medianExpression_correlation2.pdf", sep="")
pdf(outFile, width=4, height=4)

#smoothScatter(nbL1RT)
# A color palette from blue to yellow to red
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(nbL1Expr$nbL1, nbL1Expr$medianLog10Expr, n=50)

# Make the base plot
plot(nbL1Expr$nbL1, nbL1Expr$medianLog10Expr, ylim=c(1,7), xlab="Number L1 insertions / Mb", ylab="log10 (Expression level)", pch=19, cex=.4)

# Draw the colored contour lines
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

dev.off()


#### 3. Number of L1 and replication time scatterplot 
######################################################
#Read input data
nbL1GnDensity <- table[,c("nbL1","geneDensity")]

dim(nbL1GnDensity)
# [1] 3053    2

head(nbL1GnDensity)
#nbL1 geneDensity
#1    0       0.126
#2    4       0.275
#3    1       0.578
#4    3       0.217
#5    8       0.000
#6   10       0.241

# Remove rows with gene density NA values:
nbL1GnDensity <- nbL1GnDensity[complete.cases(nbL1GnDensity), ]
dim(nbL1GnDensity)
# [1] 2734    2

##Remove bins with gene density == 0
# These will mostly  correspond to telomeric, centromeric regions
nbL1GnDensity <- nbL1GnDensity[(nbL1GnDensity$geneDensity != 0), ]
dim(nbL1GnDensity)
#[1] 2634    2

#### Make the scatterplot:
#smoothScatter(nbL1RT)

#Outfile:
outDir <- "/Users/brodriguez/Research/Projects/Pancancer/Somatic/Analysis/TEIBA_0.6.5/genomicDistribution/Pictures/"
outFile <- paste(outDir, "nbL1_geneDensity_correlation2.pdf", sep="")
pdf(outFile, width=4, height=4)

# A color palette from blue to yellow to red
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(nbL1GnDensity$nbL1, nbL1GnDensity$geneDensity, n=50)

# Make the base plot
plot(nbL1GnDensity, xlab="Number L1 insertions / Mb", ylab="Gene density", pch=19, cex=.4)

# Draw the colored contour lines
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

dev.off()
