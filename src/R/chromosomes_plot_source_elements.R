
##### CUSTOM PLOTTING FUNCTION: my.prepareGenomePlot #########################
my.prepareGenomePlot = function (chrompos = NULL, cols = "grey50", paintCytobands = FALSE, 
                                 bleach = 0, topspace = 1, organism, sexChromosomes = FALSE, 
                                 units = "hg19", cytobandWidth = 0.075, ...) 
{
  par(mar = c(1, 4, 2, 3) + 0.1)
  if (!missing(organism)) {
    organism <- match.arg(organism, c("hsa", "mmu", "rno"))
    chrom.n <- switch(organism, hsa = 22, mmu = 19, rno = 20)
    if (is.null(chrompos)) {
      chroms = c(1:chrom.n, if (sexChromosomes) c(98, 99) else NULL)
      chrnames = characterCHR(chroms, "chr")
      mapinfo = rep(0, length(chroms))
      chrompos = cbind(CHR = chroms, MapInfo = mapinfo)
      rownames(chrompos) <- chrnames
    }
    chrs2 <- factor(numericCHR(chrompos[, "CHR"]), levels = c(1:chrom.n, 
                                                              if (sexChromosomes) c(98, 99) else NULL))
    if (organism %in% c("hsa", "mmu")) 
      lens <- lengthChromosome(levels(chrs2), units = units)
    else lens <- sapply(split(chrompos[, "MapInfo"], chrs2), 
                        function(x) max(c(0, x)))
    names(lens) <- characterCHR(names(lens))
    cols <- rep(cols, length.out = length(lens))
    names(cols) <- names(lens)
    dwidth <- NULL
    for (i in 1:(chrom.n%/%2)) dwidth[i] <- lens[i] + lens[chrom.n + 
                                                             1 - i]
    if (chrom.n%%2 == 1) 
      dwidth <- c(dwidth, lens[chrom.n%/%2 + 1])
    if (sexChromosomes) 
      dwidth <- c(dwidth, lens["X"] + lens["Y"])
    maxdwidth <- max(dwidth) * 1.05
    leftrow <- c(if (sexChromosomes) "X" else NULL, ((chrom.n + 
                                                        1)%/%2):1)
    rightrow <- c(if (sexChromosomes) "Y" else NULL, if (chrom.n%%2 == 
                                                         1) "" else NULL, ((chrom.n + 1)%/%2 + 1):chrom.n)
    plot(c(0, maxdwidth), c(0.5, 0.5 + length(dwidth) + topspace), 
         type = "n", ylab = "Chromosome", xlab = "", axes = FALSE, 
         las = 2, ...)
    axis(2, c(1:length(dwidth)), characterCHR(leftrow), las = 2)
    axis(4, c(1:length(dwidth)), characterCHR(rightrow), 
         las = 2)
    if (paintCytobands && organism %in% c("hsa", "mmu")) {
      for (i in 1:length(dwidth)) {
        if (lens[leftrow[i]] > 0) 
          paintCytobands(leftrow[i], c(0, i + cytobandWidth/2), 
                            units = units, width = cytobandWidth, length.out = lens[leftrow[i]], 
                            legend = FALSE, bleach = bleach)
        if (rightrow[i] != "" && lens[rightrow[i]] > 0) 
          paintCytobands(rightrow[i], c(maxdwidth - lens[rightrow[i]], 
                            i + cytobandWidth/2), units = units, width = cytobandWidth, 
                            length.out = lens[rightrow[i]], legend = FALSE, bleach = bleach)
      }
    }
    else {
      for (i in 1:length(dwidth)) {
        lines(c(0, lens[leftrow[i]]), c(i, i), col = cols[leftrow[i]], 
              lwd = 2)
        if (rightrow[i] != "") 
          lines(c(maxdwidth - lens[rightrow[i]], maxdwidth), 
                c(i, i), col = cols[rightrow[i]], lwd = 2)
      }
    }
    dchrompos <- matrix(0, nrow = nrow(chrompos), ncol = 2, 
                        dimnames = list(rownames(chrompos), c("CHR", "MapInfo")))
    for (i in 1:length(rightrow)) if (rightrow[i] != "") {
      probes <- characterCHR(chrompos[, "CHR"]) == rightrow[i]
      dchrompos[probes, 2] <- chrompos[probes, "MapInfo"] + 
        maxdwidth - lens[rightrow[i]]
      dchrompos[probes, 1] <- i
    }
    for (i in 1:length(leftrow)) {
      probes <- characterCHR(chrompos[, "CHR"]) == leftrow[i]
      dchrompos[probes, 2] <- chrompos[probes, "MapInfo"]
      dchrompos[probes, 1] <- i
    }
  }
  else {
    chrs2 <- factor(numericCHR(chrompos[, "CHR"]))
    lens <- sapply(split(chrompos[, "MapInfo"], chrs2), max)
    m <- length(lens)
    cols <- rep(cols, length.out = m)
    maxdwidth <- max(lens)
    plot(c(0, maxdwidth), c(0.5, m + 0.5 + topspace), type = "n", 
         ylab = "Chromosome", xlab = "", axes = FALSE, las = 2, 
         ...)
    axis(2, c(m:1), characterCHR(names(lens)), las = 2)
    for (i in 1:m) lines(c(0, lens[i]), c(m + 1 - i, m + 
                                            1 - i), col = cols[as.numeric(names(lens))], lwd = 2)
    dchrompos <- chrompos
    dchrompos[, 1] <- m + 1 - as.numeric(chrs2)
  }
  dchrompos
}

# 1) Input genomic positions


# 2) Plot ideograms
output.file = "/Users/brodriguez/Research/Projects/Pancancer/Germline/Analysis/TEIBA_0.5.7/SourceElements/ChromPlot/Pictures/source_elements_chromosomes_plot.pdf"  
pdf(output.file, 10, 8)  # open PDF
par(mgp=c(2.5,1,0))

library(quantsmooth)
#data = data.frame(CHR=c(CHR.A), MapInfo=c(POS.A), FreqInfo=c(FREQ.A), ActInfo=c(ACT.A))
# add CHR.C, etc. if needed
data = read.table("/Users/brodriguez/Research/Projects/Pancancer/Germline/Analysis/TEIBA_0.5.7/SourceElements/ChromPlot/germline_source_elements_input_chromPlot.tsv", header=T, sep="\t")


# 3)
# Plot positions
colours = c("firebrick", "navyblue", "forestgreen", "orange")  
# Type 1 New (red)
data1.1 = data.frame(CHR=c(), MapInfo=c())
data1.2 = data.frame(CHR=c(), MapInfo=c())
data1.3 = data.frame(CHR=c(), MapInfo=c())
data1.4 = data.frame(CHR=c(), MapInfo=c())

data2.1 = data.frame(CHR=c(), MapInfo=c())
data2.2 = data.frame(CHR=c(), MapInfo=c())
data2.3 = data.frame(CHR=c(), MapInfo=c())
data2.4 = data.frame(CHR=c(), MapInfo=c())

data3.1 = data.frame(CHR=c(), MapInfo=c())
data3.2 = data.frame(CHR=c(), MapInfo=c())
data3.3 = data.frame(CHR=c(), MapInfo=c())
data3.4 = data.frame(CHR=c(), MapInfo=c())

data4.1 = data.frame(CHR=c(), MapInfo=c())
data4.2 = data.frame(CHR=c(), MapInfo=c())
data4.3 = data.frame(CHR=c(), MapInfo=c())
data4.4 = data.frame(CHR=c(), MapInfo=c())

for(i in 1:nrow(data)){
    if (data$FreqInfo[i] > 0 & data$FreqInfo[i] <= 0.25) {
      if (data$ActInfo[i] == 0) {
        data1.1 <- rbind( data1.1, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 0 & data$ActInfo[i] <= 3) {
        data1.2 <- rbind( data1.2, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 3 & data$ActInfo[i] <= 9) {
        data1.3 <- rbind( data1.3, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else #(data$FreqInfo[i] == 4)
        data1.4 <- rbind( data1.4, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
    } else if (data$FreqInfo[i] > 0.25 & data$FreqInfo[i] <= 0.5) {
      if (data$ActInfo[i] == 0) {
        data2.1 <- rbind( data2.1, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 0 & data$ActInfo[i] <= 3) {
        data2.2 <- rbind( data2.2, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 3 & data$ActInfo[i] <= 9) {
        data2.3 <- rbind( data2.3, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else #(data$FreqInfo[i] == 4)
        data2.4 <- rbind( data2.4, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
    } else if (data$FreqInfo[i] > 0.5 & data$FreqInfo[i] <= 0.75) {
      if (data$ActInfo[i] == 0) {
        data3.1 <- rbind( data3.1, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 0 & data$ActInfo[i] <= 3) {
        data3.2 <- rbind( data3.2, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 3 & data$ActInfo[i] <= 9) {
        data3.3 <- rbind( data3.3, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else #(data$FreqInfo[i] == 4)
        data3.4 <- rbind( data3.4, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
    } else #(data$FreqInfo[i] == 4)
      if (data$ActInfo[i] == 0) {
        data4.1 <- rbind( data4.1, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 0 & data$ActInfo[i] <= 3) {
        data4.2 <- rbind( data4.2, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else if (data$ActInfo[i] > 3 & data$ActInfo[i] <= 9) {
        data4.3 <- rbind( data4.3, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
      } else #(data$FreqInfo[i] == 4)
        data4.4 <- rbind( data4.4, data.frame("CHR"=data$CHR[i], "MapInfo"=data$MapInfo[i]))
}

start.1.1 = 1            # First row of this type in the data table (see below)
if (length(data1.1) != 0){
  end.1.1 = start.1.1 + length(data1.1$CHR) - 1 
} else {
  end.1.1 = start.1.1
}
if (length(data1.2) != 0){
  start.1.2 = end.1.1 + 1          # First row of this type in the data table (see below)
  end.1.2 = start.1.2 + length(data1.2$CHR) - 1
} else{
  start.1.2 = start.1.1
  end.1.2 = end.1.1
}
if (length(data1.3) != 0){
  start.1.3 = end.1.2 + 1          # First row of this type in the data table (see below)
  end.1.3 = start.1.3 + length(data1.3$CHR) - 1
} else{
  start.1.3 = start.1.2
  end.1.3 = end.1.2
}
if (length(data1.4) != 0){
  start.1.4 = end.1.3 + 1          
  end.1.4 = start.1.4 + length(data1.4$CHR) - 1
} else {
  start.1.4 = start.1.3
  end.1.4 = end.1.3
}

if (length(data2.1) != 0){
  start.2.1 = end.1.4 + 1           # First row of this type in the data table (see below)
  end.2.1 = start.2.1 + length(data2.1$CHR) - 1
} else{
  start.2.1 = start.1.4
  end.2.1 = end.1.4
}
if (length(data2.2) != 0){
start.2.2 = end.2.1 + 1          # First row of this type in the data table (see below)
end.2.2 = start.2.2 + length(data2.2$CHR) - 1
} else {
  start.2.2 = start.2.1
  end.2.2 = end.2.1
}
if (length(data2.3) != 0){
  start.2.3 = end.2.2 + 1          # First row of this type in the data table (see below)
  end.2.3 = start.2.3 + length(data2.3$CHR) - 1
} else{
  start.2.3 = start.2.2
  end.2.3 = end.2.2
}
if (length(data2.4) != 0){
  start.2.4 = end.2.3 + 1          
  end.2.4 = start.2.4 + length(data2.4$CHR) - 1
} else{
  start.2.4 = start.2.3
  end.2.4 = end.2.3
}

if (length(data3.1) != 0){
start.3.1 = end.2.4 + 1            # First row of this type in the data table (see below)
end.3.1 = start.3.1 + length(data3.1$CHR) - 1
} else{
  start.3.1 = start.2.4
  end.3.1 = end.2.4
}
if (length(data3.2) != 0){
start.3.2 = end.3.1 + 1          # First row of this type in the data table (see below)
end.3.2 = start.3.2 + length(data3.2$CHR) - 1
} else {
  start.3.2 = start.3.1
  end.3.2 = end.3.1
}
if (length(data3.3) != 0){
start.3.3 = end.3.2 + 1          # First row of this type in the data table (see below)
end.3.3 = start.3.3 + length(data3.3$CHR) - 1
} else{
  start.3.3 = start.3.2
  end.3.3 = end.3.2
}
if (length(data3.4) != 0){
start.3.4 = end.3.3 + 1          
end.3.4 = start.3.4 + length(data3.4$CHR) - 1
} else{
  start.3.4 = start.3.3
  end.3.4 = end.3.3
}

if (length(data4.1) != 0){
start.4.1 = end.3.4 + 1            # First row of this type in the data table (see below)
end.4.1 = start.4.1 + length(data4.1$CHR) - 1
} else {
  start.4.1 = start.3.4
  end.4.1 = end.3.4
}
if (length(data4.2) != 0){
start.4.2 = end.4.1 + 1          # First row of this type in the data table (see below)
end.4.2 = start.4.2 + length(data4.2$CHR) - 1
} else{
  start.4.2 = start.4.1
  end.4.2 = end.4.1
}
if (length(data4.3) != 0){
start.4.3 = end.4.2 + 1          # First row of this type in the data table (see below)
end.4.3 = start.4.3 + length(data4.3$CHR) - 1
} else{
  start.4.3 = start.4.2
  end.4.3 = end.4.2
}
if (length(data4.4) != 0){
start.4.4 = end.4.3 + 1          
end.4.4 = start.4.4 + length(data4.4$CHR) - 1
} else{
  start.4.4 = start.4.3
  end.4.4 = end.4.3
}

dataall <- rbind(data1.1, data1.2, data1.3, data1.4, data2.1, data2.2, data2.3, data2.4, data3.1, data3.2, data3.3, data3.4, data4.1, data4.2, data4.3, data4.4)

chrom.pos = my.prepareGenomePlot(dataall, paintCytobands=T, organism="hsa", sexChromosomes=T,  # for plotting chrX and chrY
                                 cytobandWidth=0.2)


# Plot title
#title("The topography of germline L1 source elements in cancer", line=0, cex=1.1)

# Novel: Freq=4 Act=4
points(chrom.pos[start.4.4:end.4.4, 2], chrom.pos[start.4.4:end.4.4, 1] + 0.45, pch=25, bg="red", cex=3)
# Novel: Freq=4 Act=3
points(chrom.pos[start.4.3:end.4.3, 2], chrom.pos[start.4.3:end.4.3, 1] + 0.45, pch=25, bg="#fd8d3c", cex=3)
# Novel: Freq=4 Act=2
points(chrom.pos[start.4.2:end.4.2, 2], chrom.pos[start.4.2:end.4.2, 1] + 0.45, pch=25, bg="#ffffb2", cex=3)
# Novel: Freq=4 Act=1
points(chrom.pos[start.4.1:end.4.1, 2], chrom.pos[start.4.1:end.4.1, 1] + 0.45, pch=25, bg="white", cex=3)

# Novel: Freq=3 Act=4
points(chrom.pos[start.3.4:end.3.4, 2], chrom.pos[start.3.4:end.3.4, 1] + 0.35, pch=25, bg="red", cex=2)
# Novel: Freq=3 Act=3
points(chrom.pos[start.3.3:end.3.3, 2], chrom.pos[start.3.3:end.3.3, 1] + 0.35, pch=25, bg="#fd8d3c", cex=2)
# Novel: Freq=3 Act=2
points(chrom.pos[start.3.2:end.3.2, 2], chrom.pos[start.3.2:end.3.2, 1] + 0.35, pch=25, bg="#ffffb2", cex=2)
# Novel: Freq=3 Act=1
points(chrom.pos[start.3.1:end.3.1, 2], chrom.pos[start.3.1:end.3.1, 1] + 0.35, pch=25, bg="white", cex=2)

# Novel: Freq=2 Act=4
points(chrom.pos[start.2.4:end.2.4, 2], chrom.pos[start.2.4:end.2.4, 1] + 0.29, pch=25, bg="red", cex=1.5)
# Novel: Freq=2 Act=3
points(chrom.pos[start.2.3:end.2.3, 2], chrom.pos[start.2.3:end.2.3, 1] + 0.29, pch=25, bg="#fd8d3c", cex=1.5)
# Novel: Freq=2 Act=2
points(chrom.pos[start.2.2:end.2.2, 2], chrom.pos[start.2.2:end.2.2, 1] + 0.29, pch=25, bg="#ffffb2", cex=1.5)
# Novel: Freq=2 Act=1
points(chrom.pos[start.2.1:end.2.1, 2], chrom.pos[start.2.1:end.2.1, 1] + 0.29, pch=25, bg="white", cex=1.5)

# Novel: Freq=1 Act=1
points(chrom.pos[start.1.1:end.1.1, 2], chrom.pos[start.1.1:end.1.1, 1] + 0.22, pch=25, bg="white")
# Novel: Freq=1 Act=2
points(chrom.pos[start.1.2:end.1.2, 2], chrom.pos[start.1.2:end.1.2, 1] + 0.22, pch=25, bg="#ffffb2")
# Novel: Freq=1 Act=3
points(chrom.pos[start.1.3:end.1.3, 2], chrom.pos[start.1.3:end.1.3, 1] + 0.22, pch=25, bg="#fd8d3c")
# Novel: Freq=1 Act=4
points(chrom.pos[start.1.4:end.1.4, 2], chrom.pos[start.1.4:end.1.4, 1] + 0.22, pch=25, bg="red")








# 4)
# Add legend
#legend("topleft", title="VAF", legend=c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"), col="black", pch=25, pt.cex=c(1,1.5,2,3),bty="n")
#legend("topright", title="Activity", legend=c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"), col="black", pch=25, pt.bg=c("white","yellow","orange","red"), bty="n", cex=0.8)

dev.off()  # close PDF

