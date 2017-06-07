output.file = "/Users/brodriguez/Research/Projects/Pancancer/Germline/Analysis/TEIBA_0.5.7/SourceElements/MainFigure/PanelA/Pictures/source_elements_chromosomes_plot_legend.pdf" 
pdf(output.file, 10, 8)

POS.A <- c(1,2,2,3)
hist(POS.A)

legend("topleft", title="VAF(%)", legend=c("0-25", ">25-50", ">50-75", ">75-100"), col="black", pch=25, pt.cex=c(1,1.5,2,3),bty="n")
legend("topright", title="Activity", legend=c("0", ">0-3", ">3-9", ">9"), col="black", pch=25, pt.bg=c("white","yellow","orange","red"), bty="n")

dev.off()  # close PDF
