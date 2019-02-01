source("~/MAGIC/Colour_Scheme.R")
## Across Dataset Consistency ##
dat <- readRDS("TM_Data_Consistency.rds")
dat <- dat[dat$auc != 0,]
dat[is.na(dat)] <- 1;

# Average over all tissues just consistency
meths <- unique(dat$method)
aucs <- sort(unique(dat$auc))
OUT <- matrix(0, ncol=length(meths), nrow=length(aucs))
bothN_OUT <- matrix(0, ncol=length(meths), nrow=length(aucs))
dropN_OUT <- matrix(0, ncol=length(meths), nrow=length(aucs))
facsN_OUT <- matrix(0, ncol=length(meths), nrow=length(aucs))
for (m in 1:length(meths)) {
	for (t in unique(dat$tissue)) {
		OUT[,m] <- OUT[,m]+dat[dat$method==meths[m] & dat$tissue==t,"percagree"]
		bothN_OUT[,m] <- bothN_OUT[,m]+dat[dat$method==meths[m] & dat$tissue==t,"nboth"]
		dropN_OUT[,m] <- dropN_OUT[,m]+dat[dat$method==meths[m] & dat$tissue==t,"ndrop"]
		facsN_OUT[,m] <- facsN_OUT[,m]+dat[dat$method==meths[m] & dat$tissue==t,"nfacs"]
	}
}
OUT <- OUT/length(unique(dat$tissue))
rownames(OUT) <- aucs;
colnames(OUT) <- remap_names(meths);

bothN_OUT <- bothN_OUT/length(unique(dat$tissue))
rownames(bothN_OUT) <- aucs;
colnames(bothN_OUT) <- remap_names(meths);

dropN_OUT <- dropN_OUT/length(unique(dat$tissue))
rownames(dropN_OUT) <- aucs;
colnames(dropN_OUT) <- remap_names(meths);

facsN_OUT <- facsN_OUT/length(unique(dat$tissue))
rownames(facsN_OUT) <- aucs;
colnames(facsN_OUT) <- remap_names(meths);

OUT <- OUT[, match(names(Method_Legend), colnames(OUT))]
bothN_OUT <- bothN_OUT[, match(names(Method_Legend), colnames(bothN_OUT))]
facsN_OUT <- facsN_OUT[, match(names(Method_Legend), colnames(facsN_OUT))]
dropN_OUT <- dropN_OUT[, match(names(Method_Legend), colnames(dropN_OUT))]

png("Figure3_Markers_vs_AUC.png", width=8, height=8, units="in", res=300)
layout(rbind(c(1,2,3,3), c(4,5,6,6), c(7,8,8,9)), widths=c(2,2,1,1), heights=c(1,1,1))
plotlabs <- c("A", "A", "B", "C", "D", "E", "F", "G")
for ( j in 1:ncol(OUT) ) {
	if (colnames(OUT)[j] == "logcounts") {next;}
	par(new=FALSE)
	par(mar=c(4,4,3,4))
	coords <- barplot(bothN_OUT[,j], names=rownames(bothN_OUT), 
		xlab="AUC", ylab="Genes (#)", main=colnames(bothN_OUT)[j], 
		xlim=c(0.25,nrow(bothN_OUT)+1.5), ylim=c(0,max(bothN_OUT)))
	mtext(plotlabs[j], at=8500, side=2, line=2.25, font=2, las=2)
	par(new=TRUE)
	plot(coords, OUT[,j], new=FALSE, xaxt="n", yaxt="n", 
		main="", xlab="", ylab="", 
		xlim=c(0.25,nrow(bothN_OUT)+1.5), ylim=c(0.5,1), type="b")
	axis(side=4, at=seq(from=0.5, to=1, by=0.1), 
		     labels=seq(from=0.5, to=1, by=0.1)*100)
	mtext("Reproducibility (%)", side=4, line=2.25, at=0.75, cex=par()$cex)
	abline(h=0.95, lty=3)

}


this_row="0.8"
tots <- facsN_OUT[this_row,] + dropN_OUT[this_row,] - bothN_OUT[this_row,]
prop_shared <- cbind(facsN_OUT[this_row,], dropN_OUT[this_row,], bothN_OUT[this_row,])
prop_shared[,1] <- prop_shared[,1]-prop_shared[,3]
prop_shared[,2] <- prop_shared[,2]-prop_shared[,3]
prop_shared <- prop_shared/tots;
prop_shared <- prop_shared[rownames(prop_shared) != "logcounts",]

par(mar=c(4,5,1,1))
coords <- barplot(t(prop_shared)*100, col=c("cornflowerblue", "goldenrod", "grey35"), horiz=TRUE, xlab="Markers (%)", las=1)
mtext("H", at=max(coords)+0.5, side=2, line=3.25, font=2, las=2)

blank_plot <- function() {
	tmp <-  par("mar")
	par(mar=c(0,0,0,0))
	plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
	return(tmp);
}
blank_plot()
legend("left", c("SS2", "10X", "Both"), fill=c("cornflowerblue", "goldenrod", "grey35"), bty="n")

dev.off();

## Heatmap ##
require("gplots")
require("RColorBrewer")
MvM <- readRDS("TM_Method_Consistency.rds")
MvM <- round((1-MvM)*100)
MvMtext <- matrix(paste(MvM, "%", sep=""), ncol=ncol(MvM))
png("Method_vs_Method_heatmap.png", width=5, height=5, units="in", res=300)
heatmap.2(MvM, trace="none", cellnote=MvMtext, col=brewer.pal(8, "Reds"), notecol="black", margins=c(6,6), key=FALSE, dendrogram="none")
dev.off()

## pval cors ##
require("gplots")
require("RColorBrewer")
p_val_cors <- readRDS("TM_mark_cors.rds")

png("Pvalue_Cor_heatmap.png", width=5, height=5, units="in", res=300)
heatmap.2(abs(t(p_val_cors)), trace="none", cellnote=t(round(p_val_cors, digits=2)), col=brewer.pal(8, "Blues"), notecol="black", margins=c(6,6), key=FALSE, dendrogram="none", breaks=9)
dev.off()




### Absolute Reproducible Markers ###
png("TM_Abs_Repro_Markers.png", width=6, height=6, units="in", res=300)
matplot(rownames(OUT), round(OUT*bothN_OUT), type="l", lty=1, col=Method_Legend, xlab="AUC Threshold", ylab="Total Reproducible Markers", lwd=2.5)
legend("topright", names(Method_Legend), lwd=2.5, col=Method_Legend, bty="n")
dev.off()

### Correlate p-values ###
