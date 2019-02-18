# False Positive Distributions:
source("~/MAGIC/Colour_Scheme.R")
dat <- readRDS("Permuted_Fits.rds")
scores <- vector()
for (m in names(dat)) {
	score <- sum(dat[[m]][,3])/sum(dat[[m]])
	scores <- c(scores, score);
}
new_names <- Method_Colours[match(names(dat), Method_Colours[,1]),]
new_names <- names(Method_Legend)[match(new_names[,2],Method_Legend)]
new_names[is.na(new_names)] <- names(dat)[is.na(new_names)]

names(scores) <- remap_names(names(dat));
png("Supplementary_FP_Distribution.png", width=5, height=5, units="in", res=300)
barplot(scores,las=2, col=c("black", Method_Legend), ylab="Fit ZINB")
dev.off()

# False Positives by Method & Dataset
source("~/MAGIC/Colour_Scheme.R")
dat <- readRDS("Permuted_FPR.rds")
png("Figure3_Permuted_FPR.png", width=5, height=5, units="in", res=300)
layout(rbind(c(1,1), c(2,3)), widths=c(6,2), heights=c(1,1))

# FACS
dat_FACS <- unlist( dat[grep("FACS", names(dat))] )
tissue_FACS <- sub("_.*", "", names(dat_FACS))
method_FACS <- sub(".*.rds.", "", names(dat_FACS))

FACS_plot <- vector()
for(a in unique(method_FACS)) {
	FACS_plot <- rbind(FACS_plot, dat_FACS[method_FACS==a])
}
rownames(FACS_plot) <- remap_names(unique(method_FACS))
colnames(FACS_plot) <- unique(tissue_FACS)
FACS_plot <- FACS_plot[match(names(Method_Legend), rownames(FACS_plot)),]
if (grep("logcounts", rownames(FACS_plot))) {
	FACS_plot <- FACS_plot[-grep("logcounts", rownames(FACS_plot)),]
}
require("RColorBrewer")
tissue_colours <- c(brewer.pal(12, "Set3"), "#e5d8bd")
names(tissue_colours) <- colnames(FACS_plot)

par(mar=c(4.5,4,2,1))
barplot(t(FACS_plot), beside=T, ylim=c(0,1), ylab="FPR", col=tissue_colours, las=2, main="Smart-seq2")
mtext("A", side=2, line=2.5, at=1, cex=1.2, font=2, las=2)


# 10X
dat_10X <- unlist( dat[grep("10X", names(dat))] )
tissue_10X <- sub("_.*", "", names(dat_10X))
method_10X <- sub(".*.rds.", "", names(dat_10X))

tm10X_plot <- vector()
for(a in unique(method_10X)) {
	tm10X_plot <- rbind(tm10X_plot, dat_10X[method_10X==a])
}
rownames(tm10X_plot) <- remap_names(unique(method_10X))
colnames(tm10X_plot) <- unique(tissue_10X)
tm10X_plot <- tm10X_plot[match(names(Method_Legend), rownames(tm10X_plot)),]
if (grep("logcounts", rownames(tm10X_plot))) {
	tm10X_plot <- tm10X_plot[-grep("logcounts", rownames(tm10X_plot)),]
}

par(mar=c(4.5,4,2,1))
barplot(t(tm10X_plot), beside=T, ylim=c(0,1), ylab="FPR", col=tissue_colours[match(colnames(tm10X_plot), names(tissue_colours))], las=2, main="10X")
mtext("B", side=2, line=2.5, at=1, cex=1.2, font=2, las=2)

blank_plot <- function() {
        tmp <-  par("mar")
        par(mar=c(0,0,0,0))
        plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
        return(tmp);
}
blank_plot()
legend("left", fill=tissue_colours, names(tissue_colours), bty="n");

dev.off()
