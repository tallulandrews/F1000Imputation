source("/nfs/users/nfs_t/ta6/MAGIC/Splatter_Functions.R")
source("~/MAGIC/Colour_Scheme.R")
require("scater")

param_table <- Generate_param_table()

# Make Table
out_table_stats <- vector();
out_table_params <- vector();
for (i in 1:nrow(param_table)) {
	filename <- paste("vsSimParams_accuracy_", i, ".rds", sep="");
	all <- readRDS(filename)
	for (meth in names(all)) {
		x <- all[[meth]]
		names(x$stats) <- c("TP", "FP", "TN", "FN");
		out_table_stats <- rbind(out_table_stats, x$stats)
		stuff <- c(meth, param_table[i, "ngenes"], 
				param_table[i,"ncells"], 
				param_table[i, "dropouts"], 
				param_table[i, "ngroups"], 
				param_table[i, "seed"], 
				param_table[i, "propDE"], 
				param_table[i, "method"])
		names(stuff) <- c("imputation", "ngenes", "ncells", "dropouts", "ngroups", "seed", "propDE", "method")
		out_table_params <- rbind(out_table_params, stuff)
	}
}
saveRDS(list(stats=out_table_stats, params=out_table_params), file="vsSimParam_ResultsTable.rds")

tmp <- readRDS("vsSimParam_ResultsTable.rds")
out_table_params <-  tmp$params
out_table_stats <- tmp$stats

all <- cbind(out_table_params, out_table_stats)
rownames(all) <- NULL
all <- as.data.frame(all)

require("matrixStats")
FPR <- out_table_stats[,"FP"]/(out_table_stats[,"FP"]+ out_table_stats[,"TN"])
specificity <- out_table_stats[,"TN"]/(out_table_stats[,"FP"]+ out_table_stats[,"TN"])
precision <- out_table_stats[,"TP"]/(out_table_stats[,"FP"]+out_table_stats[,"TP"])
recall<- out_table_stats[,"TP"]/(out_table_stats[,"FN"]+out_table_stats[,"TP"])
all$FPR <- as.numeric(FPR)
all$specificity <- as.numeric(specificity)
all$precision <- as.numeric(precision)
all$recall <- as.numeric(recall)
all$imputation <- as.character(all$imputation)
all$imputation <- Method_Colours[match(all$imputation, Method_Colours[,1]),2]
all$imputation <- names(Method_Legend)[match(all$imputation, Method_Legend)]
all$imputation <- factor(all$imputation, levels=names(Method_Legend))

# Change NA to "0" for added dropouts.
tmp <- unique(all$dropouts[!is.na(all$dropouts)]); 
all$dropouts <- as.character(all$dropouts)
all$dropouts[is.na(all$dropouts)] <- "0"
all$dropouts <- factor(all$dropouts, levels=c("0", as.character(sort(tmp))))

require("RColorBrewer")
my_grouped_boxplot <- function(x, y, z, cols=brewer.pal(length(unique(z)), "Greys"), place_legend=NULL, legend_title=NULL, ...) {
	nx <- length(levels(x))
	ny <- length(levels(y))
	nz <- length(levels(z))
	xes <- 1:(nx*(nz+1))
	xes <- xes[-seq(from=nz+1, to=length(xes), by=nz+1)]
	thing <- boxplot(y~z+x, at=xes, col=cols, xaxt="n", names=rep("", length(xes)), ...)
	x_names <- rep(levels(x), each=nz)
	x_axis_info <- aggregate(xes, by=list(x_names), mean)
	axis(1, at=x_axis_info[,2], labels=x_axis_info[,1], las=2)
	if (!is.null(place_legend)) {
		legend("top", bty="n", fill=cols, levels(z), title=legend_title, horiz=TRUE, inset=-0.2, xpd=TRUE)
	}
	invisible(thing)
}

png("Splatter_vs_SimParam_review.png", width=8, height=8, units="in", res=300)
par(mfrow=c(2,2)) 
par(mar=c(5,4,3,1))
boxplot(recall~imputation, all, las=2, ylab="Sensitivity", col=Method_Legend)
mtext("A", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
boxplot(specificity~imputation, all, las=2, ylab="Specificity", col=Method_Legend)
mtext("B", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)


cols=brewer.pal(length(unique(all$dropouts)), "Greys")
my_grouped_boxplot(all$imputation, all$recall, all$dropouts, ylab="Sensitivity", place_legend=NULL)
mtext("C", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
legend("topleft", bty="n", col=c(NA,rep("black", length(cols))), pt.bg=c(NA, cols), 
	c("Dropout:", levels(all$dropouts)),title="", horiz=TRUE, inset=c(0,-0.2), 
	xpd=TRUE, pch=22, pt.cex=1.5, x.intersp=1, text.width=c(1,6,5,3.5,2.85,2.5))

cols=brewer.pal(length(unique(all$propDE)), "Greys")
my_grouped_boxplot(all$imputation, all$specificity, all$propDE, ylab="Specificity", place_legend=NULL, legend_title="DE")
legend("topleft", bty="n", col=c(NA,rep("black", length(cols))), pt.bg=c(NA, cols), 
	c("DE:", levels(all$propDE)),title="", horiz=TRUE, inset=c(0,-0.2), 
	xpd=TRUE, pch=22, pt.cex=1.5, x.intersp=1, text.width=c(1,2,2,2))
mtext("D", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
dev.off()

png("Splatter_Suppl1_review.png", width=8, heigh=4, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
my_grouped_boxplot(all$imputation, all$specificity, all$dropouts, ylab="Specificity", place_legend="topleft", legend_title="Dropout")
mtext("A", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
my_grouped_boxplot(all$imputation, all$recall, all$propDE, ylab="Sensitivity", place_legend="bottomleft", legend_title="DE")
mtext("B", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
dev.off()

### Review Figures ###



### ROC ###
source("~/MAGIC/CTP_Functions.R")
	
get.this.stats <- function(stats, sig_threshold=0.05, mag_threshold=0) {
	TP <- sum(
		stats$ROC.info$dir.correct & 
		stats$ROC.info$truth & 
		stats$ROC.info$magnitude >= mag_threshold & 
		stats$ROC.info$q.values <= sig_threshold)
	FP <- sum(
		((!stats$ROC.info$truth) | 
		(!stats$ROC.info$dir.correct) ) &
		(stats$ROC.info$q.values <= sig_threshold &
		stats$ROC.info$magnitude >= mag_threshold)
		)
	TN <- sum(
		(stats$ROC.info$q.values > sig_threshold |
		stats$ROC.info$magnitude < mag_threshold) &
		!stats$ROC.info$truth)
	FN <- sum(
		stats$ROC.info$truth &
		(stats$ROC.info$q.values > sig_threshold |
                stats$ROC.info$magnitude < mag_threshold)
		)
	TN <- sum(
		(!stats$ROC.info$truth) &
		(stats$ROC.info$q.values > sig_threshold |
                stats$ROC.info$magnitude < mag_threshold)
                )
	return(c(TP, FP, TN, FN))
}
		

get.ROC <- function(stats, thresholds) {
	roc <- t(sapply(thresholds, function(th) get.this.stats(stats, th)))
	colnames(roc) <- c("TP", "FP", "TN", "FN")
	roc <- data.frame(thresholds, roc, tpr=roc[,1]/(roc[,1]+roc[,4]), fpr=roc[,2]/(roc[,2]+roc[,3]), spec=roc[,3]/(roc[,3]+roc[,2]), prec=roc[,1]/(roc[,1]+roc[,2]))
	return(roc)
}

plot.ROC <- function(files, method_col_tab, n=500){
	thresholds <- seq(-100,log10(0.05),length.out=n/2)
	thresholds <- c(10^thresholds, seq(0.05, 1, length.out=n/2))
	fdr5 <- min(which(signif(thresholds, digits=1)==0.05))
	out_tpr <- list()
	out_fpr <- list()	
	labs <- vector()
	fileid <- vector()
	for (f in files) {
		obj <- readRDS(f)
	for (meth in names(obj)) {
		res <- obj[[meth]]
		roc <- get.ROC(res, thresholds);
		out_tpr <- cbind(out_tpr, as.numeric(unlist(roc$spec)))
		out_fpr <- cbind(out_fpr, as.numeric(unlist(roc$tpr)))
		labs <- c(labs, meth);
		fileid <- c(fileid, f)
	}
	}
	out_fpr <- matrix(as.numeric(out_fpr), byrow=F, ncol=ncol(out_fpr))
	out_tpr <- matrix(as.numeric(out_tpr), byrow=F, ncol=ncol(out_tpr))
	out_tpr[is.na(out_tpr)] <- 0
	out_fpr[is.na(out_fpr)] <- 0
	xes <- my_row_mean_aggregate(out_fpr, labs)
	yes <- my_row_mean_aggregate(out_tpr, labs)
	line_cols <- method_col_tab[match(colnames(xes), method_col_tab[,1]),]
	plot(1,1, col="white", xlim=c(0,1), ylim=c(0,max(yes)), xlab="Sensitivity", ylab="Specificity", bty="l")
	for (i in 1:ncol(xes)) {
if (!is.na(line_cols[i,2])) {
		lines(xes[,i], yes[,i], col=line_cols[i,2], lwd=2)
		points(xes[fdr5,i], yes[fdr5,i], pch=16, col=line_cols[i,2], cex=1.5)
}
	}
	return(list(line_dat=data.frame(thresh=thresholds, tpr=yes, fpr=xes), line_col=line_cols))
}

param_table <- Generate_param_table()
#png("Figure2_vsSimParams_ROCplot_review.png", width=10, height=4, units="in", res=300)
png("Figure2_vsSimParams_ROCplot_review.png", width=6, height=8, units="in", res=300)
#layout(rbind(c(1,2,5,6), c(3,4,5,6)), widths=c(7,7,9,4), heights=c(1,1)) # horizontal
layout(rbind(c(1,1,2,2), c(3,3,4,4), c(5,5,5,6)), widths=c(7,3,7,3), heights=c(1,1,2)) # vertical

par(mar=c(5,4,2,1))
# Sensitivity vs Dropouts
cols=brewer.pal(length(unique(all$dropouts)), "Greys")
my_grouped_boxplot(all$imputation, all$recall, all$dropouts, ylab="Sensitivity", place_legend=NULL)
mtext("A", at=1, line=2.5, las=2, cex=1.35, font=2, side=2)
legend("topleft", bty="n", col=c(NA,rep("black", length(cols))), pt.bg=c(NA, cols), 
        c("Drops:", levels(all$dropouts)),title="", horiz=TRUE, inset=c(0,-0.3), 
        xpd=TRUE, pch=22, pt.cex=1.5, x.intersp=1, text.width=c(0,6.0,3.5,3,2.85,2.5))

# Specificity vs Dropouts
cols=brewer.pal(length(unique(all$dropouts)), "Greys")
my_grouped_boxplot(all$imputation, all$specificity, all$dropouts, ylab="Specificity", place_legend=NULL)
mtext("B", at=1, line=2.5, las=2, cex=1.35, font=2, side=2)
legend("topleft", bty="n", col=c(NA,rep("black", length(cols))), pt.bg=c(NA, cols), 
        c("Drops:", levels(all$dropouts)),title="", horiz=TRUE, inset=c(0,-0.3), 
        xpd=TRUE, pch=22, pt.cex=1.5, x.intersp=1, text.width=c(0,6.0,3.5,3,2.85,2.5))

# Sensivity vs DE
cols=brewer.pal(length(unique(all$propDE)), "Greys")
my_grouped_boxplot(all$imputation, all$recall, all$propDE, ylab="Sensitivity", place_legend=NULL, legend_title="DE")
mtext("C", at=1, line=2.5, las=2, cex=1.35, font=2, side=2)
legend("topleft", bty="n", col=c(NA,rep("black", length(cols))), pt.bg=c(NA, cols),
        c("DE:", levels(all$propDE)),title="", horiz=TRUE, inset=c(0,-0.3),
        xpd=TRUE, pch=22, pt.cex=1.5, x.intersp=1, text.width=c(2.5,2.5,2.5,2.5))

# Specificity vs DE
cols=brewer.pal(length(unique(all$propDE)), "Greys")
my_grouped_boxplot(all$imputation, all$specificity, all$propDE, ylab="Specificity", place_legend=NULL, legend_title="DE")
mtext("D", at=1, line=2.5, las=2, cex=1.35, font=2, side=2)
legend("topleft", bty="n", col=c(NA,rep("black", length(cols))), pt.bg=c(NA, cols),
        c("DE:", levels(all$propDE)),title="", horiz=TRUE, inset=c(0,-0.3),
        xpd=TRUE, pch=22, pt.cex=1.5, x.intersp=1, text.width=c(2.5,2.5,2.5,2))

# ROC
par(mar=c(4,4,1,1))
files_orig <- paste("/lustre/scratch117/cellgen/team218/TA/Simulations_Temporary_Files/Imputation/vsSimParams_accuracy_", 1:60, ".rds", sep="")
par(mar=c(4,4,1,1))
res1 <- plot.ROC(files_orig, Method_Colours, n=500)
mtext("E", at=1, line=2.5, las=2, cex=1.35, font=2, side=2)
#abline(a=0, b=1, lty=2)
par(mar=c(0,0,0,0))
plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
legend("left", names(Method_Legend), col=Method_Legend, lwd=2, bty="n")
dev.off()

#### Batches ####
param_table <- Generate_param_table()

# Make Table
out_table_stats <- vector();
out_table_params <- vector();
for (i in 1:nrow(param_table)) {
        filename <- paste("Batches/vsSimParams_wBatch_accuracy_", i, ".rds", sep="");
        all <- readRDS(filename)
        for (meth in names(all)) {
                x <- all[[meth]]
                names(x$stats) <- c("TP", "FP", "TN", "FN");
                out_table_stats <- rbind(out_table_stats, x$stats)
                stuff <- c(meth, param_table[i, "ngenes"],
                                param_table[i,"ncells"],
                                param_table[i, "dropouts"],
                                param_table[i, "ngroups"],
                                param_table[i, "seed"],
                                param_table[i, "propDE"],
                                param_table[i, "method"])
                names(stuff) <- c("imputation", "ngenes", "ncells", "dropouts", "ngroups", "seed", "propDE", "method")
                out_table_params <- rbind(out_table_params, stuff)
        }
}
saveRDS(list(stats=out_table_stats, params=out_table_params), file="vsSimParam_wBatch_ResultsTable.rds")

# ROC
png("Figure2_vsSimParams_ROCplot_batch.png", width=6, height=4, units="in", res=300)
layout(rbind(c(1,2), c(1,2)), widths=c(6,2), heights=c(1,1))
par(mar=c(4,4,1,1))
files <- paste("/lustre/scratch117/cellgen/team218/TA/Simulations_Temporary_Files/Imputation/Batches/vsSimParams_wBatch_accuracy_", 1:60, ".rds", sep="")
par(mar=c(4,4,1,1))
res2 <- plot.ROC(files, Method_Colours, n=500)
#abline(a=0, b=1, lty=2)
par(mar=c(0,0,0,0))
plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
legend("left", names(Method_Legend), col=Method_Legend, lwd=2, bty="n")
dev.off()

#### vs Magnitude Threshold ####
source("~/MAGIC/Colour_Scheme.R")
param_table <- Generate_param_table()

# Make Tables
sig_thresh <- 0.05
mag_perc <- c(100, 90, 80, 70, 60, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 1)
OUT_sens <- list()
OUT_spec <- list()

for (i in 1:nrow(param_table)) {
	dat <- readRDS(paste("vsSimParams_accuracy_", i, ".rds", sep=""))
	for (meth in names(dat)) {
		res <- dat[[meth]][["ROC.info"]]
		if (i == 1) {
			OUT_sens[[meth]] <- matrix(0, nrow=length(mag_perc), ncol=nrow(param_table))
			OUT_spec[[meth]] <- matrix(0, nrow=length(mag_perc), ncol=nrow(param_table))
		}
		for (mag_t in 1:length(mag_perc)) {
			mag_thresh <- quantile(res[,2], prob=(1-mag_perc[mag_t]/100));
			Pos <- res[,1] < sig_thresh & res[,2] > mag_thresh
			TP <- sum(Pos & res[,3] & res[,4])
			FP <- sum(Pos & !(res[,3] & res[,4]))
			TN <- sum(!Pos & !res[,4])
			FN <- sum(!Pos & res[,4])

			sens <- TP/(TP+FN)
			spec <- TN/(TN+FP)
			OUT_sens[[meth]][mag_t, i] <- sens
			OUT_spec[[meth]][mag_t, i] <- spec
		}
	}
}


col_Sens <- "forestgreen"
col_Spec <- "dodgerblue"
methods <- names(OUT_sens)
method_cols <- Method_Colours[match(methods, Method_Colours[,1]),2]
png("Figure2_Splatter_MagThreshold.png", width=7, height=5, units="in", res=300)
layout(rbind(c(1, 2, 3), c(4,5,6)))
null_sens <- vector()
null_spec <- vector()
for (m in names(Method_Legend)) {
	require("matrixStats")

	meth <- methods[method_cols == Method_Legend[m]]
	meth <- meth[!is.na(meth)]
	if (grepl("counts", m)) {
		null_sens <- rowMeans(OUT_sens[[meth]])
		null_spec <- rowMeans(OUT_spec[[meth]])
		next;
	}
	
	par(mar=c(3.5, 3.5, 1.5, 0.5))
	plot(mag_perc, rowMeans(OUT_sens[[meth]]), col=col_Sens, ylim=c(0,1), xlab="", ylab="", type="b", main="")
	title(main=m, line=0.5)
	title(xlab="Magnitude Threshold (%)", line=2)
	title(ylab="Score", line=2)
	lines(mag_perc, null_sens, col="grey75", lwd=5)
	lines(mag_perc, null_spec, col="grey75", lwd=5)
	
	lines(mag_perc, rowMeans(OUT_sens[[meth]]), col=col_Sens, type="l")
	sdev <- sqrt(rowVars(OUT_sens[[meth]]))/sqrt(ncol(OUT_sens[[meth]]));
	lines(mag_perc, rowMeans(OUT_sens[[meth]])+1.96*sdev, lty=2, col=col_Sens)
	lines(mag_perc, rowMeans(OUT_sens[[meth]])-1.96*sdev, lty=2, col=col_Sens)

	lines(mag_perc, rowMeans(OUT_spec[[meth]]), type="l", col=col_Spec)
	lines(mag_perc, rowMeans(OUT_spec[[meth]]), type="b", col=col_Spec)
	sdev <- sqrt(rowVars(OUT_spec[[meth]]))/sqrt(ncol(OUT_spec[[meth]]));
	lines(mag_perc, rowMeans(OUT_spec[[meth]])+1.96*sdev, lty=2, col=col_Spec)
	lines(mag_perc, rowMeans(OUT_spec[[meth]])-1.96*sdev, lty=2, col=col_Spec)
	if (m %in% c("knn","DrImpute")) {
	legend("bottomright", c("Sensitivity", "Specificity"), col=c(col_Sens, col_Spec), lty=1, pch=1, bty="n")
	}
}
dev.off()


