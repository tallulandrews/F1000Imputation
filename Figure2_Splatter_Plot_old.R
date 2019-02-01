source("/nfs/users/nfs_t/ta6/MAGIC/Splatter_Functions.R")
source("~/MAGIC/Colour_Scheme.R")
require("scater")

param_table <- Generate_param_table()

# Make Table
out_table_stats <- vector();
out_table_params <- vector();
for (i in 1:nrow(param_table)) {
	filename <- paste("vsSimParams_accuracy_paramset_", i, ".rds", sep="");
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
saveRDS(list(stats=out_table_stats, params=out_table_params), file="vsSimParam_FinalTable.rds")

all <- cbind(out_table_params, out_table_stats)
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
all$imputation <- factor(all$imputation, levels=imputation_methods)

# Change NA to "None" for added dropouts.
tmp <- unique(all$dropouts[!is.an(all$dropouts)]); 
all$dropouts <- as.character(all$dropouts)
all$dropouts[is.na(all$dropouts)] <- "None"
all$dropouts <- factor(all$dropouts, levels=c("None", as.character(sort(tmp))))

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
		legend(place_legend, bty="n", fill=cols, levels(z), title=legend_title)
	}
}

png("Splatter_vs_SimParam.png", width=8, height=8, units="in", res=300)
par(mfrow=c(2,2)) 
par(mar=c(5,4,1,1))
boxplot(recall~imputation, all, las=2, ylab="Sensitivity", col="grey50")
mtext("A", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
boxplot(specificity~imputation, all, las=2, ylab="Specificity", col="grey50")
mtext("B", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
my_grouped_boxplot(all$imputation, all$recall, all$dropouts, ylab="Sensitivity", place_legend="topleft", legend_title="Dropout")
mtext("C", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
my_grouped_boxplot(all$imputation, all$specificity, all$propDE, ylab="Specificity", place_legend="bottomleft", legend_title="DE")
mtext("D", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
dev.off()

png("Splatter_Suppl1.png", width=8, heigh=4, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
my_grouped_boxplot(all$imputation, all$specificity, all$dropouts, ylab="Specificity", place_legend="topleft", legend_title="Dropout")
mtext("A", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
my_grouped_boxplot(all$imputation, all$recall, all$propDE, ylab="Sensitivity", place_legend="bottomleft", legend_title="DE")
mtext("B", at=1, line=2.5, las=2, cex=1.65, font=2, side=2)
dev.off()


### ROC ###
get.tpr <- function(stats, sig_threshold=0.05, mag_threshold=0) {
	Pos <- stats$ROC.info$q.values < sig_threshold & stats$ROC.info$mag > mag_threshold
	True <- stats$ROC.info$dir.correct & stats$ROC.info$truth
	sum(Pos & True)/sum(True)
}
	
get.fpr <- function(stats, sig_threshold=0.05, mag_threshold=0) {
	Pos <- stats$ROC.info$q.values < sig_threshold & stats$ROC.info$mag > mag_threshold
	False <- !stats$ROC.info$dir.correct || !stats$ROC.info$truth
	sum(Pos & False)/sum(False)
}
	

get.ROC <- function(stats, thresholds) {
	roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
	roc$tpr <- sapply(roc$threshold, function(th) get.tpr(stats, th))
	roc$fpr <- sapply(roc$threshold, function(th) get.fpr(stats, th))
	return(roc)
}

plot.ROC <- function(files, method_col_tab, n=100){
	thresholds <- seq(0,1,length.out=n)
	out_tpr <- list()
	out_fpr <- list()	
	for (f in files) {
		obj <- readRDS(f)
	for (meth in names(obj)) {
		res <- obj[[method]]
		roc <- get.ROC(res, thresholds);
		out_tpr <- cbind(out_tpr, roc$tpr)
		out_fpr <- cbind(out_fpr, roc$fpr)
	}
	}
	xes <- rowMeans(out_fpr)
	yes <- rowMeans(out_tpr)
	return(data.frame(thresh=thresholds, tpr=yes, fpr=xes))
}
