require("scater")
source("~/MAGIC/Colour_Scheme.R")

# args: rds of simulated & imputed results, name of output
do_wilcox_DE <- function(obj, mat_name, norm=FALSE) {
	require("Hmisc")
	source("/nfs/users/nfs_t/ta6/MAGIC/CTP_Functions.R")
	require("scater")
	require("Matrix")
	groups <- unique(obj$Group)
	mat <- obj@assays[[mat_name]]
	if (norm) {
		sf <- colSums(mat);
		mat <- t( t(mat)/sf*median(sf) )
	}
	Out <- vector();
	for (i in 1:(length(groups)-1)) {
		for (j in (i+1):length(groups)) {
			p.val <- apply(mat, 1, function(a) {wilcox.test(a[obj$Group==groups[i]], a[obj$Group==groups[j]])$p.value} )
			m1 <- Matrix::rowMeans(mat[,obj$Group==groups[i]])
			m2 <- Matrix::rowMeans(mat[,obj$Group==groups[j]])
			de_genes <- rownames(mat)
			de_g1 <- rep(groups[i], nrow(mat))
			de_g1[m2 > m1] <- groups[j]
			de_g2 <- rep(groups[j], nrow(mat))
			de_g2[m2 > m1] <- groups[i]
			tab <- cbind(as.character(de_genes), as.character(de_g1), as.character(de_g2), p.adjust(p.val, method="fdr"))
			Out <- rbind(Out, tab);
		}
	}
	return(Out);
}

get.this.stats <- function(DE_out, truth, sig_threshold=0.05) {
	sig <- as.numeric(DE_out[,4]) < sig_threshold;
	dir.correct <- paste(DE_out[,1], DE_out[,2], DE_out[,3]) %in% paste(truth$DE[,1], truth$DE[,2], truth$DE[,3])
	is.de <- DE_out[,1] %in% truth$DE[,1]
	is.not.de <- DE_out[,1] %in% truth$nonDE

	TP <- sum(sig & dir.correct & is.de)
	FP <- sum((sig & is.not.de) | (sig & !dir.correct & is.de))
	TN <- sum(!sig & is.not.de)
	FN <- sum(!sig & is.de)
	return(c(TP, FP, TN, FN))
}
		

get.ROC <- function(de_out, truth) {
	thresholds <- c(0, 10^seq(from=-10, to=-1, length=100), seq(from=0.1, to=1, length=100))
	roc <- t(sapply(thresholds, function(th) get.this.stats(de_out, truth, th)))
	colnames(roc) <- c("TP", "FP", "TN", "FN")
	roc <- data.frame(thresholds, roc, tpr=roc[,1]/(roc[,1]+roc[,4]), fpr=roc[,2]/(roc[,2]+roc[,3]), spec=roc[,3]/(roc[,3]+roc[,2]), prec=roc[,1]/(roc[,1]+roc[,2]))
	return(roc)
}

get_ROC_curve <- function(obj, mat_name, truth, norm=FALSE) { 
	de <- do_wilcox_DE(obj, mat_name, norm)
	roc <- get.ROC(de, truth)
	return(roc)
}

# get ROC curve for each method for each bulk dataset. & plot with appropriate colours
# Analyze DE in simulated datasets:
accuracy <- list()
obj <- readRDS(infile)
methods <- names(assays(obj))

methods <- methods[(grep("counts", methods)[1]):length(methods)]
for (meth in methods) {
	if (meth == "counts") {
		out <- get_ROC_curve(obj, mat_name=meth, norm=TRUE);
	} else {
		out <- get_ROC_curve(obj, mat_name=meth);
	}
	accuracy[[meth]] <- out;
}
method_colours <-

