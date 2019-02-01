require("scater")

# args: rds of simulated & imputed results, name of output
args <- commandArgs(trailingOnly=TRUE)
infile=args[1]
outfile=args[2]

check_multiDE_accuracy_splatter <- function(sim, mat_name="logcounts", mag_centile=NULL, norm=FALSE) {
	require("Hmisc")
	source("/nfs/users/nfs_t/ta6/MAGIC/CTP_Functions.R")
	require("scater")
	require("Matrix")
	# Ground Truth
	de_cols <- grep("DEFac", colnames(rowData(sim)))
	de_tab <- rowData(sim)[,de_cols]
	high <- apply(de_tab, 1, max)
	low <- apply(de_tab, 1, min)
	is.de <- (high-low) != 0
	colnames(de_tab) <- sub("DEFac", "", colnames(de_tab))
	de_rows <- which(is.de)
	# make paths more divergent by only considering the tips
	if (grepl("Path",sim$Group[1])){
		sim <- sim[,sim$Step > 50] 
	}
	this_mat <- assays(sim)[[mat_name]]
	if (norm) {
		sf <- Matrix::colSums(this_mat)
		this_mat <- t( t(this_mat)/sf*median(sf) )
	}

	# Non-parametric DE using Kruskal-Wallis test
	p_values <- apply(this_mat, 1, function(x) { kruskal.test(x, factor(sim$Group))$p.value})
	p_values[is.na(p_values)] <- 1;
	# FDR multiple testing correction
	q_values <- p.adjust(p_values, method="fdr")
	# Observed FC
	lvls <- my_row_mean_aggregate(this_mat, factor(sim$Group)) # Obs mean expression by group
	mag <- apply(lvls, 1, max) - apply(lvls, 1, min)

	# Check direction
	# Truth
	de_up <- high > 1
	top <- sapply(which(de_up), function(i){which(unlist(de_tab[i,]) == high[i])})
	de_dn <- low < 1
	bottom <- sapply(which(de_dn), function(i){which(unlist(de_tab[i,]) == low[i])})
	# Observed
	high_obs <- apply(lvls, 1, max)
	low_obs <- apply(lvls, 1, min)
	
	top_obs <- sapply(which(de_up), function(i){
				out <- which(unlist(lvls[i,]) == high_obs[i])
				if (length(out) > 1) { # If ties for highest obs group
					if (top[i] %in% out) { return(top[i])} 
					else {return(out[1])}
				} else { return(out)}
				})
	bottom_obs <- sapply(which(de_dn), function(i){
				out <- which(unlist(lvls[i,]) == low_obs[i])
				if (length(out) > 1) { # If ties for lowest obs group
					if (bottom[i] %in% out) {return(bottom[i])}
					else {return(out[1])}
				} else {return(out)}
				})
	# Deal with if both up & down
	dir.correct <- rep(TRUE, times=length(is.de))
	dir.correct[de_up] <- dir.correct[de_up] & top == top_obs # Is the group with highest obs expression also the group with true highest expression?
	dir.correct[de_dn] <- dir.correct[de_dn] & bottom == bottom_obs # Is the group with lowest obs expression also the group with true lowest expression?


	thresh <- 0.05 #significance
	mag_thresh <- 0 #magnitude
	if (!is.null(mag_centile)) {
		mag_thresh <- quantile(mag, probs=1-mag_centile)
	}

        TP <- sum(q_values < thresh & is.de & dir.correct & mag > mag_thresh, na.rm=T)
        FP <- sum(q_values < thresh & mag > mag_thresh, na.rm=T) - TP
        FN <- sum( (q_values >= thresh | mag <= mag_thresh) & is.de, na.rm=T)
        TN <- sum( (q_values >= thresh | mag <= mag_thresh), na.rm=T) - FN

	# Get examples of the false positives
	FPs <- which( !is.de & q_values < thresh, arr.ind=T )
	if ( sum(!is.de & q_values < thresh) > 100) {
		FPs_p <- q_values[FPs]
		FPs_r <- mag[FPs]
		most_FP <- which(abs(FPs_r) >= quantile(abs(FPs_r), prob=1-100/length(FPs_r)))	

		FPs <- data.frame(gene=FPs, p=FPs_p, r=FPs_r)
		FPs <- FPs[most_FP,]
	}
	wrong_dir <- which(q_values < thresh & is.de & ! dir.correct)
	if (length(wrong_dir) > 0) {
		wrong_dir <- data.frame(wrong_dir, q_values[wrong_dir], mag[wrong_dir])
	}

	# Collect results for making an ROC curve
	ROC_info <- data.frame(q.values=q_values, magnitude=mag, dir.correct=dir.correct, truth=is.de)

	return(list(stats=c(TP, FP, TN, FN), wrongway=wrong_dir, mostfp=FPs, ROC.info=ROC_info))
}

# Analyze DE in simulated datasets:
accuracy <- list()
sim <- readRDS(infile)
methods <- names(assays(sim))

methods <- methods[(grep("counts", methods)[1]):length(methods)]
for (meth in methods) {
	if (meth == "counts") {
		out <- check_multiDE_accuracy_splatter(sim, mat_name=meth, norm=TRUE);
	} else {
		out <- check_multiDE_accuracy_splatter(sim, mat_name=meth);
	}
	accuracy[[meth]] <- out;
}
saveRDS(accuracy, file=outfile)

# Analyze DE in simulated datasets, different magitude thresholds:
accuracy <- list()
sim <- readRDS(infile)
methods <- names(assays(sim))
mag_thresh <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

methods <- methods[(grep("counts", methods)[1]):length(methods)]
for (meth in methods) {
	if (meth == "counts") {
		out <- check_multiDE_accuracy_splatter(sim, mat_name=meth, norm=TRUE, mag_thresh=m_thresh);
	} else {
		out <- check_multiDE_accuracy_splatter(sim, mat_name=meth, mag_thresh=m_thresh);
	}
	accuracy[[meth]] <- out;
}
saveRDS(accuracy, file=outfile)
