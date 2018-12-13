Generate_param_table <- function() {
	set.seed(20194)
	nrep=10;
	my_seeds <- round(runif(nrep*2)*10000)
	my_seeds <- unique(my_seeds)
	my_seeds <- my_seeds[1:nrep]
	sim_params <- list(dropouts=c(NA, 1, 3, 5, 7),
                   propDE=c(0.01, 0.1, 0.3),
                   method=c("groups"),
                   ngroups=c(2,5,10),
                   ngenes=c(1000, 2000, 5000),
                   seeds=my_seeds[1:4])
	param_table <- data.frame(
        dropouts=rep(sim_params$dropouts, times=length(sim_params$ngroups)*length(sim_params$method)*length(sim_params$seeds)),
        ngroups=rep(sim_params$ngroups, times=length(sim_params$dropouts)*length(sim_params$method)*length(sim_params$seeds)),
        method=rep(sim_params$method, times=length(sim_params$dropouts)*length(sim_params$ngroups)*length(sim_params$seeds)),
        seed=rep(sim_params$seeds, each=length(sim_params$dropouts)*length(sim_params$method)*length(sim_params$ngroups))
        )
	set.seed(3819)
	param_table$propDE <- sample(sim_params$propDE, nrow(param_table), replace=TRUE)
	param_table$ngenes <- sample(sim_params$ngenes, nrow(param_table), replace=TRUE)
	param_table$ncells <- rep(1000, times=nrow(param_table));
	return(param_table)
}

# Fixed an error in this function 6 Dec 2018
Sim_w_splatter <- function(ngenes=1000, ncells=1000, ngroups=2, seed=101, dropouts=3, propDE=0.1, method="groups") {
	require("splatter")
	require("scater")
	if (is.null(dropouts) | is.na(dropouts)) {
		sim <- splatSimulate( nGenes = ngenes, 
				      batchCells=ncells, 
				      group.prob=rep(1/ngroups, times=ngroups), 
				      method=method, seed=seed, 
				      de.prob=propDE/ngroups)
	} else {
		sim <- splatSimulate(nGenes = ngenes, 
				     batchCells=ncells, 
				     group.prob=rep(1/ngroups, times=ngroups), 
				     method=method, seed=seed, 
				     dropout.present=TRUE, dropout.mid=dropouts, 
				     de.prob=propDE/ngroups)
	}
	sim <- normalise(sim)
	return(sim);
}

do_imputation_splatter <- function(splatter_sims, method=c("MAGIC", "scImpute", "DrImpute", "SAVER", "knn"), n.cores=1) {
	require("scater")
	source("/nfs/users/nfs_t/ta6/MAGIC/knn_smooth.R")

	if (method == "scImpute") {
		default_param=0.5;
            	saveRDS(assays(splatter_sims)[["counts"]], file="tmp.rds");
                scImpute::scimpute("./tmp.rds", infile="rds", outfile="rds",
                        type="count", drop_thre=default_param, out_dir="./",
                        Kcluster=length(unique(splatter_sims$Group)), ncores=n.cores)
                out <- readRDS("./scimpute_count.rds")
                return(out);
	} else if (method == "DrImpute") {
		default_param=0;
                out <- DrImpute::DrImpute(assays(splatter_sims)[["logcounts"]], 
			ks=length(unique(splatter_sims$Group)),
                        zerop=default_param)
                return(out);
	} else if (method == "SAVER") {
		# all genes default
		out <- saver(assays(splatter_sims)[["counts"]], do.fast=TRUE, size.factor=1)
                return(out$estimate)
	} else if (method == "MAGIC") {
		default_param=0;
                out <- Rmagic::run_magic(t(assays(splatter_sims)[["counts"]]),
                                t_diffusion=default_param, lib_size_norm=T,
                                log_transform=F, pseudo_count=1, npca=100,
                                k=12, ka=4, epsilon=1, rescale_percent=0)
                return(t(out));
	} else if (method == "knn") {
		default_param=ncol(sim)/20;
                out <- knn_smoothing(assays(splatter_sims)[["counts"]], k=default_param, d=10, seed=42)
                return(out)
        }
}

check_multiDE_accuracy_splatter <- function(sim, mat_name="logcounts", mag_centile=NULL) {
	require("Hmisc")
	source("/nfs/users/nfs_t/ta6/MAGIC/CTP_Functions.R")
	require("scater")
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
	# Non-parametric DE using Kruskal-Wallis test
	p_values <- apply(assays(sim)[[mat_name]], 1, function(x) { kruskal.test(x, factor(sim$Group))$p.value})
	p_values[is.na(p_values)] <- 1;
	# FDR multiple testing correction
	q_values <- p.adjust(p_values, method="fdr")
	# Observed FC
	lvls <- my_row_mean_aggregate(assays(sim)[[mat_name]], factor(sim$Group)) # Obs mean expression by group
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

