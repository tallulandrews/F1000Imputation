Sim <- function(n_genes=500, n_cells=500, l10fc=1, disp=0.2, min_mean=-3, max_mean=6, propDE=0.5, type=c("clust", "time"), minor_type_freq=0.5, d_rate=0) {
	mus <- runif(n_genes, min=min_mean+0.01, max=max_mean-0.01)
	up <- 1:(n_genes*propDE/2)
	down <- (max(up)+1):((n_genes*propDE))

	mus2 <- mus
	mus2[up] <- mus[up]+l10fc
	mus2[down] <- mus[down]-l10fc
	mus2[mus2 > max_mean] <- max_mean
	mus2[mus2 < min_mean] <- min_mean
	time <- runif(n_cells);
	
	n_cells1 <- n_cells * (1-minor_type_freq)
	n_cells2 <- n_cells * minor_type_freq

	mu_tab <- 10^mus %*% t(rep(1, times=n_cells))

	if (type[1]=="clust") {
		mu_tab2 <- 10^mus2 %*% t(rep(1, times=n_cells2))
		truth <- rep(0, time=n_cells);
		mu_tab[,1:n_cells2] <- mu_tab2;
		truth[1:n_cells2] <- 1
	} else if (type[1]=="time") {
		mu_tab2 <- 10^mus2 %*% t(rep(1, times=n_cells))
		mu_tab <- t( t(mu_tab2)*(time) + t(mu_tab)*(1-time) )
		truth=time
	}

	MAT <- sapply(1:n_genes, function(i) {
                        sapply(mu_tab[i,], function(m) {
                                if(runif(1) > d_rate) {rnbinom(1, mu=m, size=1/disp)}
				else{0};
                                })
                        })
	MAT <- t(MAT)
	is.up <- 1:n_genes %in% up
	is.down <- 1:n_genes %in% down

	keep_g <- rowSums(MAT) > 0
	MAT <- MAT[keep_g,]
	is.up <- is.up[keep_g]
	is.down <- is.down[keep_g]

	keep_c <- colSums(MAT) > 0
	MAT <- MAT[,keep_c]
	truth <- truth[keep_c]

	rownames(MAT) <- paste("g", 1:nrow(MAT), sep="");
	colnames(MAT) <- paste("c", 1:ncol(MAT), sep="");
	return(list(mat=MAT, g_up=is.up, g_down=is.down, cell_truth=truth))
}

gene_cor_heatmap <- function(sim, mat="mat") {
	cor_mat <- cor(t(sim[[mat]]), method="spearman")
	require("gplots")
	require("RColorBrewer")
	heat_cols <- rev(brewer.pal(8, "RdBu"))
	gene_info <- data.frame(is.up=sim$g_up, is.down=sim$g_down, expr=rowMeans(sim$mat))
	reorder <- order(gene_info[,1], gene_info[,2], gene_info[,3], decreasing=T)
	gene_info <- gene_info[reorder,]
	sideCols <- rep("grey85", nrow(sim$mat))
	sideCols[gene_info[,1]] <- "red"
	sideCols[gene_info[,2]] <- "blue"
	cor_mat <- cor_mat[reorder, reorder]
	heatmap.2(cor_mat, col=heat_cols, trace="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=sideCols, RowSideColors=sideCols)
	invisible(cor_mat)
}
	
cell_cor_heatmap <- function(sim) {
	cor_mat <- cor(sim$mat, method="spearman")
        require("gplots")
        require("RColorBrewer")
	reorder <- order(sim$cell_truth)
	cor_mat <- cor_mat[reorder,reorder]
	truth <- sim$cell_truth[reorder]

        heat_cols <- brewer.pal(8, "Greys")
	cell_cols <- colorRampPalette(c("white", "black"))(10)
	breaks <- seq(from=0, to=1, length=11)
	sideCols <- cell_cols[cut(truth, breaks=breaks)]
	heatmap.2(cor_mat, col=heat_cols, trace="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=sideCols, RowSideColors=sideCols)
}
	

make_plots_onemethod <- function(gene1, gene2, title=TRUE) {
	plot(sims$counts[gene1,], sims$counts[gene2,], xlab=rownames(sims$counts)[gene1], ylab=rownames(sims$counts)[gene2], pch=1, )
	if (title) {
		title(main="before MAGIC")
	}
	plot(sims$magic[gene1,], sims$magic[gene2,], xlab=rownames(sims$counts)[gene1], ylab=rownames(sims$counts)[gene2], pch=1, )
	if (title) {
                title(main="after MAGIC")
        }
	source("~/R-Scripts/violin_plot.R")
	vioplot(list(sims$magic[gene1,sims$cell_truth==1], sims$magic[gene1,sims$cell_truth==0]), col=c("firebrick", "dodgerblue"), drawRect=FALSE, names=c("type1", "type2"), las=2)
	if (title) {
                title(main=rownames(sims$counts)[gene1])
        }

	vioplot(list(sims$magic[gene2,sims$cell_truth==1], sims$magic[gene2,sims$cell_truth==0]), col=c("firebrick", "dodgerblue"), drawRect=FALSE, names=c("type1", "type2"), las=2)
	if (title) {
                title(main=rownames(sims$counts)[gene2])
        }

}

make_plots_allmethods <- function(gene1, gene2=NULL, title=TRUE, cols=c("firebrick", "dodgerblue")) {
	source("~/R-Scripts/violin_plot.R")
	if (is.null(gene2)) {
		vioplot(list(sims$counts[gene2,sims$cell_truth==1], 
			     sims$counts[gene2,sims$cell_truth==0]), 
			     col=cols, drawRect=FALSE, names=c("type1", "type2"), las=2)
	} else {
		plot(sims$counts[gene1,], sims$counts[gene2,], 
			xlab=rownames(sims$counts)[gene1], 
			ylab=rownames(sims$counts)[gene2], pch=1)
	}
	if (title) {
		title(main="Raw")
	}

	if (is.null(gene2)) {
		vioplot(list(sims$magic[gene2,sims$cell_truth==1], 
			     sims$magic[gene2,sims$cell_truth==0]), 
			     col=cols, drawRect=FALSE, names=c("type1", "type2"), las=2)
	} else {
		plot(sims$magic[gene1,], sims$magic[gene2,], 
			xlab=rownames(sims$counts)[gene1], 
			ylab=rownames(sims$counts)[gene2], pch=1)
	}
	if (title) {
		title(main="MAGIC")
	}

	if (is.null(gene2)) {
		vioplot(list(sims$scimpute[gene2,sims$cell_truth==1], 
			     sims$scimpute[gene2,sims$cell_truth==0]), 
			     col=cols, drawRect=FALSE, names=c("type1", "type2"), las=2)
	} else {
		plot(sims$scimpute[gene1,], sims$scimpute[gene2,], 
			xlab=rownames(sims$counts)[gene1], 
			ylab=rownames(sims$counts)[gene2], pch=1)
	}
	if (title) {
		title(main="scImpute")
	}

	if (is.null(gene2)) {
		vioplot(list(sims$drimpute[gene2,sims$cell_truth==1], 
			     sims$drimpute[gene2,sims$cell_truth==0]), 
			     col=cols, drawRect=FALSE, names=c("type1", "type2"), las=2)
	} else {
		plot(sims$drimpute[gene1,], sims$drimpute[gene2,], 
			xlab=rownames(sims$counts)[gene1], 
			ylab=rownames(sims$counts)[gene2], pch=1)
	}
	if (title) {
		title(main="DrImpute")
	}

	if (is.null(gene2)) {
		vioplot(list(sims$saver[gene2,sims$cell_truth==1], 
			     sims$saver[gene2,sims$cell_truth==0]), 
			     col=cols, drawRect=FALSE, names=c("type1", "type2"), las=2)
	} else {
		plot(sims$saver[gene1,], sims$saver[gene2,], 
			xlab=rownames(sims$counts)[gene1], 
			ylab=rownames(sims$counts)[gene2], pch=1)
	}
	if (title) {
		title(main="SAVER")
	}

	if (is.null(gene2)) {
		vioplot(list(sims$knn_sm[gene2,sims$cell_truth==1], 
			     sims$knn_sm[gene2,sims$cell_truth==0]), 
			     col=cols, drawRect=FALSE, names=c("type1", "type2"), las=2)
	} else {
		plot(sims$knn_sm[gene1,], sims$knn_sm[gene2,], 
			xlab=rownames(sims$counts)[gene1], 
			ylab=rownames(sims$counts)[gene2], pch=1)
	}
	if (title) {
		title(main="kNN Smooth")
	}

}
	
Sim_trio <- function(n_genes=500, n_cells=1000, l10fc=c(1, 0.5), disp=0.2, min_mean=-3, max_mean=4, propDE=c(0.5, 0.25), minor_type_freq=0.33, rare_type_freq=0.33, d_rate=0) {

	n_cells1 <- n_cells * (1-minor_type_freq)
	n_cells2 <- n_cells * minor_type_freq
	n_cells3 <- n_cells * rare_type_freq

	# Type 1
	mus <- runif(n_genes, min=min_mean+0.01, max=max_mean-0.01)
	mu_tab <- 10^mus %*% t(rep(1, times=n_cells))

	# Type 2
	up1 <- 1:(n_genes*propDE[1]/2)
	down1 <- (max(up1)+1):((n_genes*propDE[1]))
	mus2 <- mus
	mus2[up1] <- mus[up1]+l10fc[1]
	mus2[down1] <- mus[down1]-l10fc[1]
	mus2[mus2 > max_mean] <- max_mean
	mus2[mus2 < min_mean] <- min_mean
	time <- runif(n_cells);
	mu_tab2 <- 10^mus2 %*% t(rep(1, times=n_cells2))
	
	# Type 3
	if(is.na(l10fc[2])) {
		# doublets
		up2 <- up1
		down2 <- down1
		mu_tab3 <- cbind((10^mus2 + 10^mus)/2)  %*% t(rep(1, times=n_cells3))
	} else {
		up2 <- runif(n_genes) < propDE[2]/2
		down2 <- runif(n_genes) < propDE[2]/2
		mus3 <- mus2
		mus3[up2] <- mus3[up2] + l10fc[2]
		mus3[down2] <- mus3[down2] - l10fc[2]
		mus3[mus3 > max_mean] <- max_mean
	        mus3[mus3 < min_mean] <- min_mean
		mu_tab3 <- 10^mus3  %*% t(rep(1, times=n_cells3))
	}
	
	truth <- rep(0, time=n_cells);
	truth[1:n_cells2] <- 1
	truth[1:n_cells3+n_cells2] <- 2
	mu_tab[,truth==1] <- mu_tab2;
	mu_tab[,truth==2] <- mu_tab3;

	MAT <- sapply(1:n_genes, function(i) {
                        sapply(mu_tab[i,], function(m) {
                                if(runif(1) > d_rate) {rnbinom(1, mu=m, size=1/disp)}
				else{0};
                                })
                        })
	MAT <- t(MAT)
	is.up <- cbind(1:n_genes %in% up1, 1:n_genes %in% up2)
	is.down <- cbind(1:n_genes %in% down1, 1:n_genes %in% down2)

	keep_g <- rowSums(MAT) > 0
	MAT <- MAT[keep_g,]
	is.up <- is.up[keep_g,]
	is.down <- is.down[keep_g,]

	keep_c <- colSums(MAT) > 0
	MAT <- MAT[,keep_c]
	truth <- truth[keep_c]

	rownames(MAT) <- paste("g", 1:nrow(MAT), sep="");
	colnames(MAT) <- paste("c", 1:ncol(MAT), sep="");
	return(list(mat=MAT, g_up=is.up, g_down=is.down, cell_truth=truth))
}

Sim_w_splatter <- function(ngenes=1000, ncells=1000, ngroups=2, seed=101, dropouts=3, propDE=0.1, method="groups") {
	require("splatter")
	require("scater")
	if (is.null(dropouts) | is.na(dropouts)) {
		sim <- splatSimulate(nGenes = ngenes, batchCells=ncells, group.prob=rep(1/ngroups, times=ngroups), method=method, seed=seed, de.prob=propDE/ngroups)
	} else {
		sim <- splatSimulate(nGenes = ngenes, batchCells=ncells, group.prob=rep(1/ngroups, times=ngroups), method=method, seed=seed, dropout.present=TRUE, dropout.mid=3, de.prob=propDE/ngroups)
	}
	sim <- normalise(sim)
	return(sim);
}

do_imputation_splatter <- function(splatter_sims, method=c("MAGIC", "scImpute", "DrImpute", "SAVER", "knn"), param=NULL, n.cores=1) {

	if (method == "scImpute") {
		default_param=0.5;
		if (!is.null(param)) {
			default_param=param
		}
            	saveRDS(assays(splatter_sims)[["counts"]], file="tmp.rds");
                scImpute::scimpute("./tmp.rds", infile="rds", outfile="rds",
                        type="count", drop_thre=default_param, out_dir="./",
                        Kcluster=length(unique(splatter_sims$Group)), ncores=n.cores)
                out <- readRDS("./scimpute_count.rds")
                return(out);
	} else if (method == "DrImpute") {
		default_param=0;
		if (!is.null(param)) {
			default_param=param
		}
                lognorm <- assays(splatter_sims)[["logcounts"]]
                #out <- DrImpute::DrImpute(lognorm, ks=length(unique(splatter_sims$Group)),
                #        dropout.probability.threshold=default_param)
                out <- DrImpute::DrImpute(lognorm, ks=length(unique(splatter_sims$Group)),
                        zerop=default_param)
                return(out);
	} else if (method == "SAVER") {
		default_param=1;
		if (!is.null(param)) {
			default_param=param
		}
	   	if (default_param == 1) { 
		out <- saver(assays(splatter_sims)[["counts"]], do.fast=TRUE, size.factor=1)
		} else {
		out <- saver(assays(splatter_sims)[["counts"]], do.fast=TRUE, size.factor=1, 
			npred=nrow(splatter_sims)*default_param)
		}
                return(out$estimate)
	} else if (method == "MAGIC") {
		default_param=0;
		if (!is.null(param)) {
			default_param=param
		}
                out <- Rmagic::run_magic(t(assays(splatter_sims)[["counts"]]),
                                t_diffusion=default_param, lib_size_norm=T,
                                log_transform=F, pseudo_count=1, npca=100,
                                k=12, ka=4, epsilon=1, rescale_percent=0)
                return(t(out));
	} else if (method == "knn") {
		default_param=ncol(sim)/20;
		if (!is.null(param)) {
			default_param=param
		}
                out <- knn_smoothing(assays(splatter_sims)[["counts"]], k=default_param, d=10, seed=42)
                return(out)
        }
}

check_cor_accuracy_splatter <- function(sim, imputed) {
	require("Hmisc")
	gene_dir <- sign(rowData(sim)$DEFacGroup1 - rowData(sim)$DEFacGroup2)
	true_dir <- gene_dir %*% t(gene_dir)

	cor_out <- Hmisc::rcorr(t(imputed), type="spearman")
	thresh <- 0.05 / ( prod(dim(cor_out$P))/2 - length(diag(cor_out$P)) )

        TP <- sum(cor_out$P < thresh & true_dir != 0 &
		  sign(cor_out$r) == true_dir, na.rm=T)
        FP <- sum(cor_out$P < thresh, na.rm=T) - TP
        FN <- sum(cor_out$P >= thresh & true_dir != 0, na.rm=T)
        TN <- sum(cor_out$P >= thresh, na.rm=T) - FN
	vbad <- which( true_dir != 0 & sign(cor_out$r) != true_dir 
			& cor_out$P < thresh, arr.ind=T)
	vbad_p <- cor_out$P[vbad]
	vbad_r <- cor_out$r[vbad]

	FPs <- which( true_dir == 0 & cor_out$P < thresh, arr.ind=T )
	FPs_p <- cor_out$P[FPs]
	FPs_r <- cor_out$r[FPs]
	most_FP <- which(abs(FPs_r) >= quantile(abs(FPs_r), prob=1-100/length(FPs_r)))	

	vbad <- data.frame(row=vbad[,1], col=vbad[,2], p=vbad_p, r=vbad_r)
	FPs <- data.frame(row=FPs[,1], col=FPs[,2], p=FPs_p, r=FPs_r)

	return(list(stats=c(TP, FP, TN, FN), wrongway=vbad, mostfp=FPs[most_FP,]))
}

check_multiDE_accuracy_splatter <- function(sim, mat_name="logcounts", mag_centile=NULL) {
	require("Hmisc")
	require("CellTypeProfiles")
	require("scater")
	de_cols <- grep("DEFac", colnames(rowData(sim)))
	de_tab <- rowData(sim)[,de_cols]
	high <- apply(de_tab, 1, max)
	low <- apply(de_tab, 1, min)
	is.de <- (high-low) != 0
	colnames(de_tab) <- sub("DEFac", "", colnames(de_tab))
	de_rows <- which(is.de)
	# is.path?
	if (grepl("Path",sim$Group[1])){
		sim <- sim[,sim$Step > 50] # make paths more divergent
	}
	p_values <- apply(assays(sim)[[mat_name]], 1, function(x) { kruskal.test(x, factor(sim$Group))$p.value})
	p_values[is.na(p_values)] <- 1;
	q_values <- p.adjust(p_values, method="fdr")
	lvls <- my_row_mean_aggregate(assays(sim)[[mat_name]], factor(sim$Group))
	mag <- apply(lvls, 1, max) - apply(lvls, 1, min)

	# Check direction
	# Truth
	de_up <- high > 1
	top <- sapply(which(de_up), function(i){which(unlist(de_tab[i,]) == high[i])})
	de_dn <- low < 1
	bottom <- sapply(which(de_dn), function(i){which(unlist(de_tab[i,]) == low[i])})
	# Obs
	high_obs <- apply(lvls, 1, max)
	low_obs <- apply(lvls, 1, min)
	
	top_obs <- sapply(which(de_up), function(i){
				out <- which(unlist(lvls[i,]) == high_obs[i])
				if (length(out) > 1) {
					if (top[i] %in% out) { return(top[i])}
					else {return(out[1])}
				} else { return(out)}
				})
	bottom_obs <- sapply(which(de_dn), function(i){
				out <- which(unlist(lvls[i,]) == low_obs[i])
				if (length(out) > 1) {
					if (bottom[i] %in% out) {return(bottom[i])}
					else {return(out[1])}
				} else {return(out)}
				})
	# Deal with if both up & down
	dir.correct <- rep(TRUE, times=length(is.de))
	dir.correct[de_up] <- dir.correct[de_up] & top == top_obs
	dir.correct[de_dn] <- dir.correct[de_dn] & bottom == bottom_obs


	thresh <- 0.05 
	mag_thesh <- 0
	if (!is.null(mag_centile)) {
		mag_thresh <- quantile(mag, probs=1-mag_centile)
	}

        TP <- sum(q_values < thresh & is.de & dir.correct & mag > mag_thresh, na.rm=T)
        FP <- sum(q_values < thresh & mag > mag_thresh, na.rm=T) - TP
        FN <- sum(q_values >= thresh & is.de, na.rm=T)
        TN <- sum(q_values >= thresh, na.rm=T) - FN

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

	return(list(stats=c(TP, FP, TN, FN), wrongway=wrong_dir, mostfp=FPs))
}

