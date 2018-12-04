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
	cor_mat <- cor(as.matrix(t(sim[[mat]])), method="spearman")
	require("gplots")
	require("RColorBrewer")
	heat_cols <- rev(brewer.pal(8, "RdBu"))
	bin_edges <- seq(from=-1, to=1, length=9)
	gene_info <- data.frame(is.up=sim$g_up, is.down=sim$g_down, expr=rowMeans(sim[[mat]]))
	reorder <- order(gene_info[,1], gene_info[,2], gene_info[,3], decreasing=T)
	gene_info <- gene_info[reorder,]
	sideCols <- rep("grey85", nrow(sim$mat))
	sideCols[gene_info[,1]] <- "red"
	sideCols[gene_info[,2]] <- "blue"
	cor_mat <- cor_mat[reorder, reorder]
	heatmap.2(cor_mat, col=heat_cols, trace="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=sideCols, RowSideColors=sideCols, breaks=bin_edges)
	invisible(cor_mat)
}
	


### Generate Simulated Data ###
set.seed(20194)
sims <- Sim(n_genes=5000, n_cells=1000, max_mean=4, min_mean=-3, propDE=0.5)
sims$counts <- sims$mat;
sf <- colSums(sims$mat)
sims$norm <- t(t(sims$mat)/sf*median(sf))
sims$lognorm <- log2(sims$norm+1)

### Impute it. ###
require("Rmagic")
require("DrImpute")
require("scImpute")
require("SAVER")
source("/nfs/users/nfs_t/ta6/MAGIC/knn_smooth.R")
set.seed(28198)

# Magic
sims$magic <- t(Rmagic::run_magic(t(sims$counts), t_diffusion=0, lib_size_norm=F, log_transform=F, pseudo_count=0.1, npca=100, k=12, ka=4, epsilon=1, rescale_percent=0))

# DrImpute
sims$drimpute <- DrImpute::DrImpute(sims$lognorm, ks=2, zerop=0)

# scImpute
saveRDS(sims$counts, file="tmp.rds")
scImpute::scimpute("./tmp.rds", infile="rds", outfile="rds", type="count", drop_thre=0, out_dir="./", Kcluster=2,  ncores=16)
sims$scimpute <- readRDS("./scimpute_count.rds")

# kNN-smooth
sims$knn_sm <- knn_smoothing(sims$counts, k=50, d=10, seed=42)

# SAVER
require(doParallel)
registerDoParallel(cores = 16)
saver_out <- saver(sims$counts, do.fast=FALSE, size.factor=1)
sims$saver <- saver_out$estimate

saveRDS(sims, file="Impute_TwoClust_Example_Big.rds") # Add autoencoders to this RDS separately

### Heatmap ###
set.seed(3819)
subset <- sample(1:nrow(sims$saver), min(500, nrow(sims$saver))) #Subset genes to make visualization more legible
sims$mat <- sims$mat[subset,]
sims$counts <- sims$counts[subset,]
sims$norm <- sims$norm[subset,]
sims$lognorm <- sims$lognorm[subset,]
sims$magic <- sims$magic[subset,]
sims$drimpute <- sims$drimpute[subset,]
sims$scimpute <- sims$scimpute[subset,]
sims$knn_sm <- sims$knn_sm[subset,]
sims$saver <- sims$saver[subset,]
sims$g_up <- sims$g_up[subset]
sims$g_down <- sims$g_down[subset]

png(paste("Raw_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "counts")
dev.off()
png(paste("RawLog_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "lognorm")
dev.off()
png(paste("MAGIC_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "magic")
dev.off()
png(paste("DrImpute_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "drimpute")
dev.off()
png(paste("scImpute_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "scimpute")
dev.off()
png(paste("SAVER_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "saver")
dev.off()
png(paste("DCA_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "dca")
dev.off()
png(paste("scVI_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "scVI")
dev.off()
png(paste("Knn_GeneCor_heatmap2.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(sims, "knn_sm")
dev.off()

source("Colour_bar.R")
heat_cols <- rev(brewer.pal(8, "RdBu"))
bin_edges <- seq(from=-1, to=1, length=9)

png("GeneCor_heatmap2_colorbar.png", width=7, height=3, units="in", res=300)
blank_plot()
color.bar(heat_cols, min=-1, max=1, ticks.at=bin_edges, title="Correlation", horiz=T, add=T)
dev.off()

###
