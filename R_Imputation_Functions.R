### This is a Wrapper for imputation functions to be used in other scripts. ###
# Impute_default_all is used in analyzing real and splatter simulated data.
# Individual wrappers are used for testing effects of method parameters.

# Format required to add additional imputation methods:
#	input: a SingleCellExperiment object with "counts" and "logcounts"
#		a parameter value (with default), 
#		whether the input should be normalized (with default=TRUE)
#		number of cores to use (if applicable)
#		random seed to ensure reproducibility (with default).
#	output: the imputed matrix where rows=genes, cols=cells

require("scater")
require("Rmagic") 
require("DrImpute") 
require("scImpute") 
require("SAVER") 
source("/nfs/users/nfs_t/ta6/MAGIC/knn_smooth.R") 


scImpute_wrapper<- function(sce, param=0.5, do.norm=TRUE, n.cores=16, seed=42) {
	# normalization is never used
	tmp_dir <- paste(paste("Tmp",round(runif(1)*100000), sep="_"),"/",sep="")
	dir.create(tmp_dir)
	set.seed(seed)
	require(doParallel)
	cl <- makeCluster(n.cores)
	registerDoParallel(cl)
	param_name <- "Dropout Threshold"
	registerDoParallel(cores = n.cores)
        saveRDS(assays(sce)[["counts"]], file=paste(tmp_dir,"tmp.rds", sep=""));
        scImpute::scimpute(paste(tmp_dir,"tmp.rds", sep=""), infile="rds", outfile="rds",
                type="count", drop_thre=param, out_dir=tmp_dir,
                Kcluster=length(unique(sce$Group)), ncores=n.cores)
        out <- readRDS(paste(tmp_dir,"scimpute_count.rds", sep=""))
	stopCluster(cl)
	unlink(tmp_dir, recursive=TRUE)
        return(out);
}

DrImpute_wrapper <- function(sce, param=0, do.norm=TRUE, seed=42) {
	# uses pre-defined log-normalized matrix
	set.seed(seed)
	param_name <- "Zeros Remaining"
	out <- DrImpute::DrImpute(assays(sce)[["logcounts"]],
                  ks=length(unique(sce$Group)),
                  zerop=param)
	return(out)
}

SAVER_wrapper <- function(sce, param=1, do.norm=TRUE, seed=42, n.cores=16){
	# optional CPM-like normalization
	set.seed(seed)
	require(doParallel)
	cl <- makeCluster(n.cores)
	registerDoParallel(cl)
	param_name <- "Percent of Genes"
	sf <- 1
	if (do.norm) {
		sf <- NULL
	}
	if (param == 1) {
                out <- saver(assays(sce)[["counts"]], do.fast=TRUE, size.factor=sf, ncores=n.cores)
        } else {
                out <- saver(assays(sce)[["counts"]], do.fast=TRUE, size.factor=sf,
                        npred=nrow(sce)*param, ncores=n.cores)
        }
	stopCluster(cl)
	return(out$estimate);
}

MAGIC_k_wrapper <- function(sce, param=12, do.norm=TRUE, seed=42) {
	# optional CPM-like normalization
	set.seed(seed)
	param_name <- "K neighbours"
	out <- Rmagic::run_magic(t(assays(sce)[["counts"]]), 
		t_diffusion=3, lib_size_norm=do.norm, 
		log_transform=F, k=param)
	return(t(out))
}

MAGIC_wrapper <- function(sce, param=0, do.norm=TRUE, seed=42) {
	# optional CPM-like normalization
	set.seed(seed)
	param_name <- "Diffusion time"
	out <- Rmagic::run_magic(t(assays(sce)[["counts"]]), 
		t_diffusion=param, lib_size_norm=do.norm, 
		log_transform=F, k=12)
	out <- t(out)
	rownames(out) <- rownames(sce)
	colnames(out) <- colnames(sce)
	return(out)
}

knn_wrapper <- function(sce, param=ncol(sce)/20, do.norm=TRUE, seed=42) {
	# includes CPM-like normalization
	param_name <- "K neighbours"
	out <- knn_smoothing(assays(sce)[["counts"]], k=param, d=10, seed=seed)
        return(out)
}

Impute_default_all <- function(sce, n.cores=16) {
	res <- scImpute_wrapper(sce, n.cores=n.cores)
	assays(sce)[["sci"]] <- res;

	res <- DrImpute_wrapper(sce)
	assays(sce)[["dri"]] <- res;

	res <- MAGIC_wrapper(sce)
	assays(sce)[["magic"]] <- res;

	res <- knn_wrapper(sce)
	assays(sce)[["knn"]] <- res;

	res <- SAVER_wrapper(sce, n.cores=n.cores)
	assays(sce)[["saver"]] <- res;

	return(sce)
}
