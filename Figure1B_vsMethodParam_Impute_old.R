source("/nfs/users/nfs_t/ta6/MAGIC/Imputation_Simulation_Functions.R")
require("scater")

set.seed(20194)
nrep=10;
my_seeds <- round(runif(nrep*2)*10000)
my_seeds <- unique(my_seeds)
my_seeds <- my_seeds[1:nrep]
require("Rmagic") # t_diffusion
require("DrImpute") # dropout.probability.threshold
require("scImpute") # drop_thre
require("SAVER") # percent of genes to predict
source("/nfs/users/nfs_t/ta6/MAGIC/knn_smooth.R") # k
require(doParallel)
registerDoParallel(cores = 16)

# dropout.probability.threshold == drop_thre
args <- commandArgs(trailingOnly=TRUE)

method=args[1]

# wrapper functions
if (method == "scImpute") {
	param_name <- "Dropout Threshold"
	param_range=rev(c(0, 0.2, 0.4, 0.6, 0.8, 1))
	imputation_fxn <- function(sims, param) {
		registerDoParallel(cores = 16)
		saveRDS(sims, file="tmp.rds");
		scImpute::scimpute("./tmp.rds", infile="rds", outfile="rds",
			type="count", drop_thre=param, out_dir="./", 
			Kcluster=2, ncores=16)
		out <- readRDS("./scimpute_count.rds")
		return(out);
	}
} else if (method == "DrImpute") {
	param_name <- "Zeros Remaining"
	param_range=rev(c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1))
	imputation_fxn <- function(sims, param) {
		lognorm <- log2(sims+1)

		Tot <- ncol(lognorm) * nrow(lognorm)
		Z <- sum(lognorm==0)
		if (param >= Z/Tot) {param = (Z-1)/Tot}
		out <- DrImpute::DrImpute(lognorm, ks=2, zerop=param)
		return(as.matrix(out));
	}
} else if (method == "SAVER") {
	param_name <- "Percent of Genes"
	param_range=c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
	imputation_fxn <- function(sims, param) {
	   	if (param == 1) { 
		out <- saver(sims, do.fast=TRUE, size.factor=1)
		 } else {
		out <- saver(sims, do.fast=TRUE, size.factor=1, 
			npred=nrow(sims)*param)
		}
		return(out$estimate)
	}
} else if (method == "MAGIC") {
	param_name <- "Diffusion time"
	param_range=c(1, 2, 3, 4, 5, 6, 7, 8)
	imputation_fxn <- function(sims, param) {
		out <- Rmagic::run_magic(t(sims), 
				t_diffusion=param, lib_size_norm=F, 
				log_transform=F, pseudo_count=0.1, npca=100, 
				k=12, ka=4, epsilon=1, rescale_percent=0)
		return(t(out));
	}
} else if (method == "MAGIC_k") {
	param_name <- "K neighbours"
	param_range=c(5, 10, 15, 30, 45, 60, 100)
	imputation_fxn <- function(sims, param) {
		out <- Rmagic::run_magic(t(sims), 
				t_diffusion=3, lib_size_norm=F, 
				log_transform=F, pseudo_count=0.1, npca=100, 
				k=param, ka=4, epsilon=1, rescale_percent=0)
		return(t(out));
	}
} else if (method == "knn") {
	param_name <- "K neighbours"
	param_range=c(5, 10, 15, 30, 45, 60, 100)
	imputation_fxn <- function(sims, param) {
		out <- knn_smoothing(sims, k=param, d=10, seed=42)
		return(out)
	}
}

AllTP <- vector()
AllFP <- vector()
AllTN <- vector()
AllFN <- vector()

if (method %in% c("SAVER", "MAGIC", "knn", "MAGIC_k")) {
	param_range <- c(-1, param_range);
}

if(method=="scImpute" | method == "SAVER") {
	my_seeds <- my_seeds[1:3]
}

for (s in my_seeds) {
	sims <- readRDS(paste("vsMethodParams_seed_scVI_",s,".rds", sep=""))
	require("Hmisc")
	TP <- vector(); FP <- vector(); TN <- vector(); FN <- vector();
for (i in 1:length(param_range)) {
	val <- param_range[i]

	if (val == -1) {
		mat <- sims@assays[["counts"]]
	} else {
		mat <- imputation_fxn(sims@assays[["counts"]], val)
	}
	cor_out <- Hmisc::rcorr(t(mat), type="spearman")
	thresh <- 0.05 / ( prod(dim(cor_out$P))/2 - length(diag(cor_out$P)) )
	de <- rowData(sims)$g_up | rowData(sims)$g_down
	TP <- c(TP, sum(cor_out$P[de, de] < thresh, na.rm=T))
	FP <- c(FP, sum(cor_out$P < thresh, na.rm=T) - TP[i])
	FN <- c(FN, sum(cor_out$P[de, de] >= thresh, na.rm=T))
	TN <- c(TN, sum(cor_out$P >= thresh, na.rm=T) - FN[i])
}

AllTP <- cbind(AllTP, TP)
AllFP <- cbind(AllFP, FP)
AllTN <- cbind(AllTN, TN)
AllFN <- cbind(AllFN, FN)

}

FPR <- AllFP/(AllFP+AllTN)

out <- list(TP=AllTP, FP=AllFP, TN=AllTN, FN=AllFN, params=param_range, param_name=param_name)
saveRDS(out, file=paste(method,"vs_param.rds", sep="_"))

