source("/nfs/users/nfs_t/ta6/MAGIC/R_Imputation_Functions.R")
require("scater")

set.seed(20194)
nrep=10;
my_seeds <- round(runif(nrep*2)*10000)
my_seeds <- unique(my_seeds)
my_seeds <- my_seeds[1:nrep]

args <- commandArgs(trailingOnly=TRUE)

method=args[1]

# wrapper functions
if (method == "scImpute") {
	my_seeds <- my_seeds[1:3]
	param_name <- "Dropout Threshold"
	param_range=rev(c(0, 0.2, 0.4, 0.6, 0.8, 1))
	imputation_fxn <- scImpute_wrapper
} else if (method == "DrImpute") {
	param_name <- "Zeros Remaining"
	param_range=rev(c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1))
	imputation_fxn <- DrImpute_wrapper
} else if (method == "SAVER") {
	my_seeds <- my_seeds[1:3]
	param_name <- "Percent of Genes"
	param_range=c(-1, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
	imputation_fxn <- SAVER_wrapper
} else if (method == "MAGIC") {
	param_name <- "Diffusion time"
	param_range=c(-1, 1, 2, 3, 4, 5, 6, 7, 8)
	imputation_fxn <- MAGIC_wrapper
} else if (method == "MAGIC_k") {
	param_name <- "K neighbours"
	param_range=c(-1, 5, 10, 15, 30, 45, 60, 100)
	imputation_fxn <- MAGIC_k_wrapper
} else if (method == "knn") {
	param_name <- "K neighbours"
	param_range=c(-1, 5, 10, 15, 30, 45, 60, 100)
	imputation_fxn <- knn_wrapper
}

AllTP <- vector()
AllFP <- vector()
AllTN <- vector()
AllFN <- vector()

for (s in my_seeds) {
	sims <- readRDS(paste("vsMethodParams_seed_",s,".rds", sep=""))
	require("Hmisc")
	TP <- vector(); FP <- vector(); TN <- vector(); FN <- vector();
for (i in 1:length(param_range)) {
	val <- param_range[i]

	if (val == -1) {
		mat <- sims@assays[["counts"]]
	} else {
		mat <- imputation_fxn(sims@assays[["counts"]], param=val)
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

