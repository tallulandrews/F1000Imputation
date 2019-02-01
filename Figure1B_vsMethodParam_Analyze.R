require("scater")

set.seed(20194)
nrep=10;
my_seeds <- round(runif(nrep*2)*10000)
my_seeds <- unique(my_seeds)
my_seeds <- my_seeds[1:nrep]

#args <- commandArgs(trailingOnly=TRUE)
#method=args[1]
method <- "scVI"


AllTP <- vector();
AllFP <- vector();
AllTN <- vector();
AllFN <- vector();

files <- Sys.glob("vsMethodParams_seed_scVI*.rds")
for (f in files) {
	sims <- readRDS(f)
	require("Hmisc")
	mats <- grep(method, names(sims@assays));
	val_range <- c();
	RES <- vector();
for (i in c(-1, 1:length(mats))) {
		
	if (i == -1) {
		m <- sims@assays[["counts"]]
		val <- "Raw"

	} else {
		val <- names(sims@assays)[mats[i]];
		val <- sub(method, "", val)
		val <- sub("_", "", val)
		m <- sims@assays[[mats[i]]]
		if (val=="75"){next;}
	}
	cor_out <- Hmisc::rcorr(t(m), type="spearman")
	thresh <- 0.05 / ( prod(dim(cor_out$P))/2 - length(diag(cor_out$P)) )
	de <- rowData(sims)$g_up | rowData(sims)$g_down
	TP <- sum(cor_out$P[de, de] < thresh, na.rm=T)
	FP <- sum(cor_out$P < thresh, na.rm=T) - TP
	FN <- sum(cor_out$P[de, de] >= thresh, na.rm=T)
	TN <- sum(cor_out$P >= thresh, na.rm=T) - FN
	RES <- rbind(RES, c(TP, FP, FN, TN));
	val_range <- c(val_range, val);
}

AllTP <- cbind(AllTP, RES[,1])
AllFP <- cbind(AllFP, RES[,2])
AllTN <- cbind(AllTN, RES[,4])
AllFN <- cbind(AllFN, RES[,3])

}

FPR <- AllFP/(AllFP+AllTN)

out <- list(TP=AllTP, FP=AllFP, TN=AllTN, FN=AllFN, params=val_range, param_name="Train Percent", FPR=FPR)
saveRDS(out, file=paste(method,"vs_param.rds", sep="_"))

