source("/nfs/users/nfs_t/ta6/MAGIC/Splatter_Functions.R")
require("scater")

# args: rds of simulated & imputed results, name of output
args <- commandArgs(trailingOnly=TRUE)
infile=args[1]
outfile=args[2]

imputation_methods=c("MAGIC", "scImpute", "DrImpute", "SAVER", "knn", "dca", "scVI")

# Analyze DE in simulated datasets:
accuracy <- list()
sim <- readRDS(infile)
for (meth in c("counts", "logcounts", imputation_methods)) {
	out <- check_multiDE_accuracy_splatter(sim, mat_name=meth);
	accuracy[[meth]] <- out;
}
saveRDS(accuracy, file=outfile)

