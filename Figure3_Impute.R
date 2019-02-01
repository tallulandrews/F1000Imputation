require("SingleCellExperiment")
require("scater")
require("methods")
source("~/MAGIC/R_Imputation_Functions.R")

args <- commandArgs(trailingOnly=TRUE)
file=args[1]
n_cores=args[2]

#set.seed(9148)
obj <- readRDS(file)
obj$Group <- obj$cell_type1
#imputed <- Impute_default_all(obj, n.cores=n_cores)

resave=FALSE
if (!"dri" %in% names(obj@assays)) {
        res <- DrImpute_wrapper(obj)
        assays(obj)[["dri"]] <- res;
	resave=TRUE
}

if (!"magic" %in% names(obj@assays)) {
        res <- MAGIC_wrapper(obj)
        assays(obj)[["magic"]] <- res;
	resave=TRUE
}

if (!"knn" %in% names(obj@assays)) {
        res <- knn_wrapper(obj)
        assays(obj)[["knn"]] <- res;
	resave=TRUE
}
if (resave) {saveRDS(obj, file=file)}

if (!"saver" %in% names(obj@assays)) {
        res <- SAVER_wrapper(obj, n.cores=n_cores)
        assays(obj)[["saver"]] <- res;

saveRDS(obj, file=file)
}

if (!"sci" %in% names(obj@assays)) {
        res <- scImpute_wrapper(obj, n.cores=n_cores)
        assays(obj)[["sci"]] <- res;

saveRDS(obj, file=file)
}

