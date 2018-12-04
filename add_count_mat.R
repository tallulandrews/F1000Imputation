# args = SCE rds, matrix file, name of matrix.
args <- commandArgs(trailingOnly=TRUE)
require("scater")
require("methods")
sce <- readRDS(args[1])
check1 <- length(sce@assays);
check2 <- dim(sce)
mat <- read.table(args[2], header=T)
mat <- mat[, match(colnames(sce), colnames(mat))]
mat <- mat[match(rownames(sce), rownames(mat)), ]
mat[is.na(mat)] <- 0;
sce@assays[[ args[3] ]] <- mat

if (identical(dim(sce), check2) & length(sce@assays) > check1) {
	saveRDS(sce, file=args[1])
} else {
	print("Error: adding to SCE did not pass checks")
}

