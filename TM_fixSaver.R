args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
require("SingleCellExperiment")
require("scater")

#obj <- readRDS(file=paste(prefix,"imputed.rds", sep="_"))

#n.cores=16
#require(doParallel)
#registerDoParallel(cores = n.cores)

#require("SAVER")
#whole_saver <- saver(assays(obj)[["perm_counts_1"]], do.fast=TRUE, size.factor=1)$estimate # to normalize or not to normalize?
##saver1.cor.gene <- cor.genes(saver1)
#colnames(whole_saver) <- colnames(obj)
#rownames(whole_saver) <- rownames(obj)
#assays(obj)[["saver_1"]] <- whole_saver

#saveRDS(obj, file=paste(prefix,"imputed2.rds", sep="_"))


set.seed(38982)
obj <- readRDS(file=prefix)
require("DrImpute")
out <- DrImpute::DrImpute(as.matrix(assays(obj)[["perm_lognorm_1"]]), ks=length(unique(obj$cell_type1)), zerop=0)
colnames(out) <- colnames(obj)
rownames(out) <- rownames(obj)
assays(obj)[["dri_1"]] <- out
saveRDS(obj, file=prefix)
