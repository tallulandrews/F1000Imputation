#Set-up
require("SingleCellExperiment")
require("scater")
require("Matrix")
require("CellTypeProfiles")

n.cores=20

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
prefix <- unlist(strsplit(file, "[./]"))
prefix <- prefix[length(prefix)-1]

if (!file.exists(paste(prefix,"permuted.rds", sep="_"))) {
set.seed(1928)

obj <- readRDS(file)

types <- table(factor(obj$cell_type1));
types <- names(types)[types/sum(types) > 0.05] # 5% of cells in cell-type
obj <- obj[,obj$cell_type1 %in% types & obj$cell_type1 != "unknown"]
obj <- obj[rowSums(counts(obj) > 0) > 0.05*ncol(obj),] # detect in 5% of cells
obj <- normalise(obj)
sf <- colSums(assays(obj)[["counts"]])
norm <- t( t(assays(obj)[["counts"]] / sf*median(sf) ) )
norm <- log2(norm +1)
assays(obj)[["logcounts"]] <- norm;
colData(obj)$Size_Factor <- sf;
obj$cell_type1 <- factor(obj$cell_type1);

profiles <- my_row_mean_aggregate(logcounts(obj), obj$cell_type1)
pair_dist <- as.matrix(dist(t(profiles)))
diag(pair_dist) <- NA
pair <- which(pair_dist == min(pair_dist, na.rm=T), arr.ind=T)
pair <- rownames(pair)

obj@metadata$type_pair <- pair

#type1 vs type2
type1 <- obj@metadata$type_pair[1]
type2 <- obj@metadata$type_pair[2]
rowData(obj)$DE_p.value <- apply(assays(obj)[["logcounts"]], 1, function(x) {
				wilcox.test(x[obj$cell_type1==type1], 
					    x[obj$cell_type1==type2])$p.value})
rowData(obj)$DE_q.value <- p.adjust(rowData(obj)$DE_p.value, method="fdr")

# type1-type2 permutation
set.seed(32891) # 1
#set.seed(1046) # 2
require("scater")
require("Matrix")
require("permute")
non_DE <- rowData(obj)$DE_p.value > 0.2
non_DE[is.na(non_DE)] <- FALSE
type1 <- obj@metadata$type_pair[1]
type2 <- obj@metadata$type_pair[2]

ab_cells <- which(colData(obj)$cell_type1 == type1 | colData(obj)$cell_type1 == type2)

data <- assays(obj)[["logcounts"]][non_DE,ab_cells]
data_counts <- assays(obj)[["counts"]][non_DE,ab_cells]
shuffles <- shuffleSet( ncol(data), nset = nrow(data) )
for(i in 1:nrow(data)) {
	data[i,] <- data[i,shuffles[i,]]
#	data_counts[i,] <- data_counts[i,shuffles[i,]]
}

whole_permuted <- assays(obj)[["logcounts"]]
#whole_permuted_counts <- assays(obj)[["counts"]]
whole_permuted[non_DE,ab_cells] <- data
#whole_permuted_counts[non_DE,ab_cells] <- data_counts

sf <- colSums(assays(obj)[["counts"]])
whole_permuted_counts <- 2^whole_permuted -1
whole_permuted_counts <- round(t(t(whole_permuted_counts) * sf/median(sf)))


rowData(obj)$Permuted <- non_DE

assays(obj)[["perm_lognorm_1"]] <- whole_permuted
assays(obj)[["perm_counts_1"]] <- whole_permuted_counts

saveRDS(obj, file=paste(prefix,"permuted.rds", sep="_"))

}

if (!file.exists(paste(prefix,"imputed.rds", sep="_"))) {
set.seed(9148)

#### Imputation ####
obj <- readRDS(paste(prefix,"permuted.rds", sep="_"))

require(doParallel)
registerDoParallel(cores = n.cores)
require(SingleCellExperiment)
set.seed(2719) # I forgot this when I did it.

require("Rmagic")
whole_magic <- t(Rmagic::run_magic(t(assays(obj)[["perm_lognorm_1"]]), t_diffusion=0, 
			lib_size_norm=F, log_transform=F, pseudo_count=0.1, npca=20, 
			k=12, ka=4, epsilon=1, rescale_percent=0))
colnames(whole_magic) <- colnames(obj)
rownames(whole_magic) <- rownames(obj)
assays(obj)[["magic_1"]] <- whole_magic
# magic 1, t=4
# magic 2, t=3

require(doParallel)
registerDoParallel(cores = n.cores)

require("DrImpute")
whole_drimpute <- DrImpute::DrImpute(assays(obj)[["perm_lognorm_1"]], ks=length(unique(obj$cell_type1)), dropout.probability.threshold=0, fast=TRUE, mc.cores=n.cores)
colnames(whole_drimpute) <- colnames(obj)
rownames(whole_drimpute) <- rownames(obj)
assays(obj)[["dri_1"]] <- whole_drimpute

require(doParallel)
registerDoParallel(cores = n.cores)

source("/nfs/users/nfs_t/ta6/MAGIC/knn_smooth.R")
whole_knn <- knn_smoothing(assays(obj)[["perm_lognorm_1"]], k=50, d=10, seed=42)
colnames(whole_knn) <- colnames(obj)
rownames(whole_knn) <- rownames(obj)
assays(obj)[["knn_1"]] <- whole_knn

require(doParallel)
registerDoParallel(cores = n.cores)

require("scImpute")
saveRDS(assays(obj)[["perm_counts_1"]], file=paste(prefix,"tmp2.rds", sep="_")) 
scImpute::scimpute(paste(prefix,"tmp2.rds", sep="_"), infile="rds", outfile="rds", type="count", drop_thre=0, out_dir=prefix, Kcluster=length(unique(obj$cell_type1)),  ncores=n.cores, labeled=TRUE, labels=as.character(obj$cell_type1))
whole_sci <- readRDS(paste(prefix,"scimpute_count.rds", sep=""))
colnames(whole_sci) <- colnames(obj)
rownames(whole_sci) <- rownames(obj)
assays(obj)[["sci_1"]] <- whole_sci

require(doParallel)
registerDoParallel(cores = n.cores)

require("SAVER")
whole_saver <- saver(assays(obj)[["perm_counts_1"]], do.fast=TRUE, pred.genes=which(rowData(obj)$Permuted ), size.factor=1)$estimate # to normalize or not to normalize?
#saver1.cor.gene <- cor.genes(saver1)
colnames(whole_saver) <- colnames(obj)
rownames(whole_saver) <- rownames(obj)
assays(obj)[["saver_1"]] <- whole_saver

saveRDS(obj, file=paste(prefix,"imputed.rds", sep="_"))
}

if (!file.exists(paste(prefix,"FPR.rds", sep="_"))) {
set.seed(4817)

# Summary Statistics
require("scater")

obj <- readRDS(paste(prefix,"imputed.rds", sep="_"))
check_FPs <- function(mat, non_DE) {
	mat <- mat[non_DE, ]
	require("Hmisc")
	type1 <- obj@metadata$type_pair[1]
	type2 <- obj@metadata$type_pair[2]
	DE_p.value <- apply(mat, 1, function(x) {
				wilcox.test(x[obj$cell_type1==type1], 
					    x[obj$cell_type1==type2])$p.value})

	threshold <- 0.05/length(non_DE)
	FPs <- sum( DE_p.value < threshold, na.rm=T);
	FPR <- FPs/nrow(mat)
	return(FPR)
}

FPR <- list()
non_DE <- rowData(obj)$Permuted
FPR$magic <- check_FPs(assays(obj)[["magic_1"]], non_DE)
FPR$drimpute <- check_FPs(assays(obj)[["dri_1"]], non_DE)
FPR$scimpute <- check_FPs(assays(obj)[["sci_1"]], non_DE)
FPR$saver <- check_FPs(assays(obj)[["saver_1"]], non_DE)
FPR$knn_sm <- check_FPs(assays(obj)[["knn_1"]], non_DE)

saveRDS(FPR, paste(prefix,"FPR.rds", sep="_"))
}

