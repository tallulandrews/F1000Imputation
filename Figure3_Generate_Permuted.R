#Set-up
require("SingleCellExperiment")
require("scater")
require("Matrix")
source("~/MAGIC/CTP_Functions.R")
source("~/MAGIC/R_Imputation_Functions.R")

# Read SingleCellExperiment object from RDS passed as argument
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
prefix <- unlist(strsplit(file, "[./]"))
prefix <- prefix[length(prefix)-1]

if (!file.exists(paste(prefix,"permuted.rds", sep="_"))) {


set.seed(1928)
obj <- readRDS(file)

# Filter dataset
types <- table(factor(obj$cell_type1));
types <- names(types)[types/sum(types) > 0.05] # 5% of cells in cell-type
obj <- obj[,obj$cell_type1 %in% types & obj$cell_type1 != "unknown"] # named cell-type
obj <- obj[rowSums(counts(obj) > 0) > 0.05*ncol(obj),] # detect in 5% of cells
sf <- colSums(assays(obj)[["counts"]])
norm <- t( t(assays(obj)[["counts"]] / sf*median(sf) ) )
norm <- log2(norm +1)
assays(obj)[["logcounts"]] <- norm;
colData(obj)$Size_Factor <- sf;
obj$cell_type1 <- factor(obj$cell_type1);

# Select pair of cell-types
profiles <- my_row_mean_aggregate(logcounts(obj), obj$cell_type1)
pair_dist <- as.matrix(dist(t(profiles)))
diag(pair_dist) <- NA
pair <- which(pair_dist == min(pair_dist, na.rm=T), arr.ind=T)
pair <- rownames(pair)

obj@metadata$type_pair <- pair

#Test for DE
type1 <- obj@metadata$type_pair[1]
type2 <- obj@metadata$type_pair[2]
rowData(obj)$DE_p.value <- apply(assays(obj)[["logcounts"]], 1, function(x) {
				wilcox.test(x[obj$cell_type1==type1], 
					    x[obj$cell_type1==type2])$p.value})
rowData(obj)$DE_q.value <- p.adjust(rowData(obj)$DE_p.value, method="fdr")

# Permute logcounts
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
}

whole_permuted <- assays(obj)[["logcounts"]]
whole_permuted[non_DE,ab_cells] <- data

# revert to counts
whole_permuted_counts <- 2^whole_permuted -1
whole_permuted_counts <- round(t(t(whole_permuted_counts) * sf/median(sf)))

# save in object

saveRDS(obj, file=paste(prefix,"filtered.rds", sep="_"))

rowData(obj)$Permuted <- non_DE
assays(obj)[["logcounts"]] <- whole_permuted
assays(obj)[["counts"]] <- whole_permuted_counts

saveRDS(obj, file=paste(prefix,"permuted.rds", sep="_"))
}
