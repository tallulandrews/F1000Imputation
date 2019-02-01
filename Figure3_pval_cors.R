source("~/MAGIC/CTP_Functions.R")
source("~/MAGIC/Colour_Scheme.R")
require("scater")
library(methods)

# Matching files
file_pairs <- rbind(
	c("Bladder_10X_filtered.rds", "Bladder_FACS_filtered.rds"),
	c("Kidney_10X_filtered.rds", "Kidney_FACS_filtered.rds"),
	c("Lung_10X_filtered.rds", "Lung_FACS_filtered.rds"),
	c("Mammary_10X_filtered.rds", "Mammary_FACS_filtered.rds"),
	c("Marrow_10X_filtered.rds", "Marrow_FACS_filtered.rds"),
	c("Tongue_10X_filtered.rds", "Tongue_FACS_filtered.rds"),
	c("Muscle_10X_filtered.rds", "Muscle_FACS_filtered.rds"))
rownames(file_pairs) <- c("Bladder", "Kidney", "Lung", "Mammary", "Marrow", "Tongue", "Muscle")

p_val_cors <- matrix(0, nrow=nrow(file_pairs), ncol=length(Method_Legend))
p_val_cors_n <- matrix(0, nrow=nrow(file_pairs), ncol=length(Method_Legend))
rownames(p_val_cors) <- rownames(file_pairs);

for (f in 1:nrow(file_pairs)) {
	drop <- readRDS(file_pairs[f,1])
	facs <- readRDS(file_pairs[f,2])

	# Only cell-types that exist in both datasets
	types <- intersect( as.character(unique(drop$cell_type1)), as.character(unique(facs$cell_type1)) )
	drop <- drop[, drop$cell_type1 %in% types]
	facs <- facs[, facs$cell_type1 %in% types]

	genes <- intersect(rownames(drop), rownames(facs));
	drop <- drop[match(genes, rownames(drop)),]
	facs <- facs[match(genes, rownames(facs)),]
	# Exclude low detection genes
	exclude <- Matrix::rowSums(drop@assays[["counts"]] > 0) < ncol(drop@assays[["counts"]])*0.01 | 
		   Matrix::rowSums(facs@assays[["counts"]] > 0) < ncol(facs@assays[["counts"]])*0.01
	drop <- drop[!exclude,]
	facs <- facs[!exclude,]

	# Setup output
	methods <- names(assays(drop));
	methods_nice_names <- Method_Colours[match(methods, Method_Colours[,1]),2]
	methods_nice_names <- names(Method_Legend)[match(methods_nice_names, Method_Legend)]
	methods <- methods[!duplicated(methods_nice_names)]
	methods_nice_names <- methods_nice_names[!duplicated(methods_nice_names)]

	tissue <- rownames(file_pairs)[f]

	# Identify Markers
	for (m in 1:length(methods)) {
		meth <- methods[m];
		ctypes <- as.character(unique(drop$cell_type1))

		for (i in 1:(length(ctypes))) {
#			for (j in (i+1):length(ctypes)) {
				p_facs <- apply(as.matrix(assays(facs)[[meth]]), 1, function(x){wilcox.test(x[as.character(facs$cell_type1) == ctypes[i]], x[as.character(facs$cell_type1) != ctypes[i]], alternative="greater")$p.value})
				p_drop <- apply(as.matrix(assays(drop)[[meth]]), 1, function(x){wilcox.test(x[as.character(drop$cell_type1) == ctypes[i]], x[as.character(drop$cell_type1) != ctypes[i]], alternative="greater")$p.value})
				p_facs[!is.finite(p_facs)] <- 0.5
				p_drop[!is.finite(p_drop)] <- 0.5
				c <- cor(p_facs, p_drop, method="spearman")
				p_val_cors[f,m] <- p_val_cors[f,m]+c
				p_val_cors_n[f,m] <- p_val_cors_n[f,m]+1
				colnames(p_val_cors) <- remap_names(methods)
				colnames(p_val_cors_n) <- remap_names(methods)
#			}
		}

	}

}
p_val_cors <- p_val_cors/p_val_cors_n
saveRDS(p_val_cors, file="TM_mark_cors.rds") 

### Correlate p-values ###
require("gplots")
require("RColorBrewer")
p_val_cors <- readRDS("TM_mark_cors.rds")

png("Pvalue_Cor_heatmap.png", width=5, height=5, units="in", res=300)
heatmap.2(abs(t(p_val_cors)), trace="none", cellnote=t(round(p_val_cors, digits=2)), col=brewer.pal(8, "Blues"), notecol="black", margins=c(6,6), key=FALSE, dendrogram="none", breaks=9)
dev.off()


