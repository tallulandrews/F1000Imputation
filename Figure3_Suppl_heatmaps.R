source("~/MAGIC/R_Imputation_Functions.R")

require("Matrix")

gene_cor_heatmap <- function(obj, mat="mat") {
	require("scater")
	this_mat <- as.matrix(obj@assays[[mat]])
	cor_mat <- cor(t(this_mat), method="spearman")
	require("gplots")
	require("RColorBrewer")
	heat_cols <- rev(brewer.pal(8, "RdBu"))
	bin_edges <- seq(from=-1, to=1, length=9)
	gene_info <- data.frame(is.up=rowData(obj)$g_up, is.down=rowData(obj)$g_down, expr=Matrix::rowMeans(this_mat))
	reorder <- order(gene_info[,1], gene_info[,2], gene_info[,3], decreasing=T)
	gene_info <- gene_info[reorder,]
	sideCols <- rep("grey85", nrow(obj))
	sideCols[gene_info[,1]] <- "red"
	sideCols[gene_info[,2]] <- "blue"
	cor_mat <- cor_mat[reorder, reorder]
	heatmap.2(cor_mat, col=heat_cols, trace="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=sideCols, RowSideColors=sideCols, breaks=bin_edges)
	invisible(cor_mat)
}
	

require("scater")
### Heatmap ###
obj_sce <- readRDS("Tongue_10X_permuted.rds") # Add autoencoders to this RDS separately

obj_sce <- obj_sce[, obj_sce$cell_type1 %in% obj_sce@metadata$type_pair]
source("~/MAGIC/CTP_Functions.R")
marks <- complex_markers(as.matrix(assays(obj_sce)[["counts"]]), factor(obj_sce$cell_type1), n_max=1)
marks <- marks[marks$AUC > 0.75 & marks$q.value < 0.05,]
set.seed(3819)
mark_sample <- sample(1:nrow(marks), min(250, nrow(marks)));
mark_sample <- which(rownames(obj_sce) %in% rownames(marks[mark_sample,]))

set.seed(3819)
subset <- sample(which(rowData(obj_sce)$Permuted), min(250, nrow(obj_sce))) #Subset genes to make visualization more legible
subset <- c(mark_sample, subset)
obj_sce <- obj_sce[subset, ]
rowData(obj_sce)$g_up <- rownames(obj_sce) %in% rownames(marks)[marks$Group == obj_sce@metadata$type_pair[1]]
rowData(obj_sce)$g_down <- rownames(obj_sce) %in% rownames(marks)[marks$Group == obj_sce@metadata$type_pair[2]]

png(paste("Tongue_10X_Raw_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "counts")
dev.off()
png(paste("Tongue_10X_RawLog_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "logcounts")
dev.off()
png(paste("Tongue_10X_MAGIC_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "magic")
dev.off()
png(paste("Tongue_10X_DrImpute_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "dri")
dev.off()
png(paste("Tongue_10X_scImpute_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "sci")
dev.off()
png(paste("Tongue_10X_SAVER_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "saver")
dev.off()
png(paste("Tongue_10X_DCA_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "dca")
dev.off()
png(paste("Tongue_10X_Knn_GeneCor_heatmap3.png", sep=""), width=4, height=4, units="in", res=300);
gene_cor_heatmap(obj_sce, "knn")
dev.off()

source("~/MAGIC/Colour_bar.R")
blank_plot <- function() {
        tmp <-  par("mar")
        par(mar=c(0,0,0,0))
        plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
        return(tmp);
}

heat_cols <- rev(brewer.pal(8, "RdBu"))
bin_edges <- seq(from=-1, to=1, length=9)

png("Tongue_GeneCor_heatmap3_colorbar.png", width=7, height=3, units="in", res=300)
blank_plot()
color.bar(heat_cols, min=-1, max=1, ticks.at=bin_edges, title="Correlation", horiz=T, add=T)
dev.off()

###
