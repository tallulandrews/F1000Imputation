source("~/MAGIC/Bulk_DE_Tests.R")

get_TP <- function(de_list, sig_thresh=0.05) {
	consensus <- vector()
	for (DE in de_list) {
		DE$q.value[is.na(DE$q.value)] <- 1
		sig <- DE[DE$q.value < sig_thresh,]
		groups <- matrix(unlist(strsplit(as.character(sig$test), "_vs_")), ncol=2, byrow=TRUE)
		up <- which(sig$Mean1 > sig$Mean2)
		down <- which(sig$Mean1 < sig$Mean2)
		gene_dir <- as.character(sig$Gene)
		gene_dir[up] <- paste(gene_dir[up], groups[up,1], groups[up,2])
		gene_dir[down] <- paste(gene_dir[down], groups[down,2], groups[down,1])
		if (length(consensus)==0) {
			consensus <- gene_dir
		} else {
			consensus <- consensus[consensus %in% gene_dir]
		}
	}
	return(consensus)
}

get_TN <- function(de_list, sig_threshold=0.2) {
	exclude <- vector()
	genes <- vector();
	for (DE in de_list) {
		genes <- c(genes, DE$Gene);
		DE$q.value[is.na(DE$q.value)] <- 1
		sig <- DE[DE$q.value < sig_thresh,]
		gene <- sig[,1]
		exclude <- c(exclude, as.character(gene))
	}
	genes <- unique(genes)
	genes <- genes[! genes %in% exclude]
	return(genes)
}
	

# Kolo Bulk
mat_kolo <- read.table("/lustre/scratch117/cellgen/team218/TA/Simulations_Temporary_Files/Imputation/vsBulk/Kolo_Bulk.txt")
mat_kolo <- mat_kolo[,-1]
mat_kolo <- mat_kolo[rowSums(mat_kolo > 5) > 2,]
mat_kolo <- mat_kolo[-grep("^_", rownames(mat_kolo)),]
labels_kolo <- c("a2i", "a2i", "lif", "lif")
edgeR_DE <- do_edgeR(mat_kolo, labels_kolo, fdr_thresh=2)
DESeq2_DE <- do_DESeq2(mat_kolo, labels_kolo, fdr_thresh=2)

kolo_TP <- get_TP(list(edgeR_DE, DESeq2_DE))
kolo_TN <- get_TN(list(edgeR_DE, DESeq2_DE))

kolo_TP <- matrix(unlist(strsplit(kolo_TP, " ")), ncol=3, byrow=TRUE)
saveRDS(list(DE=kolo_TP, nonDE=kolo_TN), "Kolo_Bulk_Truth.rds")


# Tung Bulk
mat_tung <- read.table("/lustre/scratch117/cellgen/team218/TA/Simulations_Temporary_Files/Imputation/vsBulk/GSE77288_reads-raw-bulk-per-sample.txt.gz", header=TRUE, stringsAsFactors=FALSE)
mat_tung <- t(mat_tung)
ann <- mat_tung[1:3,]
mat_tung <- mat_tung[-c(1:3),]
genes <- rownames(mat_tung)
mat_tung <- matrix(as.numeric(unlist(mat_tung)), nrow=nrow(mat_tung))
rownames(mat_tung) <- genes;
mat_tung <- mat_tung[rowSums(mat_tung > 5) > 2,]
labels_tung <- factor(ann[1,])

edgeR_DE <- do_edgeR(mat_tung, labels_tung, fdr_thresh=2)
DESeq2_DE <- do_DESeq2(mat_tung, labels_tung, fdr_thresh=2)

tung_TP <- get_TP(list(edgeR_DE, DESeq2_DE))
tung_TN <- get_TN(list(edgeR_DE, DESeq2_DE))

tung_TP <- matrix(unlist(strsplit(tung_TP, " ")), ncol=3, byrow=TRUE)
saveRDS(list(DE=tung_TP, nonDE=tung_TN), "Tung_Bulk_Truth.rds")
