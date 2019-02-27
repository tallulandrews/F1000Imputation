# Examples
require("scater")
prefix <- "Pancreas"

source("~/MAGIC/CTP_Functions.R")
obj <- readRDS(paste(prefix,"FACS_permuted.rds", sep="_"))
lvls <- my_row_mean_aggregate(assays(obj)[["logcounts"]], obj$cell_type1)

check_FPs <- function(mat, non_DE, nFPs=100) {
        mat <- mat[non_DE, ]
	tmp_lvls <- lvls[non_DE,]
        require("Hmisc")
        type1 <- obj@metadata$type_pair[1]
        type2 <- obj@metadata$type_pair[2]
        DE_p.value <- apply(mat, 1, function(x) {
                                wilcox.test(x[obj$cell_type1==type1],
                                            x[obj$cell_type1==type2])$p.value})
	DE_p.value[is.na(DE_p.value)] <- 1

	threshold <- quantile(DE_p.value, probs=nFPs/length(DE_p.value))
	top_genes <- rownames(mat)[DE_p.value <= threshold]
	top_ps <- DE_p.value[DE_p.value <= threshold]
	lvl1 <- tmp_lvls[DE_p.value <= threshold, type1]
	lvl2 <- tmp_lvls[DE_p.value <= threshold, type2]
	return(data.frame(Gene=top_genes, p.value=top_ps, expr1=lvl1, expr2=lvl2))
}

top_FP <- list()
non_DE <- rowData(obj)$Permuted
top_FP$lognorm <- check_FPs(assays(obj)[["logcounts"]], non_DE, nFPs=sum(non_DE))
top_FP$magic <- check_FPs(assays(obj)[["magic"]], non_DE, nFPs=sum(non_DE))
top_FP$drimpute <- check_FPs(assays(obj)[["dri"]], non_DE, nFPs=sum(non_DE))
top_FP$knn_sm <- check_FPs(assays(obj)[["knn"]], non_DE, nFPs=sum(non_DE))
top_FP$scimpute <- check_FPs(assays(obj)[["sci"]], non_DE, nFPs=sum(non_DE))
top_FP$saver <- check_FPs(assays(obj)[["saver"]], non_DE, nFPs=sum(non_DE))
top_FP$dca <- check_FPs(assays(obj)[["dca"]], non_DE, nFPs=sum(non_DE))

all_top <- c(as.character(top_FP$magic[,1]), as.character(top_FP$drimpute[,1]), as.character(top_FP$scimpute[,1]), as.character(top_FP$saver[,1]), as.character(top_FP$knn_sm[,1]), as.character(top_FP$dca[,1]))

permed1 <- check_FPs(assays(obj)[["logcounts"]], non_DE, nFPs=1000)
permed2 <- check_FPs(assays(obj)[["counts"]], non_DE, nFPs=1000)
permed1 <- permed1[permed1[,2] < 0.1,]
permed2 <- permed2[permed2[,2] < 0.1,]
exclude <- c(as.character(permed1[,1]), as.character(permed2[,1]))


#which(table(all_top) > 2)
top_FP$magic[ !top_FP$magic[,1] %in% exclude,]

#examples <- c("Baiap2", "Snx8", "Dusp23", "Mmp19", "Dync1li1") # Heart (but lowly expressed)
#examples <- c("Ubr5", "Slc6a17", "Gimap9", "Zfp606", "Rsf1") # Pancreas (but lowly expressed)
examples <- c("Rp9", "Ubr5", "Slc6a17",  "Zfp606", "Rsf1") # Pancreas (but lowly expressed)
examples <- c("Rars2", "Dpagt1", "Tars", "Pcbp3", "Akt1", "Zfp606")


type_pair <- obj@metadata$type_pair
meth_mat_names <- c("logcounts", "magic", "dri", "knn", "sci", "saver", "dca")
meth_good_names <- c("Unimputed", "MAGIC", "DrImpute", "knn", "scImpute", "SAVER", "dca")
cell_type_good_names <- c("PP cell", "A cell")
source("~/R-Scripts/violin_plot.R")
png("Figure6_Pancreas_Examples.png", width=1.5*length(examples), height=1.3*length(meth_mat_names), units="in", res=300)
par(mfrow=c(length(examples), length(meth_mat_names)))
for(gene in examples) {
	for(i in 1:length(meth_mat_names)) {
		margins <- c(0,0,0,0);
		p <- top_FP[[i]][rownames(top_FP[[i]])==gene,2]
		if (gene==examples[1]) {
			margins[3] = 1.5;
		}
		if (i == 1) {
			margins[2] = 1.5;
		}
		par(mar=margins);
		meth <- meth_mat_names[i]
		exp <- as.numeric(assays(obj)[[meth]][rownames(obj)==gene,])
		if (meth %in% c("sci", "saver", "magic", "dca", "knn")) {
			exp <- log2(exp+1)
		}
		vioplot(list(exp[obj$cell_type1==type_pair[1]], exp[obj$cell_type1==type_pair[2]]), 
			col=c("firebrick", "dodgerblue"), names=cell_type_good_names, xaxt=FALSE, yaxt=FALSE)
		if (gene==examples[1]) {
			title(main=meth_good_names[i], cex.main=2)
		}
		if (i==1) {
			title(ylab=gene, line=0, font=2, cex.lab=2)
		}
		if (p < 0.05/sum(non_DE)) {
			text(1.5,max(exp[obj$cell_type1 %in% type_pair]), "**", cex=2)
		} else if (p < 0.05) {
			text(1.5,max(exp[obj$cell_type1 %in% type_pair]), "*", cex=2)
		}
	}
}
dev.off()

