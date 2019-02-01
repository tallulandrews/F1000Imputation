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

# Command line argument: AUC threshold 
#args <- commandArgs(trailingOnly=TRUE)
#auc_threshold <- as.numeric(args[1])
p_val_cors <- matrix(0, nrow=nrow(file_pairs), ncol=length(Method_Legend))
rownames(p_val_cors) <- rownames(file_pairs);

cross_dataset_outfile <- paste("TM", "Data_Consistency.rds", sep="_");
cross_method_outfile <- paste("TM", "Method_Consistency.rds", sep="_");

allDataOUT <- vector();
allMethOUT <- vector();
for (i in 1:nrow(file_pairs)) {
	drop <- readRDS(file_pairs[i,1])
	facs <- readRDS(file_pairs[i,2])

	# Only cell-types that exist in both datasets
	types <- intersect( as.character(unique(drop$cell_type1)), as.character(unique(facs$cell_type1)) )
	drop <- drop[, drop$cell_type1 %in% types]
	facs <- facs[, facs$cell_type1 %in% types]

	# Exclude any permuted genes if present
#	if ("Permuted" %in% colnames(rowData(drop))) {
#		perm <- c( rownames(drop)[rowData(drop)$Permuted], rownames(facs)[rowData(facs)$Permuted] )
#		drop <- drop[!rownames(drop) %in% perm, ]
#		facs <- facs[!rownames(facs) %in% perm, ]
#	}

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

	xDataOUT <- vector() #matrix(nrow=length(methods), ncol=5)
	#rownames(xDataOUT) <- methods_nice_names
	#colnames(xDataOUT) <- c("Ndrop", "Nfacs", "Nboth", "Nconsist", "perc_consist")

	Method_Markers_Drop <- list();
	Method_Markers_Facs <- list();
	
	tissue <- rownames(file_pairs)[i]
	#tissue_outfile <- paste(tissue, "_", auc_threshold, "_consistency_results.rds", sep="");
	tissue_outfile <- paste(tissue, "_consistency_results.rds", sep="");

	# Identify Markers
	for (m in 1:length(methods)) {
		meth <- methods[m];
		# Get markers
		mark_drop <- complex_markers(as.matrix(assays(drop)[[meth]]), drop$cell_type1, n_max=1) 
		mark_facs <- complex_markers(as.matrix(assays(facs)[[meth]]), factor(facs$cell_type1), n_max=1) 
		# p-value correlation
		tmp_facs <- mark_facs[rownames(mark_facs) %in% rownames(mark_drop),]
		tmp_drop <- mark_drop[match(rownames(tmp_facs), rownames(mark_drop)),]
		c <- cor(tmp_facs$p.value, tmp_drop$p.value, method="spearman")
		p_val_cors[i,m] <- c
		colnames(p_val_cors) <- remap_names(methods)

		dataout <- vector();
		for (auc_threshold in c(0, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)){
			# apply thresholds
			mark_drop_gd <- mark_drop[mark_drop$AUC > auc_threshold & mark_drop$q.value < 0.05,]
			mark_facs_gd <- mark_facs[mark_facs$AUC > auc_threshold & mark_facs$q.value < 0.05,]

			if (auc_threshold == 0) {
				Method_Markers_Drop[[meth]] <- cbind(
				rownames(mark_drop_gd), mark_drop_gd$Group)
				Method_Markers_Facs[[meth]] <- cbind(
				rownames(mark_facs_gd), mark_facs_gd$Group)
			}
			Ndrop <- nrow(mark_drop_gd)
			Nfacs <- nrow(mark_facs_gd)

			# Consistent
			mark_drop_gd <- mark_drop_gd[rownames(mark_drop_gd) %in% rownames(mark_facs_gd), ]
			mark_facs_gd <- mark_facs_gd[match(rownames(mark_drop_gd), rownames(mark_facs_gd)), ]
			Nboth <- nrow(mark_drop_gd);

			mark_drop_gd$Group <- as.character(mark_drop_gd$Group)
			mark_facs_gd$Group <- as.character(mark_facs_gd$Group)
			perc_agreement <- sum(mark_drop_gd$Group == mark_facs_gd$Group)/nrow(mark_drop_gd);
			dataout <- rbind(dataout, c(auc_threshold, Ndrop, Nfacs, Nboth, sum(mark_drop_gd$Group == mark_facs_gd$Group), perc_agreement))
		}
		dataout <- data.frame(tissue=as.character(tissue),method=as.character(meth), dataout, stringsAsFactors=FALSE)
		colnames(dataout) <- c("tissue","method", "auc", "ndrop", "nfacs", "nboth", "nagree", "percagree")
		xDataOUT <- rbind(xDataOUT, dataout);
	}

	# Cross method consistency
	DropxMethOUT <- matrix(nrow=length(methods), ncol=length(methods))
	rownames(DropxMethOUT) <- methods_nice_names
	colnames(DropxMethOUT) <- methods_nice_names

	FACSxMethOUT <- matrix(nrow=length(methods), ncol=length(methods))
	rownames(FACSxMethOUT) <- methods_nice_names
	colnames(FACSxMethOUT) <- methods_nice_names

	for (r in 1:length(methods)) {
	for (j in 1:length(methods)) {
		meth1 <- methods[r];
		meth2 <- methods[j];
		markers1 <- Method_Markers_Drop[[meth1]]
		markers2 <- Method_Markers_Drop[[meth2]]
		# Markers by both methods
		markers1 <- markers1[markers1[,1] %in% markers2[,1],]
		markers2 <- markers2[match(markers1[,1], markers2[,1]),]
		# % same group
		perc_agree1 <- sum( markers1[,2] == markers2[,2])/nrow(markers1);
		DropxMethOUT[r,j] <- perc_agree1

		markers1 <- Method_Markers_Facs[[meth1]]
		markers2 <- Method_Markers_Facs[[meth2]]
		# Markers by both methods
		markers1 <- markers1[markers1[,1] %in% markers2[,1],]
		markers2 <- markers2[match(markers1[,1], markers2[,1]),]
		# % same group
		perc_agree2 <- sum( markers1[,2] == markers2[,2])/nrow(markers1);
		FACSxMethOUT[r,j] <- perc_agree2
	}}

	saveRDS(list(crossData=xDataOUT, DropxMethOUT=DropxMethOUT, FACSxMethOUT=FACSxMethOUT), file=tissue_outfile )
	allDataOUT <- rbind(allDataOUT, xDataOUT)
	if (is.null(dim(allMethOUT))) {
		allMethOUT <- DropxMethOUT+FACSxMethOUT
	} else {
		allMethOUT <- allMethOUT+DropxMethOUT+FACSxMethOUT
	}
		
}
allMethOUT <- allMethOUT/(nrow(file_pairs)*2)
saveRDS(allMethOUT, file=cross_method_outfile) 
saveRDS(allDataOUT, file=cross_dataset_outfile) 
saveRDS(p_val_cors, file="TM_mark_cors.rds") 

