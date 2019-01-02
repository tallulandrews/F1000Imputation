source("~/MAGIC/R_Imputation_Functions.R")
require("scater")
require("methods")
# wd : /lustre/scratch117/cellgen/team218/TA/Simulations_Temporary_Files/Imputation/vsBulk/

# Create SCE - Kolo
dat <- read.table("counttable_es.csv", header=T)
dat <- dat[,-grep("_2i_", colnames(dat))]
dat <- dat[-grep("^ERCC-", rownames(dat)), ]
dat <- dat[-grep("__", rownames(dat)), ]
dat <- dat[rowSums(dat > 0) > 5,]
ann <- strsplit(colnames(dat), "_")
batch <- unlist(lapply(ann, function(a){a[4]}))
Group <- unlist(lapply(ann, function(a){a[3]}))
ann <- data.frame(batch=batch, Group=Group)
rownames(ann) <- colnames(dat)
sce_kolo <- SingleCellExperiment(assays=list(counts=as.matrix(dat)), colData=ann)
sce_kolo <- normalize(sce_kolo);

# Impute
sce_kolo <- Impute_default_all(sce_kolo)
# Save
saveRDS(sce_kolo, "Kolo_imputed.rds");

# Create SCE - Tung
dat <- read.table("GSE77288_molecules-raw-single-per-sample.txt.gz", header=T)
ann <- dat[,1:3]
dat <- dat[,-c(1:3)]
dat <- t(dat)
colnames(dat) <- paste("Cell", 1:ncol(dat));
dat <- dat[-grep("^ERCC_", rownames(dat)), ]
dat <- dat[rowSums(dat > 0) > 5,]
rownames(ann) <- colnames(dat);
ann$Group <- ann$individual
sce_tung <- SingleCellExperiment(assays=list(counts=dat), colData=ann)
sce_tung <- normalize(sce_tung);
# Impute
sce_tung <- Impute_default_all(sce_tung)
# Save
saveRDS(sce_tung, "Tung_imputed.rds");
