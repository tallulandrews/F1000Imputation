source("/nfs/users/nfs_t/ta6/MAGIC/Basic_NB_Simulations.R")

set.seed(20194)
nrep=10;
my_seeds <- round(runif(nrep*2)*10000)
my_seeds <- unique(my_seeds)
my_seeds <- my_seeds[1:nrep]

for (s in my_seeds) {
	set.seed(s)
	sims <- Sim(n_cells=2000, max_mean=4, min_mean=-3, 
			propDE=0.5, type="clust")
	sims$counts <- sims$mat;
	require("scater")
	sce <- SingleCellExperiment(assays=list(counts=sims$counts), colData=data.frame(cell_truth=sims$cell_truth), rowData=data.frame(g_up=sims$g_up, g_down=sims$g_down))
	saveRDS(sce, file=paste("vsMethodParams_seed_dca_",s,".rds", sep=""));
	saveRDS(sce, file=paste("vsMethodParams_seed_scVI_",s,".rds", sep=""));
	#write.table(sims$counts, row.names=TRUE, col.names=TRUE, file=paste("vsMethodParams_seed_", s, ".csv", sep=""), sep=",")
}
