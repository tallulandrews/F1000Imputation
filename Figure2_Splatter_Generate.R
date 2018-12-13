source("/nfs/users/nfs_t/ta6/MAGIC/Splatter_Functions.R")
require("Rmagic") # t_diffusion
require("DrImpute") # dropout.probability.threshold
require("scImpute") # drop_thre
require("SAVER") # percent of genes to predict
require("scater")
source("/nfs/users/nfs_t/ta6/MAGIC/knn_smooth.R") # k
require(doParallel)

# Get table of simulation parameters
param_table <- Generate_param_table()

n_cores=16
imputation_methods=c("MAGIC", "scImpute", "DrImpute", "SAVER", "knn")

# Generate simulated datasets:
for (i in 1:nrow(param_table)) {
	filename <- paste("vsSimParams_paramset_", i, ".rds", sep=""); # objects+imputation
	if (file.exists(filename)) {next;}
	sim <- Sim_w_splatter(ngenes=param_table[i,"ngenes"], ncells=param_table[i,"ncells"], 
				ngroups=param_table[i,"ngroups"], seed=param_table[i,"seed"],
				dropouts=param_table[i,"dropouts"], propDE=param_table[i, "propDE"],
				method=as.character(param_table[i, "method"])
				)
	for (meth in imputation_methods) {
		registerDoParallel(cores = n_cores)
		imputed <- do_imputation_splatter(sim, method=meth, n.cores=n_cores)
		assays(sim)[[meth]] <- imputed
	}
	saveRDS(sim, filename)
}

