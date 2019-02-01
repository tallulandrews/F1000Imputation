source("/nfs/users/nfs_t/ta6/MAGIC/Splatter_Functions.R")
source("~/MAGIC/R_Imputation_Functions.R")
n_cores=16

# Get table of simulation parameters
param_table <- Generate_param_table()

# Generate simulated datasets:
for (i in 1:nrow(param_table)) {
	filename <- paste("vsSimParams_wBatch_paramset_", i, ".rds", sep=""); # objects+imputation
	if (file.exists(filename)) {next;}
	sim <- Sim_w_splatter(ngenes=param_table[i,"ngenes"], ncells=param_table[i,"ncells"], 
				ngroups=param_table[i,"ngroups"], seed=param_table[i,"seed"],
				dropouts=param_table[i,"dropouts"], propDE=param_table[i, "propDE"],
				method=as.character(param_table[i, "method"]), nbatchs=3
				)
	imputed <- Impute_default_all(sim, n.cores=n_cores)
	saveRDS(imputed, filename)
}
