# Duplicated from https://github.com/tallulandrews/CellTypeProfiles for convenience.
factor_counts <- function(vec) {
        vec <- factor(vec)
        x <- split(seq(length(vec)), vec)
        result <- sapply(x, function(a) length(a))
        return(result);
}

my_row_mean_aggregate <- function(mat, groups) {
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) rowMeans(MAT[,a]))

        return(result);
}

my_row_var_aggregate <- function(mat, groups) {
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) matrixStats::rowVars(MAT[,a]))

        return(result);
}

ctp_fast_AUC <- function(expression_vec, truth) {
        # using Mann-Whitney U test

        R = rank(expression_vec);
        N1 = sum(truth)
        N2 = sum(!truth);
        #U1 = sum(R[truth])-N1*(N1+1)/2
        U2 = sum(R[!truth])-N2*(N2+1)/2
        if (N1 == 0) {return(c(0,0,0))}
        if (N2 == 0) {return(c(1,1,1))}
        AUC = 1-U2/(N1*N2); # smaller valued ranks (i.e. lower expression values for !true).
        # assumes large sample size
        #  https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_the_Area_Under_an_ROC_Curve.pdf
        # originally (Hanley and McNeil 1982)
        Q1 = AUC/(2-AUC)
        Q2 = 2*AUC^2/(1+AUC)
        SE = sqrt((AUC*(1-AUC)+ (N1-1)*(Q1-AUC^2)+(N2-1)*(Q2-AUC^2))/(N1*N2))
        return(c(max(0, AUC-1.96*SE),AUC, min(1, AUC+1.96*SE)));
}



complex_markers <- function (expr_mat, labels, n_max=length(unique(labels))-1, strict_only = FALSE) {
        if (length(labels) != length(expr_mat[1, ])) {
                stop("Length of labels does not match number of cells.");
        }
        if (n_max > (length(unique(labels))-1)) {
                stop("n_max must be less than the number of labels");
        }
	# Remove groups with a single cell
	label_counts = factor_counts(labels)
	exclude = names(label_counts)[label_counts<2]
	if (length(exclude) > 0) {
		warning(paste("Warning: Excluding",length(exclude),"groups with less than 2 samples."))
#		print(paste("Warning: Excluding",length(exclude),"groups with less than 2 samples."))
		keep = !(labels %in% exclude)
		expr_mat <- expr_mat[,keep]
		labels <- labels[keep]
		labels <- factor(labels);
	}
	if (min(label_counts) < 10) {
		print("Warning: Small groups (n < 10) may bias marker gene results.")
	}

        # Mean ranked expression all genes, each cluster (efficient)
        gene_cluster_means <- function(mat, groups) {
                MAT <- as.matrix(mat)
                x <- split(seq(ncol(MAT)), groups)
                result <- sapply(x, function(a) rowMeans(MAT[,a]))
                return(result);
        }
        ranked_matrix <- t(apply(expr_mat, 1, rank))
        gene_cluster_ranks <- gene_cluster_means(ranked_matrix, labels)
        cluster_priority <- t(apply(-1*gene_cluster_ranks, 1, rank))
	label_columns <- sort(unique(labels));
        # Consider each gene a possible marker for top 1->n_max groups
        gene_auc <- function(g) {
		# Short circuit for genes which invariant
		if (sum(cluster_priority[g,] == cluster_priority[g,1]) == length(cluster_priority[g,])) {
			group=rep(0, times=length(label_columns));
                        pval=-1
                        auc=-1
			if (n_max == 1) {
				return(c(auc, "NA", pval))
			} else {
				return(c(auc, group, pval))
			}
		}


		# Get AUC + 95% CI for this gene in top n groups
                get_auc_ci <- function(n) {
                        g_groups <- colnames(cluster_priority)[which(cluster_priority[g,] <= n)];

                        if (length(g_groups) < length(unique(labels)) & length(g_groups) > 0) {
				return(ctp_fast_AUC(expr_mat[g,], labels %in% g_groups));
                        } else {
                                return(c(-1,-1,-1))
                        }
                }
		# AUCs for this gene at all ns
                auc_tab <- sapply(1:n_max, get_auc_ci)

		# short circuit if no valid aucs
		if (max(auc_tab[2,]) < 0 ) {
			group=rep(0, times=length(label_columns));
			if (n_max == 1) { group="NA" }
                        pval=-1
                        auc=-1
			return(c(auc, group, pval))
		} 

		# Identify top set of groups this gene is a marker for & determine if good enough to return.
                top = which(auc_tab[2,] == max(auc_tab[2,]))
		
		if (n_max == 1) {
			sec <- top
		} else {
                	sec = which(auc_tab[2,] == max(auc_tab[2,-top]))
		}
                if (n_max > 1 & min(auc_tab[1,top]) <= max(auc_tab[2,sec]) & strict_only) {
                        # Not a marker : ci of top set of groups contains auc of second best set of groups
                        group=rep(0, times=length(label_columns));
                        pval=-1
                        auc=-1
                } else {
                        # Return marker info
                        n = max(top)
                        g_groups = colnames(cluster_priority)[which(cluster_priority[g,] <= n)];
			if (length(g_groups) < 1 & n_max == 1) {
				group="None/Tied"
				pval=-1
				auc=-1
			} else {
			if (n_max == 1) {
	                        group = paste(g_groups, collapse="+")
			} else {
				group = as.numeric(label_columns %in% g_groups);
			}
                        auc = auc_tab[2,n]
			# p.value from wilcox test
                        pval = wilcox.test(expr_mat[g,!(labels %in% g_groups)],expr_mat[g,(labels %in% g_groups)])$p.value
			}
                }
                return(c(auc, group, pval))
        }
	# Get best AUC across sets of groups for each gene
        out <- sapply(1:length(expr_mat[,1]),gene_auc)

	# Format output nicely
        out_matrix = as.data.frame(t(out))
        rownames(out_matrix) <- rownames(expr_mat)
	if (n_max == 1) {
	        colnames(out_matrix) <- c("AUC", "Group",  "p.value")
	} else {
		colnames(out_matrix) <- c("AUC", as.character(label_columns), "p.values")
	}
        out_matrix[,1] = as.numeric(as.character(out_matrix[,1]))
        out_matrix[,3] = as.numeric(as.character(out_matrix[,3]))
	# Apply Bonferroni correction
	n_possible_group_combos <- choose(length(unique(labels)), n_max)
        out_matrix$q.value = out_matrix$p.value*n_possible_group_combos*length(expr_mat[,1]);
        out_matrix$q.value[out_matrix$q.value < 0] = -1;
        out_matrix$q.value[out_matrix$q.value > 1] = 1;
        return(out_matrix);
}

# Turn matrix of on/off into combined names of positive or negative groups
get_combo_names <- function(marker_matrix) {
	tmp <- marker_matrix[,2:(length(marker_matrix[1,])-2)]
	out <- apply(tmp, 1, function(x){paste(colnames(tmp)[x==1], collapse="+")})
	# If gene is "off" in fewer groups than it is "on" it is a negative marker
	anti = rowSums(tmp) > 0.5*length(tmp[1,])
	out_anti <- apply(tmp, 1, function(x){paste(colnames(tmp)[x==0], collapse="+")})
	out[anti] = paste("NOT", out_anti[anti])
	return(out);
}

