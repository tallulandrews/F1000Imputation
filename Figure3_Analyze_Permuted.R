#Set-up
require("SingleCellExperiment")
require("scater")
require("Matrix")
source("~/MAGIC/CTP_Functions.R")

#args <- commandArgs(trailingOnly=TRUE)
#file <- args[1]


FILES <- Sys.glob("*_permuted.rds")
ALL_FPR <- list();
for (file in FILES) {

prefix <- unlist(strsplit(file, "[./]"))
prefix <- prefix[length(prefix)-1]

set.seed(4817)

# Summary Statistics
require("scater")

obj <- readRDS(file)
check_FPs <- function(mat, non_DE, norm=FALSE, out.type=c("FPR", "FPs")) {
	if (norm) {
		sf <- Matrix::colSums(mat)
		mat <- t( t(mat)/sf*median(sf))
	}
	mat <- mat[non_DE, ]
	require("Hmisc")
	type1 <- obj@metadata$type_pair[1]
	type2 <- obj@metadata$type_pair[2]
	DE_p.value <- apply(mat, 1, function(x) {
				wilcox.test(x[obj$cell_type1==type1], 
					    x[obj$cell_type1==type2])$p.value})
	DE_p.value[is.na(DE_p.value)] <- 1;

	threshold <- 0.05/length(non_DE)
	FPs <- sum( DE_p.value < threshold, na.rm=T);
	FPR <- FPs/nrow(mat)
	if (out.type[1] == "FPR") {
		return(FPR)
	} else if (out.type =="FPs") {
		return(rownames(mat)[DE_p.value < threshold])
	}
}

FPR <- list()
non_DE <- rowData(obj)$Permuted
methods <- names(assays(obj))
if (length(grep("_perm_", methods))) {
	methods <- methods[-grep("_perm_", methods)]
}
for(m in methods) {
	if (m == "counts") {
		out <- check_FPs(assays(obj)[[m]], non_DE, norm=TRUE)
	} else {
		out <- check_FPs(assays(obj)[[m]], non_DE, norm=FALSE)
	}
	FPR[[m]] <- out
}

#saveRDS(FPR, paste(prefix,"FPR.rds", sep="_"))
ALL_FPR[[file]] <- FPR;

}
saveRDS(ALL_FPR, file="Permuted_FPR.rds")

### Check distribution for FP vs non-FP ###
check.NB.fit <- function(x, lib.size=rep(1, length(x)), max_r=10^10) {
	# Fit
	lib.size <- lib.size/mean(lib.size)
	mus <- mean(x)*lib.size
	obs_err <- sum( (x - mus)^2 )
	rg <- sum( mus^2 )/(obs_err - sum(mus))
	if (rg <= 0) {rg <- max_r}
	# Quality of fit
	ps <- sapply(1:length(x), function(i) {pnbinom(x[i], mu=mus[i], size=rg)})
	ll <- log10(prod(ps))
	return(ll);
}

check.ZINB.fit <- function(x, lib.size=rep(1, length(x)), max_r=10^10, e=0.00001) {
	# Fit
	lib.size <- lib.size/mean(lib.size)
	d.obs <- mean(x==0)
	if (d.obs == 0) {return(check.NB.fit(x, lib.size))}
	d_curr <- d.obs;
        d_prev <- -100;
	nc <- length(x);
	while( abs(d_curr-d_prev) > e ) {
		mus <- sum(x)/(nc-d_curr*nc)*lib.size
		weights <- rep(1, times=length(x))
		weights[x == 0] <- (1-d_curr/d.obs)
		obs_err <- sum( (x - mus)^2*weights )
		rg <- sum( mus^2*weights )/(obs_err - sum(mus*weights))
                if (rg <= 0) {rg <- max_r}
		pds <- (1 + mus/rg)^(-rg)
                d.exp <- mean(pds)
		d_prev <- d_curr
                d_curr <- (d.obs - d.exp)
		if (d_curr <= 0) {d_curr <- d_prev}
	}
	# params : mus, d_prev, rg
	p0s <- d_prev+sapply(mus, function(m) {pnbinom(0, mu=m, size=rg)})
	ps <- sapply(1:length(x), function(i) {pnbinom(x[i], mu=mus[i], size=rg)})
	ps[x==0] <- p0s[x==0];
	ll <- log10(prod(ps))
	return(ll);
}

check.ZILN.fit <- function(x, lib.size=rep(1, length(x))) {
	lib.size <- lib.size/mean(lib.size)
	ln.obs <- log2(x/lib.size+1)
	dr <- mean(ln.obs==0); # dropout rate/zero inflation
	mu <- mean(ln.obs[ln.obs>0]);
	sigma <- sd(ln.obs[ln.obs>0]);
	p0s <- dr;
	ps <- pnorm(ln.obs, mean=mu, sd=sigma)
	ps[ln.obs==0] <- p0s;
	ll <- log10(prod(ps))
        return(ll);
}

FILES <- Sys.glob("*_permuted.rds")
OUT <- list();
OUT[["all"]] <- vector()

for (file in FILES) {

# Set up
obj <- readRDS(file)
obj <- obj[,obj$cell_type1 %in% obj@metadata$type_pair];
obj <- obj[Matrix::rowSums(obj@assays[["counts"]] > 0) > 5,]
lib.size <- Matrix::colSums(obj@assays[["counts"]])
non_DE <- rowData(obj)$Permuted
methods <- names(assays(obj))
if (length(grep("_perm_", methods))) {
	methods <- methods[-grep("_perm_", methods)]
}
FPs <- list(); 
for(m in methods) {
        if (m == "counts") {
                out <- check_FPs(assays(obj)[[m]], non_DE, norm=TRUE, out.type="FPs")
        } else {
                out <- check_FPs(assays(obj)[[m]], non_DE, norm=FALSE, out.type="FPs")
        }
        FPs[[m]] <- out
}

best_fits <- apply(obj@assays[["counts"]][non_DE,], 1, function(x) { c(check.NB.fit(x, lib.size), check.ZINB.fit(x,lib.size), check.ZILN.fit(x,lib.size))})
best_fits <- t(best_fits)
colnames(best_fits) <- c("NB", "ZINB", "ZILN")

best_model <- apply(best_fits, 1, function(a){out <- colnames(best_fits)[which(a==max(a))]; if (length(out) > 1) {return("None")} else {return(out)}})

best_model <- factor(unlist(best_model), levels=c("None", "NB", "ZINB", "ZILN"))

# Save results
rnames <- c(rownames(OUT[["all"]]), file)
OUT[["all"]] <- rbind(OUT[["all"]], table(best_model))
rownames(OUT[["all"]]) <- rnames;
for(m in methods) {
if (is.null(OUT[[m]])){
	OUT[[m]] <- vector()
}
rnames <- c(rownames(OUT[[m]]), file)
OUT[[m]] <- rbind(OUT[[m]], table(best_model[names(best_model) %in% FPs[[m]]]) )
rownames(OUT[[m]]) <- rnames;
}

}
saveRDS(OUT,  file="Permuted_Fits.rds")

