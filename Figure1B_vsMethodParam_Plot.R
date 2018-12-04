require("matrixStats")
# Better Figures
methods <- c("scImpute", "DrImpute", "SAVER", "MAGIC", "MAGIC_k", "knn");
methods_titles <- c("scImpute", "DrImpute", "SAVER", "MAGIC", "MAGIC", "KNN");
#default_vals <- c(0.5, 0, 1, 3, NA);
png("vsMethodParam_Figure1.png", width=4*1.8, height=4, units="in", res=300)
par(mfrow=c(2,3))
for(i in 1:length(methods)) {
	m <- methods[i]
	file <- paste(m, "vs_param.rds", sep="_")
	out <- readRDS(file)
	FPR <- out$FP/(out$FP+out$TN)
	par(mar=c(3.5,3,1.5,1))
	param_range <- out$params
	names(param_range) <- param_range;
	names(param_range)[1] <-"Raw"
	plot(1:length(param_range), rowMeans(FPR), main=methods_titles[i], xlab="", 
	ylab="FPR", type="b", lwd=2, col="black", ylim=c(0,1), xaxt="n", las=3)
	if (i == 2) {
	title(xlab="Zeros Remaining", line=2.2)
	} else {
	title(xlab=out$param_name, line=2.2)
	}
	title(ylab="FPR", line=2)
	axis(1, at=1:length(param_range), labels=names(param_range))
	row_sd <- sqrt(rowVars(FPR))
	lines(1:length(param_range), rowMeans(FPR)+row_sd*2, lty=2)
	lines(1:length(param_range), rowMeans(FPR)-row_sd*2, lty=2)

}
dev.off()


png("vsMethodParam_Figure1Suppl1.png", width=4*1.8, height=4, units="in", res=300)
par(mfrow=c(2,3))
for(i in 1:length(methods)) {
        m <- methods[i]
        file <- paste(m, "vs_param.rds", sep="_")
        out <- readRDS(file)
        sens <- out$TP/(out$TP+out$FN)
        spec <- out$TN/(out$TN+out$FP)
        par(mar=c(3.5,3,1.5,1))
        param_range <- out$params
        names(param_range) <- param_range;
        names(param_range)[1] <-"Raw"
        plot(1:length(param_range), rowMeans(sens), main=methods_titles[i], xlab="", 
        ylab="", type="b", lwd=2, ylim=c(0,1), xaxt="n", las=3, col="black")
        if (i == 2) {
        title(xlab="Zeros Remaining", line=2.2)
        } else {
        title(xlab=out$param_name, line=2.2)
        }
        title(ylab="TPR", line=2)
        axis(1, at=1:length(param_range), labels=names(param_range))
        row_sd <- sqrt(rowVars(sens)/ncol(sens))
        lines(1:length(param_range), rowMeans(sens)+row_sd*2, lty=2, col="black")
        lines(1:length(param_range), rowMeans(sens)-row_sd*2, lty=2, col="black")
}
dev.off()

