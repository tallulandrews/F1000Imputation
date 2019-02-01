require("matrixStats")
source("~/MAGIC/Colour_Scheme.R")

files <- Sys.glob("*_vs_param.rds")
methods <- sub("_vs_param.rds", "", files)
m_order <- Method_Colours[ match(methods, Method_Colours[,1]),2]
files <- files[!is.na(m_order)]
methods <- methods[!is.na(m_order)]
m_order <- m_order[!is.na(m_order)]
m_names <- names(Method_Legend)[match(m_order, Method_Legend)]

m_order <- order(as.numeric(factor(m_names, levels=names(Method_Legend))))
files <- files[m_order]
methods <- methods[m_order]
m_names <- m_names[m_order]

FPR_col <- "black"
TPR_col <- "dodgerblue"


png("Figure1_SensSpec_Lines.png", width=8, height=4, units="in", res=300)
par(mfrow=c(2,4))
for (i in 1:length(files)) {
	f <- files[i]
	m <- strsplit(f, "_")[[1]]
	res <- readRDS(f);
	FPR <- res$FP/(res$FP+res$TN)
        TPR <- res$TP/(res$TP+res$FN)
	
	par(mar=c(3.5,3,1.5,1))
	# xlab
	param_range <- res$params
	names(param_range) <- param_range;
	names(param_range)[1] <-"Raw"
	# FPR
	plot(1:length(param_range), rowMeans(FPR), main=m_names[i], xlab="", 
	ylab="Score", type="b", lwd=2, col=FPR_col, ylim=c(0,1), xaxt="n", las=3)
	title(xlab=res$param_name, line=2.2)
	axis(1, at=1:length(param_range), labels=names(param_range))
	row_sd <- sqrt(rowVars(FPR))
	lines(1:length(param_range), rowMeans(FPR)+row_sd*2, lty=2, col=FPR_col)
	lines(1:length(param_range), rowMeans(FPR)-row_sd*2, lty=2, col=FPR_col)

	# TPR
	par(new=TRUE)
	plot(1:length(param_range), rowMeans(TPR), main=m_names[i], xlab="", 
	ylab="", type="b", lwd=2, col=TPR_col, ylim=c(0,1), xaxt="n", las=3)
	row_sd <- sqrt(rowVars(TPR))
	lines(1:length(param_range), rowMeans(TPR)+row_sd*2, lty=2, col=TPR_col)
	lines(1:length(param_range), rowMeans(TPR)-row_sd*2, lty=2, col=TPR_col)
	par(new=FALSE)

}
blank_plot <- function() {
        tmp <-  par("mar")
        par(mar=c(0,0,0,0))
        plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
        return(tmp);
}
blank_plot()
legend("center", bty="n", c("TPR", "FPR", "Mean", "95% CI"), lty=c(NA,NA,1,2), pch=c(1,1,NA,NA), col=c(TPR_col, FPR_col, "black", "black"), cex=1.25, title="Quality Scores", lwd=c(2,2,2,1))
dev.off()



