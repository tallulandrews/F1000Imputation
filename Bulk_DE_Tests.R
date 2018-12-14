do_edgeR <- function(data, labels, fdr_thresh=2) {
        require("edgeR")
        mygroups = as.character(unique(labels));
        edgerdata = DGEList(counts=data,group=labels)
        edgerdata = calcNormFactors(edgerdata)
        edgerdata = estimateCommonDisp(edgerdata)
        edgerdata = estimateTagwiseDisp(edgerdata)
        normalized = cpm(edgerdata)
	TABLE <- c("JUNK",NA,NA,NA,NA,NA);
        for (i in 1:(length(mygroups)-1)) {
                group1 = mygroups[i];
                for (j in (i+1):length(mygroups)) {
                        group2 = mygroups[j];
                        print (paste(group1,"vs",group2));
                        res = exactTest(edgerdata, pair=c(group1,group2));
			out=cbind(res$table$PValue,p.adjust(res$table$PValue, method="BH"));
			rownames(out) = rownames(res$table);
			out = out[!is.na(out[,1]),]
	
               	        output =  data.frame(Gene=rownames(out),
				    group1=rowMeans(data[match(rownames(out), rownames(data)),labels==group1]),
            	                    group2=rowMeans(data[match(rownames(out), rownames(data)),labels==group2]),
				    pval=out[,1],qval=out[,2], 
				    test=rep(paste(group1,"_vs_",group2,sep=""),times=length(out[,1])))
			TABLE = rbind(TABLE,output)
                }
        }
        colnames(TABLE) <- c("Gene", "Mean1", "Mean2", "p.value","q.value","test")
	TABLE$q.value = p.adjust(TABLE$p.value, method="fdr")
        return(TABLE[2:length(TABLE[,1]),])
}

do_DESeq2 <- function(data, labels, fdr_thresh=2, ncores=1) {
	colnames(data) = 1:length(data[1,]) # prevent duplicate rownames when this is transposed by DESeq
        require("DESeq2")
        if (length(data[1,]) > 1000) {
                set.seed(101);
                subset = subsample_forDE(labels, 1000);
                data = data[,subset]
                labels = labels[subset]
        }

        mygroups = as.character(unique(labels));
        labels = data.frame(labels); rownames(labels) = colnames(data);
	TABLE <- c("JUNK",NA,NA,NA,NA,NA);
        for (i in 1:(length(mygroups)-1)) {
                group1 = mygroups[i];
                for (j in (i+1):length(mygroups)) {
                        group2 = mygroups[j];
                        deseqdata = DESeqDataSetFromMatrix(countData = round(data[,labels==group1 | labels==group2]), colData = data.frame(group = labels[labels==group1 | labels==group2]), design = ~ group);
			res=results(DESeq(deseqdata), parallel=n.cores)
		
			res$padj[is.na(res$padj)] = 1;
	
               	        output =  data.frame(Gene=rownames(res),
				    group1=rowMeans(data[match(rownames(res), rownames(data)),labels==group1]),
            	                    group2=rowMeans(data[match(rownames(res), rownames(data)),labels==group2]),
				    pval=res$pvalue,qval=res$padj, 
				    test=rep(paste(group1,"_vs_",group2,sep=""),times=length(res$padj)))
			TABLE = rbind(TABLE,output)
                }
        }
	colnames(TABLE) <- c("Gene", "Mean1", "Mean2", "p.value","q.value","test")
	TABLE$q.value = p.adjust(TABLE$p.value, method="fdr")
        return(TABLE[2:length(TABLE[,1]),])
}

do_DESeq <- function(data, labels, fdr_thresh=2) {
	colnames(data) = 1:length(data[1,]) # prevent duplicate rownames when this is transposed by DESeq
        require("DESeq")
        mygroups = unique(as.character(labels));
        deseqdata = newCountDataSet(round(data), as.character(labels));
        deseqdata = estimateSizeFactors(deseqdata);
        deseqdata = estimateDispersions(deseqdata); #gene-est-only can only be used if many replicates.
	TABLE <- c("JUNK",NA,NA,NA,NA,NA);
        for (i in 1:(length(mygroups)-1)) {
                group1 = mygroups[i];
                for (j in (i+1):length(mygroups)) {
                        group2 = mygroups[j];
                        print (paste(group1,"vs",group2));
                        res = nbinomTest(deseqdata, group1, group2)
                        output =  data.frame(Gene=res$id, group1=res$baseMeanA, group2=res$baseMeanB,pval=res$pval,qval=res$padj, test=rep(paste(group1,"_vs_",group2,sep=""),times=length(res$padj)))
			TABLE = rbind(TABLE,output)
                }
        }
	colnames(TABLE) <- c("Gene", "Mean1", "Mean2", "p.value","q.value","test")
	TABLE$q.value = p.adjust(TABLE$p.value, method="fdr")
	return(TABLE[2:length(TABLE[,1]),])
}

do_LimmaVoom <- function(data, labels, tag="Dataset", fdr_thresh=2, suppress.plot=TRUE) {
	require("limma")
	require("edgeR")
	# make design matrix
	ncol = length(unique(labels));
	labels = factor(labels);
	design = model.matrix(~labels, data=labels)

	dge <- DGEList(counts=data);
	dge <- calcNormFactors(dge, method="TMM");
	voomed <- voom(dge, design, plot=!suppress.plot)
	fit <- lmFit(voomed, design)
	fit <- eBayes(fit)

	TABLE <- topTable(fit, number=length(data[,1]), p.value=fdr_thresh)
	last = length(TABLE[1,])
	if (ncol > 2) {
		TABLE <- cbind(rownames(TABLE), rep(1, times=length(TABLE[,1])), TABLE[,c(1:(length(unique(labels))-1))],TABLE[,c(last-1,last)])
		colnames(TABLE) <- c("Gene", levels(labels),"p.value","q.value")
	} else {
		TABLE <- cbind(rownames(TABLE), TABLE[,1], TABLE[,2], TABLE[,4], TABLE[,5])
		colnames(TABLE) <- c("Gene", "logFC", "AvgExpr", "p.value", "q.value")
	}
        return(TABLE)
}

do_wilcox <- function(data,labels,fdr_thresh=2) {
        mygroups = as.character(unique(labels));
	norm <- t(t(data)/colSums(data)*1000000)

        TABLE <- c("JUNK",NA,NA,NA,NA,NA);
        for (i in 1:(length(mygroups)-1)) {
                group1 = mygroups[i];
                for (j in (i+1):length(mygroups)) {
                        group2 = mygroups[j];
			
			pvals <- unlist(apply(norm, 1, function(x) {wilcox.test(x[labels==group1], x[labels==group2])$p.value}))
			pvals[is.na(pvals)] = 1;
			

               	        output =  data.frame(Gene=rownames(norm),
				    group1=rowMeans(norm[,labels==group1]),
            	                    group2=rowMeans(norm[,labels==group2]),
				    pval=pvals,qval=p.adjust(pvals, method="fdr"), 
				    test=rep(paste(group1,"_vs_",group2,sep=""),times=length(pvals)))
			TABLE = rbind(TABLE,output)
                }
        }
        colnames(TABLE) <- c("Gene", "Mean1", "Mean2", "p.value","q.value","test")
	TABLE$q.value = p.adjust(TABLE$p.value, method="fdr")
        return(TABLE[2:length(TABLE[,1]),])
}
