makePCAPlotByPhenoType <- function(se, title="Title", figDir, threshold.count=5) {
    require(DESeq2)
    require(ggplot2)
    dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > threshold.count, ],
                        design = ~ pheno_type)
    rlg <- rlog(dds)
    data <- plotPCA(rlg, intgroup=c("pheno_type"), returnData=TRUE)
    percentVar <- round(100 * attr(data, "percentVar"))
    ggplot(data, aes(x=PC1, y=PC2, color=group, label=rownames(data))) +
               geom_point() + 
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                        labs(list(title=title)) +
                            geom_text(hjust=0, vjust=0.5)

    ggsave(file=file.path(figDir, paste0(title, "_SampleDistancePCAplot.png")))
}

