simpleDESeqReport <- function(dds, OrgDb, lfcThreshold=0,
                              gene_id_type=c("ENTREZID", "ENSEMBL"),
                              thres.padj=NULL, filename) {
     ### dds is the returned values from DESeq
    gene_id_type <- match.arg(gene_id_type)
    require(DESeq2)

    df <- results(dds, lfcThreshold=lfcThreshold)
    if (!is.null(thres.padj))
        df <- subset(df, df$padj < thres.padj)

    df <- df[order(df$padj, decreasing=FALSE), ]
    ## annotation (EntrzID, ENSEMBL and GENENAME)
    if (gene_id_type == "ENTREZID") {
        anno <- select(x=OrgDb, keys=rownames(df),
                       keytype="ENTREZID", columns=c("ENSEMBL", "SYMBOL", "GENENAME"),
                       multiVals="first")
        anno <- anno[!duplicated(anno$ENTREZID), ]
        rownames(anno) <- anno$ENTREZID
        anno <- anno[rownames(df), ]
    }

    if (gene_id_type == "ENSEMBL") {
        keys <- sapply(strsplit(rownames(df), ".", fix=TRUE), "[[", 1)
        anno <- select(x=OrgDb, keys=keys,
                       keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"),
                       multiVals="first")
        anno <- anno[!duplicated(anno$ENSEMBL), ]
        rownames(anno) <- rownames(df)
    }
    
    cnt <- counts(dds[rownames(df)], normalized=TRUE)
    output <- cbind(anno, as.data.frame(df), as.data.frame(cnt))
    output <- output[order(output$padj, decreasing=FALSE), ]
    write.csv(output, file=filename)
    output 
}

simpleDESeqReport.EnsDb <- function(dds, null.hypothesis.lfcThreshold=0,
                              thres.padj=NULL, thres.lfc=NULL, filename=NULL) {
     ### dds is the returned values from DESeq
    require(DESeq2)

    df <- results(dds, lfcThreshold=null.hypothesis.lfcThreshold)
    if (!is.null(thres.padj))
        df <- subset(df, df$padj < thres.padj)
    
    if (!is.null(thres.lfc)) {
        thres.lfc <- max(0, thres.lfc) # in case it is negative
        df <- subset(df, abs(df$log2FoldChange) > thres.lfc)
    }

    df <- df[order(df$padj, decreasing=FALSE), ]

    anno <- mcols(rowRanges(dds[rownames(df)]))[ c("gene_id", "gene_name", "gene_biotype")]
    cnt <- counts(dds[rownames(df)], normalized=TRUE)
    output <- cbind(anno, df, as(cnt, "DataFrame"))

    if (!is.null(filename)) write.csv(as.data.frame(output), file=filename)
    
    output 
}

simpleLog10Scatter <- function(dds) {
    require(ggplot2)
}
