makeSEwrapper <- function(bamFiles, pkgDir, OrgDb,
                          TxDb=NULL, features=NULL, workers=1L,
                          figDir=NULL,
                          pheno_type=NULL, title="Project",
                          plot.PCA=TRUE, gene_id_type=c("ENTREZID", "ENSEMBL"),
                          do.DESeq=FALSE, lfcThreshold=0,
                          print.report=FALSE) {
    require(DESeq2)
    require(GenomicFeatures)
    require(BiocParallel)
    require(Rsamtools)
    require(ggplot2)
    require(GenomicAlignments)
    require(pheatmap)
    require(RColorBrewer)
    mparam <- MulticoreParam(workers = workers, type = "SOCK")
    register(mparam, default = TRUE)
    registered()
    

    #' check if all bamFiles and package directory exists
    if (!all(file.exists(bamFiles)))
        stop("Some of the bam files do not exists.")

    if (!file.exists(pkgDir))
        stop(pkgDir, " does not exist")

    #' TxDb and features cannot both be NULL
    if (is.null(TxDb) & is.null(features))
        stop("TxDb and features cannot be both NULL")

    #' if TxDb is not null, features will be replaced by exonsBy
    if (!is.null(TxDb))
        features <- GenomicFeatures::exonsBy(TxDb, by="gene")

    #' create inst and data directory if not yet created
    instDir <- file.path(pkgDir, "inst")
    if (!file.exists(instDir)) system(paste("mkdir", instDir))
    dataDir <- file.path(pkgDir, "data")
    if (!file.exists(dataDir)) system(paste("mkdir", dataDir))
    
    #' Make summarizedExperiment object
    message("Read counts and make summarizedExperiment object")
    bamFiles <- Rsamtools::BamFileList(bamFiles)
    se <- GenomicAlignments::summarizeOverlaps(features=features,
                                               reads=bamFiles,
                                               mode="IntersectionStrict",
                                               inter.feature=TRUE,
                                               singleEnd=TRUE,
                                               ignore.strand=TRUE,
                                               BPPARAM=mparam)
    colnames(se) <- sub(".bam", "", colnames(se))
        
    #' Get sample Information
    message("Get sample information")
    sampleInfo <- data.frame(sample_name=colnames(se),
                             file_bam=path(bamFiles),
                             stringsAsFactors=FALSE)
    si <- SGSeq::getBamInfo(sampleInfo, cores=workers)
    rownames(si) <- colnames(se)
    colData(se) <- as(si, "DataFrame")
        
    if (!is.null(pheno_type)) {
        if (!is.factor(pheno_type)) {
            pheno_type <- factor(pheno_type)
        }
        se$pheno_type <- pheno_type
    }
    colnames(se) <- sub(".bam", "", colnames(se))
    ##save(se, file=file.path(dataDir, "se.rda"))
    
    #' Visualization of sample distance by PCA
    if (plot.PCA) {
        message("Plot sample distance by PCA and clustering.")
        if (is.null(figDir)) system(paste("makdir", file.path(pkgDir, "figures")))
        if (!file.exists(figDir)) system(paste("makdir", file.path(pkgDir, "figures")))
        
        dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), ], design = ~ pheno_type)
        rlg <- rlog(dds)
        data <- plotPCA(rlg, intgroup=c("pheno_type"), returnData=TRUE) 
        percentVar <- round(100 * attr(data, "percentVar"))

        ##png(file.path(figDir, paste0(title, "_SampleDistancePCAplot.png")))
        qplot(PC1, PC2, color=pheno_type, data=data) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    labs(list(title=title))
        ggsave(file=file.path(figDir, paste0(title, "_SampleDistancePCAplot.png")))
        ##dev.off()

        ## heatmap
        sampleDists <- dist(t(assay(rlg)))
        colors <- colorRampPalette(brewer.pal(9, "RdYlBu"))(255)
        sampleDistMatrix <- as.matrix(sampleDists)
        #colnames(sampleDistMatrix) <- rownames(sampleDistMatrix) <-
        #    paste0(colnames(rlg), ":", rlg$pheno_type)
        fname <- file.path(figDir, paste0(title, "_SampleDistanceHeatmap.png"))
        annotation_col <- data.frame(pheno_type=factor(rlg$pheno_type))
        rownames(annotation_col) <- colnames(rlg)
        pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists, border_color=NA,
                 treeheight_row=0, main=paste0(title, " - sample distance"),
                 annotation_col=annotation_col,
                 col=colors, silent=TRUE, filename=fname)
    }


    if (do.DESeq) {
        if (is.null(pheno_type)) stop("The formula (pheno_type) for DESeq is not given.")
        
        dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), ], design = ~ pheno_type)
        message("DESeq (~ pheno_type) ...")
        dds <- DESeq(dds)
        ##save(dds, file=file.path(dataDir, "dds.rda"))
        if (print.report) {
            statsDir <- file.path(pkgDir, "inst", "stats")
            filename <- file.path(statsDir, paste0(title,"_results_dds.csv"))
            message("Generate report: ", filename)
            if (!file.exists(statsDir)) system(paste("mkdir", statsDir))
            res <- simpleDESeqReport(dds=dds, OrgDb=OrgDb, lfcThreshold=lfcThreshold,
                                     gene_id_type=gene_id_type,
                                     filename=filename)
        }
    }
    
    if (!do.DESeq) dds <- "No DESeq result returned"
    return(list(se=se, dds=dds))
}
