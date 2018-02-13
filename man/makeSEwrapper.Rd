\name{makeSEwrapper}
\alias{makeSEwrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make SummarizedExpriment instance wrapper
}
\description{
A wrapper funtion that performed several tasks including counting hit, visualization of sample distance, DESeq analysis and report with annotation.
}
\usage{
makeSEwrapper(bamFiles, pkgDir, OrgDb, TxDb = NULL, features = NULL,
              workers = 1L, pheno_type = NULL,
              title = "Project", plot.PCA = TRUE,
              gene_id_type = c("ENTREZID", "ENSEMBL"),
	      do.DESeq = FALSE, lfcThreshold = 0,
	      print.report = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
%%     ~~Describe \code{bamFiles} here~~
}
  \item{pkgDir}{
%%     ~~Describe \code{pkgDir} here~~
}
  \item{OrgDb}{
%%     ~~Describe \code{OrgDb} here~~
}
  \item{TxDb}{
%%     ~~Describe \code{TxDb} here~~
}
  \item{features}{
%%     ~~Describe \code{features} here~~
}
  \item{workers}{
%%     ~~Describe \code{workers} here~~
}
  \item{pheno_type}{
%%     ~~Describe \code{pheno_type} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{plot.PCA}{
%%     ~~Describe \code{plot.PCA} here~~
}
  \item{gene_id_type}{
%%     ~~Describe \code{gene_id_type} here~~
}
  \item{do.DESeq}{
%%     ~~Describe \code{do.DESeq} here~~
}
  \item{lfcThreshold}{
%%     ~~Describe \code{lfcThreshold} here~~
}
  \item{print.report}{
%%     ~~Describe \code{print.report} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bamFiles, pkgDir, OrgDb, TxDb = NULL, features = NULL, 
    workers = 1L, pheno_type = NULL, title = "Project", plot.PCA = TRUE, 
    gene_id_type = c("ENTREZID", "ENSEMBL"), do.DESeq = FALSE, 
    lfcThreshold = 0, print.report = FALSE) 
{
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
    if (!all(file.exists(bamFiles))) 
        stop("Some of the bam files do not exists.")
    if (!file.exists(pkgDir)) 
        stop(pkgDir, " does not exist")
    if (is.null(TxDb) & is.null(features)) 
        stop("TxDb and features cannot be both NULL")
    if (!is.null(TxDb)) 
        features <- GenomicFeatures::exonsBy(TxDb, by = "gene")
    instDir <- file.path(pkgDir, "inst")
    if (!file.exists(instDir)) 
        system(paste("mkdir", instDir))
    dataDir <- file.path(pkgDir, "data")
    if (!file.exists(dataDir)) 
        system(paste("mkdir", dataDir))
    message("Read counts and make summarizedExperiment object")
    bamFiles <- Rsamtools::BamFileList(bamFiles)
    se <- GenomicAlignments::summarizeOverlaps(features = features, 
        reads = bamFiles, mode = "IntersectionStrict", inter.feature = TRUE, 
        singleEnd = TRUE, ignore.strand = TRUE, BPPARAM = mparam)
    colnames(se) <- sub(".bam", "", colnames(se))
    message("Get sample information")
    sampleInfo <- data.frame(sample_name = colnames(se), file_bam = path(bamFiles), 
        stringsAsFactors = FALSE)
    si <- SGSeq::getBamInfo(sampleInfo, cores = workers)
    colData(se) <- as(si, "DataFrame")
    if (!is.null(pheno_type)) {
        if (!is.factor(pheno_type)) {
            pheno_type <- factor(pheno_type)
        }
        se$pheno_type <- pheno_type
    }
    colnames(se) <- sub(".bam", "", colnames(se))
    if (plot.PCA) {
        message("Plot sample distance by PCA and clustering.")
        figDir <- file.path(pkgDir, "inst", "figures")
        if (!file.exists(figDir)) 
            system(paste("mkdir", figDir))
        dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), 
            ], design = ~pheno_type)
        rlg <- rlog(dds)
        data <- plotPCA(rlg, intgroup = c("pheno_type"), returnData = TRUE)
        percentVar <- round(100 * attr(data, "percentVar"))
        qplot(PC1, PC2, color = pheno_type, data = data) + xlab(paste0("PC1: ", 
            percentVar[1], "\% variance")) + ylab(paste0("PC2: ", 
            percentVar[2], "\% variance")) + labs(list(title = title))
        ggsave(file = file.path(figDir, paste0(title, "_SampleDistancePCAplot.png")))
        sampleDists <- dist(t(assay(rlg)))
        colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
        sampleDistMatrix <- as.matrix(sampleDists)
        fname <- file.path(figDir, paste0(title, "_SampleDistanceHeatmap.png"))
        annotation_col <- data.frame(pheno_type = factor(rlg$pheno_type))
        rownames(annotation_col) <- colnames(rlg)
        pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, 
            clustering_distance_cols = sampleDists, border_color = NA, 
            treeheight_row = 0, main = paste0(title, " - sample distance"), 
            annotation_col = annotation_col, col = colors, silent = TRUE, 
            filename = fname)
    }
    if (do.DESeq) {
        if (is.null(pheno_type)) 
            stop("The formula (pheno_type) for DESeq is not given.")
        dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), 
            ], design = ~pheno_type)
        message("DESeq (~ pheno_type) ...")
        dds <- DESeq(dds)
        if (print.report) {
            statsDir <- file.path(pkgDir, "inst", "stats")
            filename <- file.path(statsDir, paste0(title, "_results_dds.csv"))
            message("Generate report: ", filename)
            if (!file.exists(statsDir)) 
                system(paste("mkdir", statsDir))
            res <- simpleDESeqReport(dds = dds, OrgDb = OrgDb, 
                lfcThreshold = lfcThreshold, gene_id_type = gene_id_type, 
                filename = filename)
        }
    }
    if (!do.DESeq) 
        dds <- "No DESeq result returned"
    return(list(se = se, dds = dds))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
