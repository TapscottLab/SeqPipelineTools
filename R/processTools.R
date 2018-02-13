preprocessFastq <- function(fls, output_name=NULL, destDir) {
    ## this filter remove reads from FASTQ with at least one N in the sequence
    ## and trim tail with 2 to 5 low quality alignment
    ## must be one file at a time
    if (!file.exists(destDir)) stop(destDir, " does not exist!")
    if (!file.exists(fls)) stop(fls, " does not exist!")
    
    message("Processing ", fls)
    require(ShortRead)
    fq <- readFastq(fls)
    org_length <- length(fq)
    ## filter reads containing N and trim tail with 2 to 5 low quality alignment
    fq <- fq[nFilter()(fq)]
    fq <- trimTailw(fq, 2, "4", 2)
    fq <- fq[width(fq) >= 36]
    if (is.null(output_name)) fname <- file.path(destDir, basename(fls))
    else fname <- file.path(destDir, output_name)

    if (file.exists(fname)) system(sprintf("rm %s", fname))

    writeFastq(fq, file=fname)
    c(org_length=org_length, length=length(fq), wrange=range(width(fq)))
}

simpleFilterFastq <- function(fls, out_fname) {
    ## just simply remove read from FASTQ with at least one N in the sequence
    require(Rsamtools)
    require(ShortRead)
    message("Filtering ", fls, "to ", out_fname)
    fq <- readFastq(fls)
    org_length <- length(fq)
    if (file.exists(out_fname)) system(sprintf("rm %s", out_fname))
    fq <- fq[nFilter()(fq)]
    writeFastq(fq, file=out_fname)
    c(org_length=org_length, length=length(fq))
}

useSubread <- function(fq_file) {
}

filterByMAPQ <- function(bam_files, destination,  cores=1) {
    ## This filter is for filtering bad quality reads (mapq < 12) from BAM files
    require(Rsamtools)
    require(parallel)
            
    if (!all(file.exists(dirname(destination))))
        stop("Some destination directory does not exist.")

    if (length(bam_files) != length(destination))
        stop("The lengths of bam_files and destination files must equal to each other")

    filterLowQ <- function(read) {
        read$mapq > 12
    }

    params <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                           what=c("qname", "mapq"),
                           reverseComplement=TRUE, simpleCigar=TRUE)
    mcmapply(filterBam, file=bam_files, destination=destination,
             MoreArgs=list(param=params, filter=FilterRules(filterLowQ)),
             mc.preschedule=FALSE, SIMPLIPY=FALSE, mc.cores=cores)
    
    invisible()
}

.filter_peaks <- function(gr, qthres=0.01, FC=NULL, pileup=NULL,
                         n=NULL) {
    value <- -log10(qthres)
    gr <- gr[gr$X.log10.qvalue. > value, ]

    if (!is.null(FC) & is.numeric(FC))
        gr <- gr[gr$fold_enrichment > FC, ]

    if (!is.null(pileup) & is.numeric(pileup))
        gr <- gr[gr$pileup > pileup, ]

    if (!is.null(n)) {
        gr <- gr[order(gr$X.log10.qvalue., decreasing=TRUE)]
        gr <- gr[1:n]
    }

    gr
}

.aroundSummit <- function(x, L) {
    #' x must be an GRanges; this scripts is unawared of strand
    #' This script expend up/down for L bps.
    Lhalf <- round(L/2)
    start(x) <- start(x) - Lhalf
    end(x) <- end(x) + Lhalf - 1
    x
}

makePeaksSeq <- function(peaks, BS, n=NULL, width=100, seq=TRUE) {
    summits <- .filter_peaks(peaks, n=n)
    start(summits) <- end(summits) <- summits$abs_summit
    peaks <- .aroundSummit(summits, width) ## 50 around the summit
    names(peaks) <- 1:length(peaks)
    if (seq) {
        peaks.seq <- BSgenome::getSeq(BS, peaks)
        names(peaks.seq) <- 1:length(peaks.seq)
        return(peaks.seq)
    }
    else
        return(peaks)
}

rpkm <- function(bam_file) {
}
