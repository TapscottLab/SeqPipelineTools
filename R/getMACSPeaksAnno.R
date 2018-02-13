#' these scripts are used to get the annotation of the peaks called by macs2
#' DATE: 2/10/16 by Chao-Jen Wong
#'
#' Scripts:
#'      getMACSPeaksAnn


dfToGr <- function(df) {
    ## convert MACS' peak.xls to GRanges
    gr <- GRanges(seqnames=df$chr,
                  ranges=IRanges(start=df$start, end=df$end))
    mcols(gr) <- df[5:10]
    gr
}


.getBestHit <- function(query, subject) {
    ## We want one query hits ony one subject. If there is 1-to-n mapping
    ## where n > 1, the function picks the best subject based upon the overlapping
    ## ranges. The function does allow n-to-1 mapping.
    ol <- findOverlaps(query, subject, ignore.strand=TRUE)
    w <- width(pintersect(query[queryHits(ol)], subject[subjectHits(ol)]))
    best <- tapply(w, queryHits(ol), function(x) {
        keep <- rep(FALSE, length(x))
        keep[which.max(x)[1]] <- TRUE
        keep
    })
    best <- unlist(best)
    ol.best <- ol[best]
}
    
olRMSK <- function(peaks.gr, rmsk) {
    ## what if it is not one-to-one mapping? pick one
    ## base upon the range of overlap
    ol <- .getBestHit(query=peaks.gr, subject=rmsk)
    df <- DataFrame(overlap.LTR=rep(FALSE, length(peaks.gr)),
                    repRange=rep(NA, length(peaks.gr)),
                    repName=rep(NA, length(peaks.gr)),
                    repClass=rep(NA, length(peaks.gr)),
                    repFamily=rep(NA, length(peaks.gr)))
    rng <-  rmsk[subjectHits(ol)]
    df$repRange[queryHits(ol)] <- as.character(rng)
    df$repName[queryHits(ol)] <- as.character(rng$repName)
    df$repClass[queryHits(ol)] <- as.character(rng$repClass)
    df$repFamily[queryHits(ol)] <- as.character(rng$repFamily)
    df$overlap.LTR[queryHits(ol)] <- df$repClass[queryHits(ol)] == "LTR"
    df
}

## add one column at the time to a data frame 
addColToDF <- function(df, addColumn, n) {
    ## one column at a time
    df[, names(addColumn)] <- rep(NA, nrow(df))
    df[n, names(addColumn)] <- as.character(addColumn[, names(addColumn)])
    df
}

getMACSPeaksAnno <- function(peakFile, genome="mm", txdb=NULL,
                             annoDb=NULL, rmsk=NULL, skip=26,
                             annotate.repeat=TRUE) {
    ## only support mm now
    #require(ChIPseeker)
    #require(rtracklayer)
    #require(xlsx)
    
    if (genome == "mm" & is.null(txdb))  {
        message("The genome is default to mm10")
        require(TxDb.Mmusculus.UCSC.mm10.knownGene)
        require(org.Mm.eg.db)
        annoDb <- "org.Mm.eg.db"
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
        rmsk <- get(load("~/tapscott/mm10/mm10.rmsk.rda"))
    }

    if (is.null(txdb)) stop("txdb is missing ...")
    #if (is.null(annoDb)) stop("annoDb is missing ...")
    
    message(basename(peakFile))
    peaks <- read.csv(peakFile, sep="\t", skip=skip, stringsAsFactors=FALSE)
    ##peaks.gr <- dfToGr(peaks)
    peaks.gr <- as(peaks, "GRanges")
    message("Total ", length(peaks.gr), " peaks.")
    ## append annotation: the nrow != nrow(peaks.gr)
    peakAnno <- ChIPseeker::annotatePeak(peaks.gr, tssRegion=c(-3000, 3000),
                                         TxDb=txdb, annoDb=annoDb)
    tmp <- as.data.frame(peakAnno)
    n <- strsplit(tmp$name, "_")
    l <- elementNROWS(n)[1]
    n <- as.numeric(sapply(n, "[[", l)) ## the lable of peaks
    
    df <- mcols(peaks.gr)
    for (i in c(13:ncol(tmp))) ## add peakAnno to mcols. peakAnno has missing rows
        df <- addColToDF(df=df,  n=n, addColumn=tmp[, i, drop=FALSE])

    mcols(peaks.gr) <- as(df, "DataFrame")

    if (annotate.repeat) {
        anno_rmsk <- olRMSK(peaks.gr, rmsk)
        mcols(peaks.gr) <- append(mcols(peaks.gr), anno_rmsk)
    }
    
    peaks.gr
}


    



