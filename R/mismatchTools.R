#.fetchSNPs <- function(injectedGenome, tally) {
#    ## get nucleotide form genome injected with snps
#    tally$seqnames <- factor(tallyf$seqnames)
#    gr <- GRanges(seqnames=tally$seqnames,
#                  IRanges(start=tally$pos, width=1))
#    ref_w_snp <- getSeq(injectedGenome, tally$seqnames,
#                        start=tally$pos, width=1)
#    
#             }

#.fetchSNPs <- function(snps, tally){
#    gr <- GRanges(seqnames=tally$seqnames,
#                      IRanges(start=tally$pos, width=1))
#    seqlevels(gr) <- gsub("chr", "ch", seqlevels(gr))
#    message("Finding SNPs overlap with mismatched sites")
#    snps_gr <- snpsBySeqname(snps, seqlevels(snps))
#    seqlevels(snps_gr) <- gsub("ch", "chr", seqlevels(snps_gr))
#    ov <- findOverlaps(gr, snps_gr, ignore.strand=TRUE)
    
    #snps_ov <- snpsByOverlaps(snps, range=gr, minoverlap=1L)
    #seqlevels(snps_ov) <- gsub("ch", "chr", seqlevels(snps_ov))
#    df <- data.frame(SNP=rep(FALSE, nrow(tally)),
#                     alleles_as_ambig=rep(NA, nrow(tally)))
#}
    
.getTally <- function(pup, genome) {
    ## pup: pileup
    ## BSgenome
    pup$seqnames <- factor(pup$seqnames)
    tally <- lapply(split(pup, pup$seqnames), function(x) {
        if (identical(nrow(x), 0L)) return(x)
        
        CHR <- as.character(x$seqnames)[1]
        message(CHR)
        rf <- DNAStringSet(genome[[CHR]],
                           start=as.integer(x$pos), width=1)
        x$ref <- as.character(rf)
        x[x$nucleotide != x$ref, ]
    })
    tally <- do.call(rbind, unname(tally))
}

.getTally.advance <- function(pup, genome) {
    ## pup: pileup
    ## BSgenome
    pup$seqnames <- factor(pup$seqnames)
    tally <- lapply(split(pup, pup$seqnames), function(x) {
        if (identical(nrow(x), 0L)) return(x)
        
        CHR <- as.character(x$seqnames)[1]
        message(CHR)
        rf <- DNAStringSet(genome[[CHR]],
                           start=as.integer(x$pos), width=1)
        x$ref <- as.character(rf)
        #x[x$nucleotide != x$ref, ]
        #' add percentage
        tmp <- x
        tmp$name <- rownames(x)
        tmp_list <- tapply(tmp$pos, tmp$count, function(cnt) {
            cnt / sum(cnt)
        })
        tmp <- unlist(tmp_list, use.name=FALSE)
        #' reorder to the same order as x
        tmp <- tmp[order(as.integer(tmp$name), decreasing=FALSE), ]
        
    })
    tally <- do.call(rbind, unname(tally))
}

catchMismatch <- function(bam_file, BSgenome, snps=NULL,
                          paired=FALSE,
                          BSgenomeInjectSNP=NULL, include_deletions=FALSE,
                          include_insertions=FALSE,
                          distinguish_strand=FALSE,
                          distinguish_nucleotides=FALSE,
                          min_base_quality=10, min_nucleotide_depth=1,
                          min_mapq=10, max.mismatch=10) {
    require(rtracklayer)
    require(Rsamtools)
    message(bam_file)
    params <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                           tagFilter=list(NM=c(1:max.mismatch)),
                           what=c("qname", "seq"), tag=c("MD", "NM"))
    pp <- PileupParam(min_base_quality=min_base_quality,
                      min_mapq=min_mapq,
                      min_nucleotide_depth=min_nucleotide_depth,
                      max_depth=1000,
                      distinguish_strand=distinguish_strand,
                      distinguish_nucleotides=distinguish_nucleotides,
                      include_deletions=include_deletions,
                      include_insertions=include_insertions)
    
    pup <- pileup(bam_file, scanBamParam=params,
                  pileupParam=pp)
    
    ## sanitize check: seq levels of pup and BSgenome
    pup$seqnames <- factor(pup$seqnames)
    tally <- .getTally(pup=pup, genome=BSgenome)
    tally$mismatch_type <- paste0(tally$ref, tally$nucleotide)

    ## SNP locs: (1) change the seqnames to UCSC standard
    ## (2) overlap with mismatch sides
    if (!is.null(snps)) {
        tally$SNP <- FALSE
        gr <- GRanges(seqnames=tally$seqnames,
                      IRanges(start=tally$pos, width=1))
        message("Finding SNPs overlap with mismatched sites")
        ov <- findOverlaps(gr, snps, minoverlap=1L, ignore.strand=TRUE)
        tally$SNP[queryHits(ov)] <- TRUE 
    }
    tally
}


freqAIMismatch <- function(tally, exclude.SNP=TRUE) {
    message("The AI mismatch includes A-G and T-C.")
    if (exclude.SNP) {
        tally <- subset(tally, !tally$SNP)
    }
    tt <- table(tally$mismatch_type)
    (tt["AG"]+tt["TC"])/sum(tt)
}

tallyWrapper <- function(bamFiles, genome="hg38", destination=".",
                         cores=1L) {

    require(parallel)
    require(rtracklayer)
    require(Rsamtools)
    if (genome=="hg38") {
        message("loading snp147 common track")
        require(BSgenome.Hsapiens.UCSC.hg38)
        BSgenome <- BSgenome.Hsapiens.UCSC.hg38
        snps <- get(load(file="/fh/fast/tapscott_s/CompBio/hg38/FDb.UCSC.snp147Common.hg38/data/snp147common.GR.rda"))
    } else {stop("does not support non-hg38 genome")}

    tally <- mclapply(bamFiles, catchMismatch,
                  BSgenome=BSgenome, snps=snps,
                  min_base_quality=10, min_nucleotide_depth=2,
                  min_mapq=10, max.mismatch=10,
                  mc.cores=cores, mc.preschedule=FALSE)
    names(tally) <- sub(".bam", "", basename(bamFiles))
    
    ## save tally to a csv file
    if (!file.exists(destination)) {
        message(destination, "does not exist.")
        return(tally)
    }

    lapply(names(tally), function(x) {
        filename <- file.path(destination, paste0("tally_", x, ".csv"))
        message("Exporting ", filename)
        write.csv(tally[[x]], file=filename)
    })

    tally
}
    
