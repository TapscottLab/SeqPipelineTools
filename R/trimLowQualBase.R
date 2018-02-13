trimLowQualBase <- function(fls, output_name=NULL, destDir,
                            CleanNFilter=FALSE) {
    #' The main function is to trim low-quality reads from the right end using a
    #' sliding window with 5 nucleotides and only keep the reads that's longer than
    #' 36L long
    #' nucleoties.
    #' Use:
    #' (1) trimTailw(x, 2, "4", 2)
    #'  Has option to remove reads from FASTQ with at least one N in the sequence
    #' (2) nFilter(): filter reads containing N in the nucleotides
    #' (3) retain reads longer than 36N.

    #' fls: character or a vector of characters of FASTQ file name
    #' The base quality threshold is set to "4", equivalent to probability of (Q=19)
    #' 1% that the base is incorrect (for Phred+33 standard).

    ## must be one file at a time? 

    ## require library
    require(ShortRead)
    #' sanity check
    if (!file.exists(destDir)) stop(destDir, " does not exist!")
    if (!all(file.exists(fls))) stop(fls, " does not exist!")
    
    message("Trimming ", fls)
    fq <- ShortRead::readFastq(fls)
    org_length <- length(fq)
    #' filter reads containg N in the nucleotides
    if (CleanNFilter)  fq <- fq[ShortRead::nFilter()(fq)]
    
    #'  trim tail with 2 to 5 low quality (4) alignment
    fq <- ShortRead::trimTailw(fq, 2, "4", 2)
    fq <- fq[width(fq) >= 36]

    if (is.null(output_name)) fname <- file.path(destDir, basename(fls)[1])
    else fname <- file.path(destDir, output_name)

    if (file.exists(fname)) system(sprintf("rm %s", fname))

    message("Exporting ", fname)
    ShortRead::writeFastq(fq, file=fname)
    c(org_length=org_length, length=length(fq), wrange=range(width(fq)))
}
