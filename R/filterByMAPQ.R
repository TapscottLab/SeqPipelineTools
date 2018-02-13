filterByMAPQ <- function(bam_files, destination, mapq=12, cores=1) {
    ## Rsamtools::bamMapqFilter or (SacnBamParam(..., mapqFilter=12))
    ## This filter is for filtering bad quality reads (mapq < 12) from BAM files
    require(Rsamtools)
    require(parallel)
            
    if (!all(file.exists(dirname(destination))))
        stop("Some destination directory does not exist.")

    if (length(bam_files) != length(destination))
        stop("The lengths of bam_files and destination files must equal to each other")

    filterLowQ <- function(read) {
        read$mapq > mapq
    }

    params <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                           what=c("qname", "mapq"),
                           reverseComplement=TRUE, simpleCigar=TRUE)
    mcmapply(filterBam, file=bam_files, destination=destination,
             MoreArgs=list(param=params, filter=FilterRules(filterLowQ)),
             mc.preschedule=FALSE, SIMPLIPY=FALSE, mc.cores=cores)
    ## return reports?? how many reads have been filtered? 
    invisible() 
}
