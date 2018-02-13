\name{trimLowQualBase}
\alias{trimLowQualBase}

\title{
Trim Low Quality Base and more
}
\description{
This function is a wrapper performing three major preprocessing steps: trimming low
quality base, filtering reads containing N in the nucleotides (optional) and
filtering reads less than 36N long. 
}

\usage{
trimLowQualBase(fls, output_name = NULL, destDir, CleanNFilter = FALSE)
}

\arguments{
  \item{fls}{Character strings specifying names and location of the FASTQ files.}
  \item{output_name}{A character string specifying the resulting FASTQ file.}
  \item{destDir}{A character string specifying the destination of the resulting FASTQ
  file.}
  \item{CleanNFilter}{A logical indicating whether to filter reads containing N in the
  nucleotides.}
}

\details{
This wrapper uses "ShortRead::trimTailw(fq, 2, "4", 2)" to remove low-quality tails when 2 consecutive bases from the right of a 5-nucleotide window fall lower than "4", correspoding to Q=19 (one percent of the chance that the base is incorrect). "ShortRead::nFilter()" is used to filter reads containing N in the nucleotides. This filter is optional. 
}

\value{
A FASTQ file is saved in the designated location. The return value is a vector indicating the number of reads before and after the trimming and filtering process as well as the resulted range of reads length. 
}

\author{
Chao-Jen Wong
}


\seealso{
\code{\link{trimTailw}}, \code{\link{nFilter}}
}
\examples{

\dontrun{
ngsDir <- "/shared/ngs/illumina/jwhiddon/160122_D00300_0236_AHJ7JNBCXX"
fqDir <- file.path(ngsDir,"Unaligned", "Project_jwhiddon")
pkgDir <- "/fh/fast/tapscott_s/CompBio/ChIP-Seq/mm10.HuMoDoxy"

## step by step
fqSampleDir <- dir(fqDir, full=TRUE)[1]
fqFiles <- list.files(fqSampleDir, pattern="\\.fastq.gz$", full.names=TRUE)
output_name <- paste0(basename(fqSampleDir), ".fastq.gz")
destDir <- file.path(pkgDir, "inst", "filtered.fq")
l <- trimLowQualBase(fls=fqFiles, output_name=output_name, destDir=destDir,
    CleanNFilter=TRUE)
}
\dontrun{
## do all the samples
library(parallel)
fqSampleDir <- dir(fqDir, full=TRUE)
l <- mclapply(fqSampleDir, function(x) {
    fls <- list.files(pattern="\\.fastq.gz$", full.names=TRUE)
    output_name <- paste0(basename(x), ".fastq.gz")
    destDir <- file.path(pkgDir, "inst", "filtered.fq")
    trimLowQualBase(fls=fls, output_name=output_name, destDir=destDir,
        CleanNfilter=TRUE)
}, mc.cores=2L, mc.preschedule=FALSE)
l
}


## The function is currently defined as
function (fls, output_name = NULL, destDir, CleanNFilter = FALSE) 
{
    require(ShortRead)
    if (!file.exists(destDir)) 
        stop(destDir, " does not exist!")
    if (!file.exists(fls)) 
        stop(fls, " does not exist!")
    message("Trimming ", fls)
    fq <- ShortRead::readFastq(fls)
    org_length <- length(fq)
    if (CleanNFilter) 
        fq <- fq[(ShortRead::nFilter())(fq)]
    fq <- ShortRead::trimTailw(fq, 2, "4", 2)
    fq <- fq[width(fq) >= 36]
    if (is.null(output_name)) 
        fname <- file.path(destDir, basename(fls)[1])
    else fname <- file.path(destDir, output_name)
    if (file.exists(fname)) 
        system(sprintf("rm \%s", fname))
    ShortRead::writeFastq(fq, file = fname)
    c(org_length = org_length, length = length(fq), wrange = range(width(fq)))
  }
}

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
