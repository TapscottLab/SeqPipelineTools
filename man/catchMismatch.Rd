\name{catchMismatch}
\alias{catchMismatch}

\title{
Catch mismatch positions from BAM files
}
\description{

}
\usage{
catchMismatch(bam_file, BSgenome, snps = NULL,
    BSgenomeInjectSNP = NULL,
    include_deletions = FALSE, include_insertions = FALSE,
    min_base_quality = 10, min_nucleotide_depth = 1,
    min_mapq = 10, max.mismatch = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_file}{A character string specifying the bam file.
}
  \item{BSgenome}{BSGenome object of the reference genome.
}
  \item{snps}{A GRange object of SNP position on the reference genome.
}
  \item{BSgenomeInjectSNP}{BSgenome object injected with SNP. See details.
}
  \item{include_deletions}{A logical specifying whether the pileup should include deletion.
}
  \item{include_insertions}{A logical specifying whether the pileup should include insertions.
}
  \item{min_base_quality}{A locial indicating the minimum base quality score.
}
  \item{min_nucleotide_depth}{An integer indicating the mininum depth of nucleotide.
}
  \item{min_mapq}{An integer indicating the mininum MAPQ score.
}
  \item{max.mismatch}{an integer indicating the max. mismatch of a read.
}
}
\details{
This function generates a data.frame of mismatch nucleotides and the type of mismatch relative to the reference genome. See example.
}
\value{A data.frame includes columns of mismatch position, nucleotides, nucleotides depth, reference nucleotides and mismatch type.
}

\author{
Chao-Jen Wong <cwon2@fredhutch.org>
}

\seealso{
\code{\link{tallyWrapper}}
}
\examples{
\dontrun{
mainDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.dsRNA"
bamDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/Sean/MB135_dsRNA/hg38_bwa"
bam_file <- file.path(bamDir, "Sample_mb135_wDox_k1.sorted.bam")

library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
library(rtracklayer)
library(Rsamtools)
snps <- get(load(file="/fh/fast/tapscott_s/CompBio/hg38/FDb.UCSC.snp147Common.hg38/data/snp147common.GR.rda"))
BSgenome <- BSgenome.Hsapiens.UCSC.hg38
tally <- catchMismatch(bam_file=bam_file, BSgenome=BSgenome, snps=snps,
                       min_base_quality=10,
                       min_nucleotide_depth=2,
                       include_deletions=FALSE,
                       include_insertions=FALSE,
                       min_mapq=10, max.mismatch=10)
SeqPipelineTools:::freqAIMismatch(tally, exclude.SNP=FALSE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
