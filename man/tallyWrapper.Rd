\name{tallyWrapper}
\alias{tallyWrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A wrapper for catching mismatch nucleotides
}
\description{

}
\usage{
tallyWrapper(bamFiles, genome = "hg38", destination = ".", cores = 1L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
%%     ~~Describe \code{bamFiles} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{destination}{
%%     ~~Describe \code{destination} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
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
See Also as \code{\link{catchMismatch}}
}

\examples{
\dontrun{
cores <- 4L
mainDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.dsRNA"
bamDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/Sean/MB135_dsRNA/hg38_bwa"
bamFiles <- list.files(bamDir, full.name=TRUE, pattern=".bam$")[2:5]
source("/fh/fast/tapscott_s/CompBio/R_package/SeqPipelineTools/R/mismatchTools.R"
)

tally <- tallyWrapper(bamFiles, destination=file.path(mainDir, "test"),
                      cores=cores)
lapply(tally, SeqPipelineTools:::freqAIMismatch, exclude.SNP=TRUE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
