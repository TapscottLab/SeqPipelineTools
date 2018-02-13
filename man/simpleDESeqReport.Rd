\name{simpleDESeqReport}
\alias{simpleDESeqReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple DESeq report
}
\description{
This function generates a simple CSV file of the DESeq analysis results and annotation.
}
\usage{
simpleDESeqReport(dds, OrgDb, lfcThreshold = 0,
    gene_id_type = c("ENTREZID", "ENSEMBL"),
    thres.padj = NULL, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dds}{A `DESeqDataSet` object returned by DESeq.}
  \item{OrgDb}{An orgDb object in which the annotatio is obtained.}
  \item{lfcThreshold}{A numeric specifying the log fold change that the hypothesis tests
  are based upon. }
  \item{gene_id_type}{A character string that is either "ENTREZID" or "ENSEMBL".}
  \item{thres.padj}{Default to NULL. If NULL, returns statistical results and annotation
  for the whole set of genes. If given a fractional number, the results will be
  truncated accordingly.}
  \item{filename}{A character string specifying the destination and the filename of the
  return. The extension is preferred to be ".csv".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A DataFrame instance containing the results of DESeq with annotation provided by the organic database, OrgDb.
}

\author{
Chao-Jen Wong
}


\seealso{
\code{\link{results}}
}

\examples{
\dontrun{
library(org.Hs.eg.db)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
dds <- get(load(file.path(pkgDir, "data", "HinH.ens.dds.rda")))
filename <- file.path(pkgDir, "inst", "stats", "testing.csv")
df <- simpleDESeqReport(dds, OrgDb=org.Hs.eg.db, gene_id_type="ENTREZID",
    lfcThreshold=1, thres.padj=0.1, 
    filename=filename)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
