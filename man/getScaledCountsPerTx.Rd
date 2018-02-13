\name{getScaledCountsPerTx}
\alias{getScaledCountsPerTx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getScaledCountsPerTx(SE, features = rowRanges(SE), counts = assays(SE)$counts, read_length = colData(SE)$read_length, lib_size = colData(SE)$lib_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SE}{
%%     ~~Describe \code{SE} here~~
}
  \item{features}{
%%     ~~Describe \code{features} here~~
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{read_length}{
%%     ~~Describe \code{read_length} here~~
}
  \item{lib_size}{
%%     ~~Describe \code{lib_size} here~~
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
function (SE, features = rowRanges(SE), counts = assays(SE)$counts, 
    read_length = colData(SE)$read_length, lib_size = colData(SE)$lib_size) 
{
    if (is.null(lib_size)) 
        stop("Lib size is essential")
    if (is.null(read_length)) 
        stop("Read length is essential")
    if (is.null(features)) 
        stop("Feature range is essential")
    n_features <- length(features)
    n_samples <- ncol(SE)
    E <- matrix(NA_real_, nrow = n_features, ncol = n_samples)
    R <- rep(read_length, each = n_features)
    L <- rep(sum(width(features)), n_samples)
    E <- matrix(L + R - 1, nrow = n_features, ncol = n_samples)
    X <- counts/(E * 0.001)
    X <- sweep(X, 2, lib_size * 1e-06, FUN = "/")
    return(X)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
