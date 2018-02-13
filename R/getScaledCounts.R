#' This convertion work only for single-ended read
#' 
getScaledCountsPerTx <- function(SE,
    features = rowRanges(SE),
    counts = assays(SE)$counts,
    read_length = colData(SE)$read_length,
    lib_size = colData(SE)$lib_size)
{
    if (is.null(lib_size)) stop("Lib size is essential")
    if (is.null(read_length)) stop("Read length is essential")
    if (is.null(features)) stop("Feature range is essential")
    
    n_features <- length(features)
    n_samples <- ncol(SE)

    E <- matrix(NA_real_, nrow = n_features, ncol = n_samples)
    
    R <- rep(read_length, each=n_features)
    L <- rep(sum(width(features)), n_samples)
    E <- matrix(L + R - 1, nrow = n_features, ncol = n_samples)
    X <- counts / (E * 1e-3)
    X <- sweep(X, 2, lib_size * 1e-6, FUN = "/")

    return(X)

}
