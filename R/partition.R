partition <- function(x, ...)
{
    UseMethod("partition")
}

partition.partana <- function (x,dist,...)
{
    if (!inherits(x,'partana')) {
        stop("You must supply an object of class partana")
    }
    if (class(dist) != 'dist') {
        stop("You must specify an object of class dist as the second argument")
    }
    out <- list()
    attr(out,"call") <- match.call()
    out$dist <- dist
    out$clustering <- x$clustering
    out$silinfo <- silhouette(x$clustering,dist)
    attr(out,'class') <- 'partition'
    return(out)
}

partition.clustering <- function (x, dist, ...)
{
    if (!inherits(x,'clustering')) {
        stop("You must supply an object of class clustering as the first argument")
    }
    if (class(dist) != 'dist') {
        stop("You must supply an object of class dist as the second argument")
    }
    out <- list()
    attr(out,"call") <- match.call()
    out$dist <- dist
    out$clustering <- x$clustering
    out$silinfo <- silhouette(x$clustering,dist)
    attr(out,'class') <- 'partition'
    return(out)
}
