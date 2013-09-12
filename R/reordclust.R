reordclust <- function (clustering,from,to) 
{
    if (inherits(clustering,c('clustering','partana','slice'))) clustering <- clustering$clustering
    if (is.numeric(clustering)) {
        if (min(clustering)< 0 || (length(table(clustering)) != max(clustering))) {
            cat('WARNING: renumbering clusters to consecutive integers\n')
            clustering <- match(clustering,sort(unique(clustering)))
        }
    }

    nfrom <- length(table(from))
    nto <- length(table(to))
    if (length(from)!=length(to)) stop("membership vectors must be the same size")

    out <-rep(0,length(clustering))
    for (i in 1:nfrom) {
        out[clustering==from[i]] <- to[i]
    }
    out <- as.numeric(factor(out))
    out <- list(clustering=out)
    class(out) <- 'clustering'
    out
}
