mergeclust <- function(clustering,from,to) 
{
    if (inherits(clustering,c('clustering','partana','slice'))) clustering <- clustering$clustering
    if (is.numeric(clustering)) {
        if (min(clustering)< 0 || (length(table(clustering)) != max(clustering))) {
            cat('WARNING: renumbering clusters to consecutive integers\n')
            clustering <- match(clustering,sort(unique(clustering)))
        }
    }

    clustering[clustering==from] <- to
    out <- list()
    out$clustering <- as.numeric(factor(clustering))
    out
}
