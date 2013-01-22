mergeclust <- function(clustering,from,to) 
{
    if (inherits(clustering,c('clustering','partana','slice'))) clustering <- clustering$clustering
    clustering[clustering==from] <- to
    out <- list()
    out$clustering <- as.numeric(factor(clustering))
    out
}
