bestfit <- function(x,cluster) {

    if (inherits(x,'partana')) {
       clustid <- x$clustering
       rows <- x$ptc[clustid==cluster,]
       vals <- rows[,cluster]
       names <- x$names[clustid==cluster]
    }
    if (inherits(x,'silhouette')) {
       rows <- x[,1] == cluster
       vals <- as.numeric(x[rows,3])
       names <- dimnames(x)[[1]]
       names <- names[rows]
    }
    fit <- data.frame(names[rev(order(vals))],sort(vals,decreasing=TRUE))
    names(fit) <- c('ID','fit')
    fit
}

