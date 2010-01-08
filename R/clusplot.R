clusplot.partana <- function (x,dist,...)
{
    cluster:::clusplot.default(dist,x$clustering,diss=TRUE)
}

clusplot.clustering <- function (x,dist,...)
{
    cluster:::clusplot.default(dist,x$clustering,diss=TRUE)
}
