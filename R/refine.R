refine <- function (x, clustering, ...) 
{
   UseMethod("refine")
}

refine.pco <- function (x, clustering, ax=1, ay=2, ...) 
{
    if (inherits(clustering,c('clustering','partana','partition'))) clustering <- clustering$clustering
    if (is.numeric(clustering)) {
        if (min(clustering)< 0 || (length(table(clustering)) != max(clustering))) {
            cat('WARNING: renumbering clusters to consecutive integers\n')
            clustering <- match(clustering,sort(unique(clustering)))
        }
    }

    for (i in 1:max(clustering)) {
        plot(x, ax, ay)
        cat(paste("Refining cluster # ", i, "\n"))
        hilight(x, clustering, ax, ay)
        chullord(x, clustering == i, ax, ay, col = i + 1)
        new <- plotid(x)
        clustering[new] <- i
        points(x, clustering == i, ax, ay, col = i + 1)
    }
    plot(x, ax, ay)
    hilight(x, clustering, ax, ay)
    for (i in 1:max(clustering)) {
        chullord(x, clustering == i, ax, ay, col = i + 1)
    }
    out <- list(clustering=clustering)
    attr(out, "class") <- "clustering"
    return(out)
}

refine.nmds <- function (x, clustering, ax=1, ay=2, ...) 
{
    if (inherits(clustering,c('clustering','partana','partition'))) clustering <- clustering$clustering
    if (min(clustering)< 0 || (length(table(clustering)) != max(clustering))) {
        cat('WARNING: renumbering clusters to consecutive integers\n')
        clustering <- match(clustering,sort(unique(clustering)))
    }

    for (i in 1:max(clustering)) {
        plot(x, ax, ay)
        cat(paste("Refining cluster # ", i, "\n"))
        hilight(x, clustering, ax, ay)
        chullord(x, clustering == i, ax, ay, col = i + 1)
        new <- plotid(x)
        clustering[new] <- i
        points(x, clustering == i, ax, ay, col = i + 1)
    }
    plot(x, ax, ay)
    hilight(x, clustering, ax, ay)
    for (i in 1:max(clustering)) {
        chullord(x, clustering == i, ax, ay, col = i + 1)
    }
    out <- list(clustering=clustering)
    attr(out, "class") <- "clustering"
    return(out)
}

refine.default <- function (x,clustering, ...) 
{
    if (inherits(clustering,c('partana','clustering'))) {
        clustering <- clustering$clustering
    }
    if (min(clustering)< 0 || (length(table(clustering)) != max(clustering))) {
        cat('WARNING: renumbering clusters to consecutive integers\n')
        clustering <- match(clustering,sort(unique(clustering)))
    }

    repeat {
        plots <- readline(' enter the plots    : ')
        if (plots == "") break
        new <- as.numeric(readline(' New cluster        : '))
        for (i in strsplit(plots,",")[[1]]){
            ord <- 1:nrow(x)
            y <- match(i,row.names(x))
            if (!is.na(y)) {
                clustering[y] <- new
            }
            else print('no such plot')
        }
    }
    out <- list(clustering=clustering)
    class(out) <- 'clustering'
    out 
}

