disdiam <- function(x,dist,digits)
{
    UseMethod("disdiam")
}

disdiam.default <- function (x, dist, digits = 3)
{
    clustering <- x
    if (is.numeric(clustering)) {
        if (min(clustering)< 0 || (length(table(clustering)) != max(clustering))) {
            cat('WARNING: renumbering clusters to consecutive integers\n')
            clustering <- match(clustering,sort(unique(clustering)))
        }
    }

    if (class(dist) != "dist")
        stop("The second argument must an object of class dist")
    dist <- as.matrix(dist)
    nclustering <- as.numeric(clustering)
    cname <- names(table(clustering))
    csize <- as.numeric(table(clustering))
    diam <- rep(0, length(table(clustering)))
    for (i in 1:length(table(clustering))) {
        pnt <- nclustering == i
        subdis <- dist[pnt, pnt]
        diam[i] <- max(subdis)
    }
    diam <- as.numeric(format(diam, digits = digits, nsmall = digits))
    res <- data.frame(cname, csize, diam)
    names(res) <- c("cluster", "N", "diameter")
    mean <- sum(res$N[res$N > 1] * res$diameter[res$N > 1])/sum(res$N[res$N > 1])
    out <- list(diameters = res, mean = mean)
    class(out) <- "disdiam"
    out
}

disdiam.clustering <- function (x,dist,digits=3) 
{
    clustering <- x$clustering
    if (min(clustering)< 0 || (length(table(clustering)) != max(clustering))) {
        cat('WARNING: renumbering clusters to consecutive integers\n')
        clustering <- match(clustering,sort(unique(clustering)))
    }

    if (class(dist) != 'dist') stop('The second argument must an object of class dist') 
    dist <- as.matrix(dist)
    nclustering <- as.numeric(clustering)
    cname <- names(table(clustering))
    csize <- as.numeric(table(clustering))
    diam <- rep(0,length(table(clustering)))

    for (i in 1:length(table(clustering))) {
        pnt <- nclustering == i
        subdis <- dist[pnt,pnt]
        diam[i] <- max(subdis)
    }
    diam <- as.numeric(format(diam,digits=digits,nsmall=digits))
    res <- data.frame(cname,csize,diam)
    names(res) <- c('cluster','N','diameter')
    mean <- sum(res$N[res$N>1]*res$diameter[res$N>1])/sum(res$N[res$N > 1])
    out <- list(diameters=res,mean=mean)
    class(out) <- 'disdiam'
    out
}

disdiam.partition <- function(x,dist,digits=3)
{    
    clustering <- x$clustering
    if (class(dist) != 'dist') stop('The second argument must an object of class dist')
    dist <- as.matrix(dist)
    nclustering <- as.numeric(clustering)
    cname <- names(table(clustering))
    csize <- as.numeric(table(clustering))
    diam <- rep(0,length(table(clustering)))

    for (i in 1:length(table(clustering))) {
        pnt <- nclustering == i
        subdis <- dist[pnt,pnt]
        diam[i] <- max(subdis)
    }
    diam <- as.numeric(format(diam,digits=digits,nsmall=digits))
    res <- data.frame(cname,csize,diam)
    names(res) <- c('cluster','N','diameter')
    mean <- sum(res$N[res$N>1]*res$diameter[res$N>1])/sum(res$N[res$N > 1])
    out <- list(diameters=res,mean=mean)
    class(out) <- 'disdiam'
    out
}

disdiam.stride <- function (x, dist, digits = 3)
{
    res <- rep(NA, ncol(x$clustering))
    dist <- as.matrix(dist)
    for (i in 1:ncol(x$clustering)) {
        members <- table(x$clustering[,i])
        sum <- 0
        for (j in 1:x$seq[i]) {
            if (members[j] > 1) {
                pnt <- x$clustering[,i] == j
                subdis <- dist[pnt,pnt]
                diam <- max(subdis)
                sum <- sum + diam * members[j]
            }
        }
    res[i] <- round(sum/sum(members[members>1]),digits)
    }
    clusters <- x$seq
    diameters <- res
    out <- data.frame(clusters, diameters)
    out
}

print.disdiam <- function(x, ...)
{
    print(x$diameters)
    cat(paste('\nMean = ',format(x$mean,digits=4),"\n"))
}

