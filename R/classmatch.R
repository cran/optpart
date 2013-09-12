classmatch <- function (x,y,type='full') 
{
    if (inherits(x,c('partana','clustering','partition'))) x <- x$clustering
    if (min(x)< 0 || (length(table(x)) != max(x))) {
        cat('WARNING: renumbering clusters to consecutive integers\n')
        x <- match(x,sort(unique(x)))
    }

    if (inherits(y,c('partana','clustering','partition'))) y <- y$clustering 
    if (min(y)< 0 || (length(table(y)) != max(y))) {
        cat('WARNING: renumbering clusters to consecutive integers\n')
        y <- match(y,sort(unique(y)))
    }

    tab <- table(x,y)
    total <- sum(tab)

    if (type == 'full') {
        match <- sum(tab>0)
        pairs <- matrix(0,nrow=match,ncol=3)
        partial <- rep(0,match)
        combo <- rep(0,length(x))
        ord <- matrix(0,nrow=nrow(tab),ncol=ncol(tab))
        running <- 0
    
        for (i in 1:match) {
            find <- max(tab)
            for (j in 1:nrow(tab)) {
                test <- max(tab[j,])
                if (test == find) {
                     col <- which.max(as.vector(tab[j,]))
                     pairs[i,] <- c(j,col,tab[j,col])
                     tab[j,col] <- 0
                     ord[j,col] <- i 
                     break
                }
            }
        }
    
        for (i in 1:length(x)) {
            for (j in 1:nrow(pairs)) {
                if (x[i] == pairs[j,1] && y[i] == pairs[j,2]) {
                    combo[i] <- j
                    break
                }
            }
        }
        partial <- cumsum(pairs[,3])/total
        pairs <- data.frame(pairs)
        names(pairs) <- c('row','column','n')
        res <- list(tab=table(x,y),pairs=pairs,partial=partial,ord=ord,combo=combo)
        return(res)
    }
    else {
        fulltab <- tab
        match <- rep(0,nrow(tab))
        ord <- matrix(0,nrow=nrow(tab),ncol=ncol(tab))
        pairs <- matrix(0,nrow=nrow(tab),ncol=3)
    
        for (i in 1:nrow(tab)) {
            cols <- max.col(tab)
            vals <- rep(NA,ncol(tab))
            for (j in 1:nrow(tab)){
                vals[j] <- tab[j,cols[j]]
            }
            row <- which.max(vals)
            col <- which.max(as.vector(tab[row,]))
            match[i] <- tab[row,col]
            pairs[i,] <- c(row,col,match[i])
            ord[row,col] <- i
            tab[row,] <- 0
            tab[,col] <- 0
        }
	res <- list(tab=fulltab,pairs=pairs,partial=cumsum(match)/total,ord=ord)
        return(res)
    }
}


