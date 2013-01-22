testpart <- function (part, ord = TRUE, minval = 0, digits=4) 
{
    if (!inherits(part, "partana")) 
        stop("you must pass an object of class partana")
    clustering <- part$clustering
    ptc <- part$ptc
    best <- apply(ptc, 1, which.max)
    worst <- clustering != best
    assigned <- clustering[worst]
    better <- best[worst]
    old <- rep(0, sum(worst))
    new <- rep(0, sum(worst))
    ptc <- part$ptc[worst, ]
    for (i in 1:nrow(ptc)) {
        old[i] <- round(ptc[i, assigned[i]],digits=digits)
        new[i] <- round(ptc[i, better[i]],digits=digits)
    }
    tmp <- data.frame(assigned, better, new)
    row.names(tmp) <- part$names[worst]
    if (ord) {
        tmp <- tmp[order(tmp$assigned, tmp$better), ]
    }
    tmp <- tmp[tmp$new - tmp$old > minval, ]
    badtypes <- rep(0, length(table(clustering)))
    for (i in clustering) {
        badtypes[i] <- round(sum(tmp[, 1] == i)/sum(clustering == i),digits=digits)
    }
    badtypes <- data.frame(badtypes)
    names(badtypes) <- "fraction of misfits"
    out <- list(plots = tmp, types = badtypes)
    class(out) <- "testpart"
    out
}

print.testpart <- function (x, ...) 
{
    print(x$plot)
    cat("\n")
    print(x$types)
}

