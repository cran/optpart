testpart <- function (part, ord = TRUE)
{
    if (!inherits(part, "partana"))
        stop("you must pass an object of class partana")
    clustering <- part$clustering
    ptc <- part$ptc
    best <- apply(ptc, 1, which.max)
    worst <- part$clustering != best

    assigned <- part$clustering[worst]
    better <- best[worst]
    old <- rep(0,sum(worst))
    new <- rep(0,sum(worst))

    ptc <- part$ptc[worst,]

    for (i in 1:nrow(ptc)) {
        old[i] <- ptc[i,assigned[i]]
        new[i] <- ptc[i,better[i]]
    }

    tmp <- data.frame(assigned,better,old,new)
    row.names(tmp) <- part$names[worst]

    if (ord) {
        tmp <- tmp[order(tmp$assigned,tmp$better),]
    }

    tmp
}

