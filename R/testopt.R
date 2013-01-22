testopt <- function (part, ord=TRUE, digits=4, format='f') 
{
    if (!inherits(part,'partana'))
        stop("you must pass an object of class partana")
    clustering <- part$clustering
    ptc <- part$ptc
    best <- apply(ptc, 1, which.max)
    out <- NULL
    for (i in 1:nrow(ptc)) {
        if (clustering[i] != best[i]) {
            out <- rbind(out, c(i, part$names[i], clustering[i], 
                best[i], formatC(ptc[i, ], digits = digits, format=format)))
        }
    }
    out <- data.frame(out, row.names = 1)
    if (ord) out <- out[order(out[, 2], out[, 1]), ]
    names(out) <- c("sample","cluster", "better", as.character(1:ncol(ptc)))
    out
}

