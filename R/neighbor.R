neighbor <- function (x,all=FALSE) 
{
    if (inherits(x,'pam')) {
        tmp <- silhouette(x)
        if (!all) tmp <- tmp[tmp[,3]<0,]
        out <- table(tmp[,1],tmp[,2])
    }
    else if (inherits(x,'partana')) {
        if (all) {
            y <- x$clustering
            x <- x$ptc
            for (i in 1:nrow(x)) {
                x[i,y[i]] <- NA
            }
            out <- table(y,apply(x,1,which.max))
        }
        else {
            tmp <- testpart(x)
            out <- table(tmp[,1],tmp[,2])
        }
    }
    else stop('The first argument must be of class pam or optpart')
    out
}
