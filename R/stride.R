stride <- function(seq,arg2,type='pam',numrep=10,maxitr=100)
{
    ncol <- length(seq)
    pnt <- 0

    if (class(arg2) == 'dist') {
        if (type == 'pam') {
            res <- matrix(NA,nrow=attr(arg2,'Size'),ncol=ncol)
            for (i in seq) {
                pnt <- pnt + 1
                res[,pnt] <- pam(arg2,i)$clustering
            }
            out <- data.frame(res)
            names(out) <- as.character(seq)
            row.names(out) <- attr(arg2,'Labels')
            out <- list(clustering=out,seq=seq)
            class(out) <- 'stride'
            out
        }
        else {
            res <- matrix(NA,nrow=attr(arg2,'Size'),ncol=ncol)
            for (i in seq) {
                pnt <- pnt + 1
                res[,pnt] <- bestopt(arg2,i,numrep=numrep,maxitr=maxitr)$clustering
            }
            out <- data.frame(res)
            names(out) <- as.character(seq)
            row.names(out) <- attr(arg2,'Labels')
            out <- list(clustering=out,seq=seq)
            class(out) <- 'stride'
            out
        }        
    }
    else if (class(arg2) == 'hclust') {
        res <- matrix(NA,nrow=length(arg2$order),ncol=ncol)
        for (i in seq) {
            pnt <- pnt + 1
            res[,pnt] <- cutree(arg2,i)
        }
        out <- data.frame(res)
        names(out) <- as.character(seq)
        row.names(out) <- arg2$labels
        out <- list(clustering=out,seq=seq)
        class(out) <- 'stride'
        out
    }
    else print('you must enter a dissimilarity or a tree')
}

plot.stride <- function (x,dist,col2=4, ...) 
{
    oldpar <- par(no.readonly=TRUE)
    par(mar=c(5,4,4,4)+0.1,bg='white',fg='black')
    plot(partana(x,dist),xlab='Number of Clusters',
        ylab='Partana Ratio',type='b')
    par(new=TRUE)
    tmp <- silhouette(x,dist)
    plot(tmp,axes=FALSE,xlab='',ylab='',col=col2,type='b')
    axis(4,at=pretty(range(tmp$sil_width)),col.axis=col2,col.tick=col2)
    mtext('Silhouette Width',4,2.5,col=col2)
    par(oldpar)
}

extract <- function(stride,k)
{
    UseMethod("extract")
}

extract.stride <- function(stride,k)
{
   pnt <- which(stride$seq==k)
   out <- list(clustering=stride$clustering[,pnt])
   class(out) <- 'clustering'
   out
}
 
