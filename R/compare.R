compare <- function (const,left,right,thresh=0.2) 
{
#   if (!inherits(const,'const')) 
#       stop("The first argument must be an object of class 'const'")
#   if (!is.integer(left) || !is.integer(right))
#       stop('The second and third arguments must be integers')
    new <- const[,c(left,right)]
    either <- as.numeric(new[,1])>thresh | as.numeric(new[,2])>thresh
    new <- new[either,]
    diff <- as.numeric(new[,1])-as.numeric(new[,2])
    new$diff <- diff
    a <- new[diff > thresh,]
    a <- a[rev(order(a$diff)),]
    b <- new[diff < (-1*thresh),]
    b <- b[order(b$diff),]
    b$diff <- -1 * b$diff
    c <- new[abs(diff) < thresh,]
    c <- c[order(c$diff),]
    out <- list(left=a,right=b,both=c)
    out
}
