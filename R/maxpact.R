maxpact <- function (dist,size,alphac=0.01,mean=FALSE)
{
    if (class(dist) != 'dist') stop('The first argument must be an object of class dist')
    sim <- 1 - as.matrix(dist)
    numplt <- nrow(sim)
    if (mean) morf <- 1
    else morf <- 0
    musuby <- matrix(0,nrow=numplt,ncol=numplt)
    membry <- matrix(0,nrow=numplt,ncol=numplt)
    numset <- 0
    used <- rep(0,numplt)
    musubx <- matrix(0,nrow=numplt,ncol=numplt)
    membrx <- matrix(0,nrow=numplt,ncol=numplt)
    mnsimi <- rep(0,numplt)

    tmp <- .Fortran("maxpact",
                as.double(sim),
                as.integer(numplt),
                as.integer(size),
                as.double(alphac),
                as.integer(morf),
                musuby=as.double(musuby),
                membry=as.integer(membry),
                numset=as.integer(numset),
                as.integer(used),
                as.double(musubx),
                as.integer(membrx),
                as.double(mnsimi),
                PACKAGE='optpart')
    member <- matrix(tmp$membry,nrow=numplt)
    member <- member[1:tmp$numset,1:size]
    musuby <- matrix(tmp$musuby,nrow=numplt)
    musubx <- musuby[1:tmp$numset,1:size]
    distname <- deparse(substitute(dist))
    out <- list(musubx=musubx,member=member,numset=tmp$numset, size=size, distname=distname)
    class(out) <- "mps"
    return(out)
}

mps.test <- function(mps,env,main=deparse(substitute(env)),...)
{
    if (class(mps) != "mps") {
        stop("You must pass an object of class mps from maxpact()")
    }
    if (!is.numeric(env)) {
        stop("You must pass only numeric vectors as environment variables")
    }
    size <- ncol(mps$member)
    nset <- mps$numset
    null <- rep(0,nset)
    for (i in 1:nset) {
        tmp <- sample(1:length(env),size)
        nullmin <- min(env[tmp])
        nullmax <- max(env[tmp])
        null[i] <- nullmax - nullmin
    }
        
    low <- apply(mps$member,1,function(x){min(env[x])})
    high <- apply(mps$member,1,function(x){max(env[x])})
    observed <- high - low
    plot(sort(null),ylim=c(0,max(null)),main=main,ylab="Within-Set Difference")
    points(sort(observed),col=2)
    readline("Hit return\n")
    boxplot(null,observed,names=c("null","observed"),ylab="Within-Set Difference",main=main)
    print(wilcox.test(null,observed))
}

plot.mps <- function(x,...)
{
    plot(x$musubx[1,],ylim=c(0,1),xlab="Size",ylab="Similarity",type="n")
    for (i in 1:x$numset) {
        lines(x$musubx[i,])
    }
}
