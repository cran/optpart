\name{testopt}
\alias{testopt}
\title{Identify Misclassified Plots in an Optpart Classification}
\description{Optimal partitioning is an
iterative re-allocation algorithm to maximize the ratio of within-cluster 
similarity/among-cluster similarity for a given number of clusters.  testopt analyzes 
the mean similarity of each sample to the cluster to which it is assigned to all 
other clusters, and lists those samples which have similarity higher to another
cluster than to the one to which they are assigned.}
\usage{testopt(part,ord=TRUE)}
\arguments{
	\item{part}{a object of class \sQuote{partana} from \code{optpart} or  \code{partana}} 
        \item{ord}{a switch to control whether the output is ordered}
}
\details{Simply examines each sample plot, comparing the mean similarity of that 
sample to all other samples in the cluster to which it is assigned as compared to its
mean similarity to all other clusters.  Samples which are more similar to other
clusters than to the one to which they are assigned are listed in a table which gives
their current cluster assignment, the cluster to which they are more similar, and 
the mean similarities of that sample to all clusters.

If \sQuote{ord=TRUE} then the output is ordered t reflect target clusters.
}

\value{a table of values}
\author{
	David W. Roberts
	\email{droberts@montana.edu}
}
\seealso{
\code{partana}, \code{partana}, \code{\link[cluster]{silhouette}}
}
\examples{
    data(shoshveg)
    dis.bc <- dsvdis(shoshveg,'bray/curtis')
    opt.5 <- optpart(5,dis.bc)
    testopt(opt.5)
}
\keyword{cluster}

