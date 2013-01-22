\name{testopt}
\alias{testopt}
\title{Identify Misclassified Plots in an Optpart Classification}
\description{Optimal partitioning is an iterative re-allocation algorithm to maximize the ratio of within-cluster similarity/among-cluster similarity for a given number of clusters.  testopt analyzes the mean similarity of each sample to the cluster to which it is assigned to all other clusters, and lists those samples which have similarity higher to another cluster than to the one to which they are assigned.
}
\usage{testopt(part,ord=TRUE,digits=4,format='f')
}
\arguments{
    \item{part}{a object of class \sQuote{partana} from \sQuote{optpart}
 or \sQuote{partana}}
    \item{ord}{a switch to control whether the output is ordered}
    \item{digits}{how many digits to print}
    \item{format}{number of digits or number of significant digits}
}
\details{
Simply examines each sample plot, comparing the mean similarity of that sample to all other samples in the cluster to which it is assigned as compared to its mean similarity to all other clusters. Samples which are more similar to other clusters than to the one to which they are assigned are listed in a table which gives their current cluster assignment, the cluster to which they are more similar, and the mean similarities of that sample to all clusters.

If ord=TRUE then the output is ordered to reflect target
clusters.

The format argument controls how the argument digits is
interpreted.  See \code{formatC} for more information.
}

\value{a table of values
}
\author{
    David W. Roberts
    \email{droberts@montana.edu}
}
\seealso{
\code{partana}, \code{silhouette}
}
\examples{
    data(shoshveg)
    dis.bc <- dsvdis(shoshveg,'bray/curtis')
    opt.5 <- optpart(5,dis.bc)
    testopt(opt.5)
}
\keyword{cluster}
