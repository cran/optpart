\name{partition}
\alias{partition}
\alias{partition.partana}
\alias{partition.clustering}
\title{Convert Object to Partition Object}
\description{Convert an object of class \sQuote{partana} or class \sQuote{clustering}
to an object of class \sQuote{partition}.}

\usage{
\method{partition}{partana}(x, dist, \dots)
\method{partition}{clustering}(x, dist, \dots)
}

\arguments{
    \item{x}{an object of class \sQuote{partana} or \sQuote{clustering}}
    \item{dist}{an object of class \sQuote{dist}}
    \item{\dots}{ancillary arguments to pass to \sQuote{partition}}
}
\details{A \sQuote{partition} object is the output of several functions
in package \sQuote{cluster}.  This utility function converts objects
from package \sQuote{optpart} to \sQuote{partitions} so that
functions in that library are available.
}

\value{an object of class \sQuote{partition} with components (and possibly others):
\item{clustering}{vector of numeric cluster assignments}
\item{silinfo}{a list with all silhouette information, only available when
          the number of clusters is non-trivial, i.e., \eqn{1 < k < n}.
          See \sQuote{silhouette}}
}
\references{\url{http://ecology.montana.msu.edu/optpart/}}
\author{
    David W. Roberts
    \email{droberts@montana.edu}
    \url{http://ecology.msu.montana.edu/labdsv/R/optpart}
}
\seealso{\code{\link{silhouette}}, \code{\link{partition}},
\code{\link{optpart}}}
\examples{
    data(shoshveg)
    dis.bc <- dsvdis(shoshveg,'bray/curtis')
    \dontrun{part <- partana(sample(1:5,nrow(shoshveg),replace=TRUE),dis.bc)}
    \dontrun{result <- partition(part,dis.bc)}
}
\keyword{cluster}
