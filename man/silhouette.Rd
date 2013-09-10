\name{silhouette.partana}
\alias{silhouette.partana}
\alias{silhouette.clustering}
\alias{testsil}
\title{Produce a Silhouette Object From a Partana or Clustering Object}
\description{Extracts components from a \code{partana}
or \code{clustering} object,
and passes the values to the \code{\link[cluster]{silhouette}} function to
produce an object of class silhouette.
}
\usage{
\method{silhouette}{partana}(x, dist, \dots)
\method{silhouette}{clustering}(x, dist, \dots)
testsil(sil)
}
\arguments{
    \item{x}{an object of class partana}
    \item{dist}{an object of class dist}
    \item{\dots}{miscellaneous arguments to pass to function silhouette}
    \item{sil}{an object of class \sQuote{silhouette'}}
}
\details{Extracts the component \sQuote{clustering} from the partana or clustering object,
and calls function \code{\link[cluster]{silhouette}} in library \sQuote{cluster}.

\code{testsil} identifies \sQuote{misfits} in a partition, defined as plots with
a negative silhouette width, and prints them out in a sorted list.
}

\value{
An object of class \code{\link[cluster]{silhouette}}
}
\references{\url{http://ecology.montana.msu.edu/optpart}}
\note{This is a a simple conversion routine to allow plotting a silhouette plot
for an object of class partana.}
\author{
    David W. Roberts
    \email{droberts@montana.edu}
}
\seealso{\code{\link[cluster]{silhouette}}}

\examples{
    data(shoshveg) # produces a data frame of vegetation data, 
                   #samples as rows, attributes as columns
    dis.bc <- dsvdis(shoshveg,'bray/curtis') # produces a Bray/Curtis dissimilarity matrix
    opt.5 <- optpart(5,dis.bc) # produces an optimal partitioning into 5 clusters
    silhouette(opt.5,dis.bc) # calculates the silhouette values
    \dontrun{plot(silhouette(opt.5,dis.bc))} # produce silhouette plot on current device
}
\keyword{cluster}
