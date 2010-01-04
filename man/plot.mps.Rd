\name{plot.mps}
\alias{plot.mps}
\title{Plotting Maximally Compact Sets}
\description{Produces an analytical plot of a maximally compact sets
The \code{plot} function plots the membership values for each maximally compact set in the order
samples are added to the set.  The initial slope is indicative of the connectedness of the
similarity relation, as it represent the range of \sQuote{nearest neighbor} distances for
the initial seeds.  The subsequent slope is indicative of the range of values  or \sQuote{local density}
of the similarity relation.
}
\usage{
\method{plot}{mps}(x, \dots)
}
\arguments{
    \item{x}{an object of class \sQuote{mps}}
    \item{\dots}{ancillary arguments to pass to \sQuote{plot}}
}
\value{
Produces a plot on the current graphic device.
}
\author{
	David W. Roberts
	\email{droberts@montana.edu}
}
\seealso{
\code{\link[optpart]{maxpact}},\code{\link[optpart]{mps.test}},
}
\examples{
    data(shoshveg)
    dis.bc <- dsvdis(shoshveg,'bray/curtis')
    mps.10 <- maxpact(dis.bc,10)
    \dontrun{plot(mps.10)}
}
\keyword{cluster}

