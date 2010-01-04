\name{plot.partana}
\alias{plot.partana}
\title{Plotting a Partition Analysis}
\description{Plots an object of class \sQuote{partana}}

\usage{
\method{plot}{partana}(x, zlim = range(x$ptc), col = heat.colors(12), ...)
}
\arguments{
    \item{x}{an object of class \sQuote{partana}}
    \item{zlim}{the minimum and maximum values to define the color range of the
             \code{image} function}
    \item{col}{the color map to use in the plot}
    \item{\dots}{ancillary arguments to function \code{plot}}
}
\details{Plots two panels in sequence in the current device.  The first shows the mean similarity
of every object to each cluster, sorted by mean similarity to the other
members of its own cluster, with objects as columns and clusters as rows.
The second panel shows the mean similarity of every cluster to every other cluster
and mean within-cluster similarity, ignoring cluster size.  These plots are known as
\sQuote{Mondriaan} plots, where the similarities are given by lines colored from min to max.

If the \sQuote{partana} object was produced by \code{\link[optpart]{optpart}},
a third panel is plotted showing the trace of the optimization.
}

\value{two (or three) plots in sequence on the current graphic device}

\author{
    David W. Roberts
    \email{droberts@montana.edu}
}
\seealso{\code{optpart}, \code{\link[optpart]{partana}},
             \code{\link[cluster]{partition}}}
\examples{
    \dontrun{data(shoshveg)}
    \dontrun{dis.bc <- dsvdis(shoshveg,'bray/curtis')}
    \dontrun{demo.part <- partana(sample(1:5,nrow(shoshveg),replace=TRUE),dis.bc)}
    \dontrun{plot(demo.part)}
}
\keyword{hplot}
