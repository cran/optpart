\name{clusplot.partana}
\alias{clusplot.partana}
\alias{clusplot.clustering}
\title{Bi-variate Clusplot of a Partana Object}
\description{
   Clusplot (Clustering Plot) method for an object of class \code{partana}.
}
\usage{
\method{clusplot}{partana}(x, dist, \dots)
\method{clusplot}{clustering}(x, dist, \dots)
}
\arguments{
    \item{x}{an object of class partana or class \sQuote{clustering},
         e.g. created by the functions \code{partana}, \code{optpart} or \code{slice}.}
    \item{dist}{an object of class \code{dist}, e.g. produced by function
        \code{dist} or \code{\link[labdsv]{dsvdis}} or
        \code{\link[vegan]{vegdist}}}
    \item{\dots}{ancillary arguments to satisfy the generic function}
}

\value{
    An invisible list with components
  \item{Distances}{
    When option lines is 1 or 2 we obtain a k by k matrix (k is the number of
    clusters). The element at row j and column s is the distance between
    ellipse j and ellipse s.  If lines=0, then the value of this
    component is NA.
  }
  \item{Shading}{
    A vector of length k (where k is the number of clusters), containing
    the amount of shading per cluster.  Let y be a vector where element
    i is the ratio between the number of objects in cluster i and the
    area of ellipse i.  When the cluster i is a line segment, y[i] and
    the density of the cluster are set to NA.  Let z be the sum of all
    the elements of y without the NAs.  Then we put shading = y/z *37 + 3.
  }
}

\details{
  This code is simply a conversion routine to pass the appropriate values of a
  partana or clustering object to the \code{\link[cluster]{clusplot}} function.  Clusplots
  are two-dimensional representations of cluster structure superimposed on an
  eigenanalysis of a dissimilarity/distance matrix.  This function extracts components
  from the partana object and calls \code{partana.default} which performs
  all calculations and does the plotting.

}

\note{\code{\link[cluster]{clusplot}} is a graphic function from package \sQuote{cluster}
to graphically represent clusters in a reduced dimensional space.  \sQuote{clusplot.partana}
is simply a wrapper function to employ that function on objects of class \code{partana}}
\author{
	David W. Roberts
	\email{droberts@montana.edu}
	\url{http://ecology.msu.montana.edu/labdsv/R/optpart}
}
\seealso{
\code{\link[cluster]{clusplot}}
}
\examples{
    \dontrun{library(labsv)}
    \dontrun{data(shoshveg)}
    \dontrun{dis.bc <- dsvdis(shoshveg,'bray/curtis')}
    \dontrun{opt.5 <- optpart(5,dis.bc)}
    \dontrun{clusplot(opt.5)}
}
\keyword{cluster}
\keyword{hplot}

