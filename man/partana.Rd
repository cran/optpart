\name{partana}
\alias{partana}
\alias{partana.default}
\alias{partana.partition}
\alias{partana.clustering}
\alias{partana.stride}
\alias{summary.partana}
\alias{testpart}
\title{Partition Analysis}
\description{Partition analysis evaluates the within-cluster to among-cluster
similarity of classifications as a measure of cluster validity}

\usage{
\method{partana}{default}(c, dist)
\method{partana}{partition}(c, dist)
\method{partana}{clustering}(c, dist)
\method{summary}{partana}(object, \dots)
testpart(part,ord=TRUE)
}

\arguments{
    \item{c}{an integer or factor vector, or an object of class \sQuote{clustering},
            \sQuote{partana}, \sQuote{partition}, \sQuote{stride}}
    \item{dist}{an object of class \sQuote{dist}.  \sQuote{dist} is always necessary
               unless you pass \sQuote{partana} an objects of class \sQuote{partition}
                which contains an embedded \sQuote{dist} object.}
    \item{object}{an object of class \sQuote{partana}}
    \item{part}{an object of class \sQuote{partana}}
    \item{ord}{a logical switch to control ordering of output}
    \item{\dots}{ancillary arguments to pass to summary}
}
\details{Calculates mean object-to-cluster similarity, mean cluster-to-cluster
similarity, and mean within-cluster to among-cluster similarity.  partana operates
on partitions or clustering produced by a wide range of algorithms, including specific 
methods for the products of functions \code{\link[optpart]{optpart}}, \code{\link[optpart]{slice}}, 
\code{\link[cluster]{pam}} and \code{\link[cluster]{diana}}.

\code{summary} produces a matrix of the mean cluster-to-cluster similarities,
and the overall within-cluster/among-cluster similarity ratio.

\code{testpart} identifies \sQuote{misfits} in a partition, defined as items with
lower dissimilarity to another cluster than to the one to which they belong.  
The \sQuote{ord} switch re-orders the output in order of poor fit.
}

\value{an object of class \sQuote{partana} with components:
\item{ptc}{matrix of mean object-to-cluster similarity}
\item{ctc}{matrix of mean cluster-to-cluster similarity}
\item{clustering}{vector of numeric cluster assignments}
\item{ratio}{within-cluster to among-cluster similarity ratio}
}
\author{
    David W. Roberts
    \email{droberts@montana.edu}
}
\seealso{\code{\link[cluster]{partition}}, \code{optpart},
         \code{plot.partana}}
\examples{
    data(shoshveg)
    data(shoshsite)
    dis.bc <- dsvdis(shoshveg,'bray/curtis')
    demo.part <- partana(cut(shoshsite$elev,5),dis.bc)
    summary(demo.part)
}
\keyword{cluster}
