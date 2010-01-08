\name{maxpact}
\alias{maxpact}
\alias{mps.test}
\title{Maximally Compact Sets Analysis}
\description{Maximally compact sets is an approach to deriving relatively homogeneous
subsets of objects as determined by similarity of the composition of the objects.  The sets
so derived can be tested against random sets of the same size to determine whether a
vector of independent data exhibits a restricted distribution within the sets.
}
\usage{
maxpact(dist,size,alphac=0.01,mean=FALSE)
mps.test(mps, env, main = deparse(substitute(env)), \dots)
}
\arguments{
    \item{dist}{a dist object from \code{dist}, \code{dsvdis}, or \code{vegdist}}
    \item{size}{the size of desired compact sets}
    \item{alphac}{the alpha-cut to specify maximum dissimilarity for inclusion in a set}
    \item{mean}{if mean is FALSE (the default), the algorithm uses a furthest neighbor
	criterion; if mean is TRUE, it uses a mean similarity criterion}
    \item{mps}{an object of class \sQuote{mps}}
    \item{env}{a quantitative environmental variable for analysis}
    \item{main}{a title for the plot of mps.test}
    \item{\dots}{ancillary arguments for \sQuote{plot}}
}
\details{maxpact starts with each sample as a seed, and adds the most similar
plot to the set.  Plots are added in turn to the set (up to the size specified)
in order of maximum similarity.  If mean is FALSE, the sample most similar to
set is the sample with the max-min similarity, that is, the sample whose minimum
similarity to the set if highest, equivalent to furthest-neighbor or
complete-linkage in cluster analysis.  If mean is TRUE, the sample most similar
to a set is the sample with highest mean similarity to the set.  Once the sets
are determined for each seed, the list is examined for duplicate sets, which are
deleted, to return the list of unique sets.

If \sQuote{alphac} is specified, sets are grown to maximum size, or to maximum
dissimilarity as specified by alphac, whichever is smaller.

The \code{\link[optpart]{mps.test}} function analyzes within-set variability in attributes
of the objects other than those used to calculate the similarity relation.  If maximally compact
sets exhibit a narrower range of values than expected at random it may be that the
variable analyzed has an underlying role in determining the attributes on which the similarity
is calculated.  The function \sQuote{plot} plots the sorted within-set range of values in red,
and the sorted
range of values of random sets of the same size in black.  This followed by a boxplot of
within-set values for the random replicates versus the observed sets, and calculates a
Wilcoxon rank sum test of the difference.
}

\value{an object of class \sQuote{mps}, a list with elements:
\item{musubx}{a matrix of sample membership in the sets where membership is given by the
similarity with which a sample joined the set}
\item{member}{a list of set members in the order they were added to the set}
\item{numset}{the number of unique sets derived}
\item{size}{the number of members in each set}
\item{distname}{the name of the dissimilarity/distance object employed}
}
\author{
	David W. Roberts
	\email{droberts@montana.edu}
}
\seealso{
\code{\link[optpart]{plot.mps}}
}
\examples{
    data(shoshveg)
    library(labdsv) # to get the dsvdis() function
    dis.bc <- dsvdis(shoshveg,'bray/curtis')
    mps.10 <- maxpact(dis.bc,10)
    data(shoshsite)
    \dontrun{mps.test(mps.10,shoshsite$elevation)} # plots graph and produces summary
}
\keyword{cluster}

