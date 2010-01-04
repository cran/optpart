\name{isamic}
\alias{isamic}

\title{Indicator Species Analysis to Minimize Intermediate Constancy}
\description{Calculates the indicator value (extreme constancies across types)
of species in clusters or types}

\usage{
isamic(taxa,clustering,sort=FALSE)
}

\arguments{
    \item{taxa}{a matrix or data.frame of samples with species as columns and
               samples as rows}
    \item{clustering}{a vector of numeric cluster memberships for samples, or a
                 classification object returned from \code{\link[cluster]{pam}}, 
                 or \code{\link[optpart]{partana}}}
    \item{sort}{a switch to control sorting the results table}
}
\details{Calculates the indicator value \eqn{d} of species based on constancy within
types or clusters.  Specifically,

where:\cr
\eqn{c_{i,j} = } constancy of species \eqn{i} in 
   cluster \eqn{j} of \eqn{k} clusters; \cr

\deqn{d_i = 2 * \sum_j abs(c_{i,j}-0.5) / k}
}

\value{a vector of indictor values}

\note{Indicator value analysis is a set of techniques designed to identify
species of special interest in clusters or types.  The most widely used indicator species
analysis was proposed by Dufrene and Legendre (1997), and is included in package
\sQuote{labdsv} as \code{\link[labdsv]{indval}}. isamic differs significantly from
indval in assumption and objective, seeking to identify species that largely present or 
absent with types, and down-weighting species with intermediate constancies.}

\author{
    David W. Roberts
    \email{droberts@montana.edu}
    \url{http://ecology.msu.montana.edu/labdsv}
}

\seealso{\code{\link[labdsv]{indval}}, \code{\link[optpart]{tabdev}}}

\examples{
    data(shoshveg) # returns a vegetation dataframe
    dis.bc <- dsvdis(shoshveg,'bray/curtis') # returns a dissimilarity matrix
    clust <- optpart(5,dis.bc)
    isamic(shoshveg,clust)
}
\keyword{cluster}
