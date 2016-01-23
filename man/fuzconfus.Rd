\name{fuzconfus}
\alias{fuzconfus}
\title{Fuzzy Confusion Matrix}
\description{A confusion matrix is a cross-tabulation of actual class
membership with memberships predicted by a discriminant function,
classification tree, or other predictive model.  A fuzzy confusion
matrix is a confusion matrix that corrects for \sQuote{near misses} 
in prediction by comparing the similarity of the predicted type to
the actual type and giving credit for the similarity.}
\usage{fuzconfus(part,fitted,dis)}
\arguments{
    \item{part}{an object of class \sQuote{partana}, \sQuote{clustering}, or \sQuote{pam}
a vector of (integer) class membership values}
    \item{fitted}{a matrix of probabilities of membership in classes, such
           as produced by \code{fitted} functions}
    \item{dis}{an object of class \sQuote{dist} such as produced by 
             \sQuote{dist} \code{\link[labdsv]{dsvdis}}, or \code{\link[vegan]{vegdist}}.
             if \sQuote{part} is an object of class \sQuote{clustering} or \sQuote{'pam'}
             argument \sQuote{dis} must be provided.}
}
\details{Cross-classifies each sample by actual class membership and
predicted membership,  and then calculates the similarity of the predicted type to the actual.
Correct predictions are assigned values of 1.0, and other predictions are given the value of the 
similarity of the two types on the diagonal.  The dissimilarity of the two types is added off 
the diagonal as fuzzy error.
}

\value{produces a list with elements
\item{confuz}{the fuzzy cross-tabulation matrix as a data.frame}
\item{correct}{the fraction of fuzzy correctly predicted samples}
\item{percent}{the percent of fuzzily correct predictions}
\item{kappa}{the value of the Kappa statistic}
}
\references{\url{http://ecology.msu.montana.edu/labdsv/R}}
\note{Confusion matrices are commonly computed in remote sensing
applications, but are equally suited to the evaluation of any
predictive methods of class membership or factors.}
\author{
    David W. Roberts
    \email{droberts@montana.edu}
    \url{http://ecology.msu.montana.edu/labdsv/R}
}
\examples{
    data(shoshveg) # returns a data frame of vegetation data
    data(shoshsite) # returns a data frame of site data
    dis.bc <- dsvdis(shoshveg,'bray/curtis')
    opt.5 <- optpart(5,dis.bc)
    library(tree)
    mod <- tree(factor(opt.5$clustering)~ elevation+slope+av,data=shoshsite)
    pred <- predict(mod,newdata=shoshsite)
    fuzconfus(opt.5,pred,dis.bc)
}
\keyword{htest}
