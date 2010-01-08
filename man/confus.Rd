\name{confus}
\alias{confus}
\title{Confusion Matrix}
\description{A confusion matrix is a cross-tabulation of actual class
membership with memberships predicted by a discriminant function,
classification tree, or other predictive model.}
\usage{confus(clustering,fitted)}
\arguments{
    \item{clustering}{an object of class \sQuote{clustering} or
a vector of (integer) class membership values}
    \item{fitted}{a matrix of probabilities of membership in classes, such
           as produced by \code{fitted} functions}
}
\details{Cross-classifies each sample by actual class membership and
predicted membership, computing overall accuracy, and the Kappa
statistic of agreement.}

\value{produces a list with elements
\item{matrix}{the cross-tabulation matrix as a data.frame}
\item{correct}{the fraction of correctly predicted samples}
\item{kappa}{the value of the Kappa statistic}
\item{legend}{the text legend for the cross-tabulation matrix}
}
\references{\url{http://ecology.msu.montana.edu/labdsv/R}}
\note{Confusion matrices are commonly computed in remote sensing
applications, but are equally suited to the evaluation of any
predictive methods of class membership or factors.}
\author{
    David W. Roberts
    \email{droberts@montana.edu}
    \url{http://ecology.msu.montana.edu/labdsv/R/optpart}
}
\examples{
    data(shoshveg) # returns a data frame of vegetation data
    data(shoshsite) # returns a data frame of site data
    library(tree)
    mod <- tree(factor(shoshveg$PSEMEN>0)~ elevation+slope+av,data=shoshsite)
    pred <- predict(mod,newdata=shoshsite)
    confus(shoshveg$PSEMEN>0,pred)
}
\keyword{htest}
