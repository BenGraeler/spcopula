\name{composeSpCopula}
\alias{composeSpCopula}
\title{
Composing a bivariate Spatial Copula
}
\description{
A bivariate spatial copula is composed out of a set of bivariate copulas. These are combined using a convex linear combination with weights based on distances where for copulas with a 1-1 correspondence of Kendall's tau or Spearman's rho a dependence function providing measures of association based on distances might be used.
}
\usage{
composeSpCopula(bestFit, families, bins, calcCor, range=max(bins$meanDists))
}
\arguments{
  \item{bestFit}{
A vector indicating the best fitting copula family index (as in \code{families}) per lag-class. If this vector is shorter than the actual number of bins, only these first bins will be used. It is an alternative way of setting a range.
}
  \item{families}{
The list of investigated families.
}
  \item{bins}{
The spatial bins used.
}
  \item{calcCor}{
Optional, some spatially dependent function providing values of Kendall's tau for distances.
}
  \item{range}{
  The maximal distance to be used for the spatial copula. By default, the maximum distance provided by bins will be used. 
  }
}
\value{
A (\code{\linkS4class{spCopula}}) object.
}
\author{
Benedikt Graeler
}

\examples{
composeSpCopula(c(1,1,2,3),families=list(frankCopula(.4), gumbelCopula(1.6),gumbelCopula(1.4)),
                bins=data.frame(meanDists=c(500,1000,1500,2000,2500)),range=2250)
}

\keyword{ spatial }
\keyword{ multivariate }
\keyword{ distribution }
