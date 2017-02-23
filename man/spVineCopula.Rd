\name{spVineCopula}
\alias{spVineCopula}

\title{
Constructor of the virtual super class \code{\linkS4class{spVineCopula}}.
}
\description{
Construct an instance of the \code{\linkS4class{mixedSpVineCopula}} or \code{\linkS4class{pureSpVineCopula}} grouped as the \code{\linkS4class{spVineCopula}} super class.
}
\usage{
spVineCopula(spCop, topCop=NULL)
}
\arguments{
  \item{spCop}{Object of class \code{\linkS4class{spCopula}}. A bivariate spatially varying copula. }
  \item{topCop}{Object of class \code{\linkS4class{copula}}. A bivariate, multivariate or R-Vine copula joining the spatial components (neighbours). In case all trees are spatial copulas, this argument might be \code{NULL}.}
}
\value{
An instance of the \code{\linkS4class{spVineCopula}} class (either \code{\linkS4class{mixedSpVineCopula}} or \code{\linkS4class{pureSpVineCopula}}).
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}
\examples{
# a spatial C-vine copula (with independent dummy copulas in the upper vine)
library("VineCopula")
data("spCopDemo")

calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

library("VineCopula")
RVM <- RVineMatrix(matrix(c(1,0,0,2,2,0,3,3,3),3,byrow=TRUE))
spVineCopula(spCop, vineCopula(RVM))
}
\keyword{ mulitvariate }
\keyword{ distribution }