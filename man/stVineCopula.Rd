\name{stVineCopula}
\alias{stVineCopula}
\title{
Constructor of the class \code{\linkS4class{stVineCopula}}.
}
\description{
Construct an instance of the \code{\linkS4class{stVineCopula}}.
}
\usage{
stVineCopula(stCop, topCop)
}
\arguments{
  \item{stCop}{Object of class \code{\linkS4class{stCopula}}. A bivariate spatio-temporal copula. }
  \item{topCop}{Object of class \code{\linkS4class{copula}}. A bivariate, multivariate or R-Vine copula joining the spatio-temporal components (neighbours). In case all trees are spatio-temporal copulas, this argument might be \code{NULL}.}
}
\value{
An instance of the \code{\linkS4class{stVineCopula}} class.
}
\author{
Benedikt Graeler
}
\examples{
# a spatio-temporal C-vine copula (with independent copulas in the upper vine)
library("VineCopula")

spCopT0 <- spCopula(components=list(claytonCopula(8), claytonCopula(4), 
                                    claytonCopula(2), claytonCopula(1),
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300,400,500),
                    unit="km")
spCopT1 <- spCopula(components=list(claytonCopula(4), claytonCopula(2), 
                                    claytonCopula(1), claytonCopula(0.5),
                                    indepCopula()),
                    distances=c(0,100,200,300,400),
                    unit="km")
spCopT2 <- spCopula(components=list(claytonCopula(2), claytonCopula(1), 
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300),
                    unit="km")

stCop <- stCopula(components=list(spCopT0, spCopT1, spCopT2),
                  tlags=-(0:2))

stVineCopula(stCop, vineCopula(9L))
}
\keyword{ mulitvariate }
\keyword{ distribution }