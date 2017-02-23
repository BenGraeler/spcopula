\name{condStCoVarVine}
\alias{condStCoVarVine}

\title{
conditional distribution function of spatio-temporal covariate vine copula
}
\description{
Returns a conditional distribution function of spatio-temporal covariate vine copula
}
\usage{
condStCoVarVine(condVar, dists, stCVVC, stInd, n = 1000)
}

\arguments{
  \item{condVar}{
the conditioning variables
}
  \item{dists}{
spatio-temporal distances to the conditioning variables
}
  \item{stCVVC}{
the spatio-temporal covariate vine copula of the model
}
  \item{stInd}{
spatio-temporal index pair to be used with covariate copula (which is in first place a function taking a pair of indices and returns a copula object)
}
  \item{n}{
number of approximation points
}
}
\value{
a univariate distribution function over [0,1]
}
\author{
Benedikt Graeler
}
\note{
The distribution is linearly approximated at a limited number (\code{n}) of points.
}

\seealso{
\code{\link{condStVine}}, \code{\link{condSpVine}}
}

\examples{
library("VineCopula")
spCopT0 <- spCopula(components=list(claytonCopula(8), claytonCopula(4), 
                                    claytonCopula(2), claytonCopula(1),
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(100,200,300,400,500,600),
                    unit="km")
spCopT1 <- spCopula(components=list(claytonCopula(4), claytonCopula(2), 
                                    claytonCopula(1), claytonCopula(0.5),
                                    indepCopula()),
                    distances=c(100,200,300,400,500),
                    unit="km")
spCopT2 <- spCopula(components=list(claytonCopula(2), claytonCopula(1), 
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(100,200,300,400),
                    unit="km")

stCop <- stCopula(components=list(spCopT0, spCopT1, spCopT2),
                  tlags=-(0:2))

# only a constant copula ius used for the covariate
stCVVC <- stCoVarVineCopula(function(x) gumbelCopula(7), stCop, vineCopula(5L))

dists <- array(c(150, 250, 150, 250,0,0,-1,-1),dim=c(1,4,2))
condVar <- c(0.95, 0.29, 0.55, 0.05, 0.41)

condDensity <- condStCoVarVine(condVar, dists, stCVVC, c(1,1))
curve(condDensity)
}