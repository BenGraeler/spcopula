\name{condStVine}
\alias{condStVine}

\title{
Conditions a spatio-temporal vine copula for conditional prediction
}
\description{
A spatio-temporal vine copula is conditioned under the observations of all but one neighbour generating a conditional univariate distribution used for prediction.
}
\usage{
condStVine(condVar, dists, stVine, n = 1000)
}
\arguments{
  \item{condVar}{
the conditional variables
}
  \item{dists}{
the distances between the neighbours to adjust the spatial vine copula parameters.
}
  \item{stVine}{
the spatio-temporal vine copula
}
  \item{n}{
a proxy to the number of intervals to be used in the numerical integration. The intervals will not be split uniform with shorter intervals at the corners of the copula.
}
}
\value{
A function describing the conditional density.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{stVineCopula}}, \code{\link{condSpVine}}
}
\examples{
# a spatio-temporal C-vine copula (with independent copulas in the upper vine)

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

stCop <- stCopula(components=list(spCopT0, spCopT1),
                  tlags=-(0:1))

stVineCop <- stVineCopula(stCop, vineCopula(4L))

dists <- array(c(150, 250, 150, 250,0,0,-1,-1),dim=c(1,4,2))
condVar <- c(0.29, 0.55, 0.05, 0.41)

condDensity <- condStVine(condVar,dists,stVineCop)
curve(condDensity)
}

\keyword{ distribution }