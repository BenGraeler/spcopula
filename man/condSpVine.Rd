\name{condSpVine}
\alias{condSpVine}

\title{
Conditions a spatial vine copula for conditional prediction
}
\description{
A spatial vine copula is conditioned under the observations of all but one neighbour generating a conditional univariate distribution used for prediction.
}
\usage{
condSpVine(condVar, dists, spVine, n = 1000)
}
\arguments{
  \item{condVar}{
the conditional variables
}
  \item{dists}{
the distances between the neighbours to adjust the spatial vine copula parameters.
}
  \item{spVine}{
the spatial vine copula
}
  \item{n}{
a proxy to the number of intervals to be used in the numerical integration. The intervals will not be split uniform with a shorter intervals at the corners of the copula.
}
}
\value{
A function describing the conditional density.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{spVineCopula}}
}
\examples{
library("VineCopula")
data("spCopDemo")

calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

spVineCop <- spVineCopula(spCop, vineCopula(4L))

dists <- list(c(473, 124, 116, 649))
condVar <- c(0.29, 0.55, 0.05, 0.41)
condDensity <- condSpVine(condVar,dists,spVineCop)

curve(condDensity)
mtext(paste("Dists:",paste(round(dists[[1]],0),collapse=", ")),line=0)
mtext(paste("Cond.:",paste(round(condVar,2),collapse=", ")),line=1)
}

\keyword{ distribution }
