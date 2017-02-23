\name{spGaussLogLik}
\alias{spGaussLogLik}

\title{
Density evaluation for a spatial Gaussian Copula
}
\description{
Evaluates the density for a spatial Gaussian Copula.
}
\usage{
spGaussLogLik(corFun, neigh, dataLocs, log = T)
}

\arguments{
  \item{corFun}{A valid correlogram (i.e. producing a valid correlation matrix; e.g. based on a variogram).}
  \item{neigh}{A \code{\linkS4class{neighbourhood}} object to be evaluated.}
  \item{dataLocs}{The same \code{\linkS4class{Spatial}} object used to generate \code{neigh}.}
  \item{log}{Should the log-likelihood be returned?}
}
\details{
Evaluates the density for all neighbourhoods in \code{neigh} and returns the (log)-likelihood.
}
\value{
The (log)-likelihood value.
}
\author{
Benedikt Graeler
}

\examples{
library("spcopula")

# load data from the Meuse demo
data("spCopDemo")

# calculate the correlation function based on Kendall's tau
calcKTauPol <- fitCorFun(bins, degree=1)

# translate Kendall's tau correlation function into Gaussian Copula parameters 
# using a linear variogram
meuseGaussCorFun <- function(h) {
  res <- pmax(iTau(normalCopula(0), calcKTauPol(0))/658*(658-h),0)
  res[h ==0] <- 1
  return(res)
}

# get the neighbours
library("sp")
data("meuse")
coordinates(meuse) <- ~x+y

meuse$rtZinc <- rank(meuse$zinc)/(length(meuse)+1)

meuseNeigh <- getNeighbours(meuse, size=5L, var="rtZinc",
                            prediction=FALSE)

# calculate the log-likelihood
spGaussLogLik(meuseGaussCorFun, meuseNeigh, meuse)
}

\keyword{ distribution }
