\name{spGaussCopPredict}
\alias{spGaussCopPredict}
\title{
spatial prediction using a Gaussian Copula
}
\description{
This function allows to do spatial prediction using a Gaussian Copula based on the spatial parametrization.
}
\usage{
spGaussCopPredict(corFun, predNeigh, dataLocs, predLocs, margin, p = 0.5, ..., n = 1000)
}

\arguments{
  \item{corFun}{
A valid correlogram (i.e. producing a valid correlation matrix; e.g. based on a variogram).
}
  \item{predNeigh}{
A \code{\linkS4class{neighbourhood}} object used for prediction.
}
  \item{dataLocs}{a spatial object providing the data}
  \item{predLocs}{a spatial object providing the prediction locations}
  \item{margin}{
a list containing the marginal distribution. Currently only the entry \code{q} is required defining the quantile function.
}
  \item{p}{
the fraction the quantile function shall be evaluated for. This can be used to calculate besides the median estimate confidence estimates as well.
}
  \item{\dots}{
currently unused
}
  \item{n}{
the approximate number of points used in the linear approximation of the conditional distribution function.
}
}
\details{
Based on \code{corFun} provided with a distance matrix a Gaussian copula (\code{\linkS4class{normalCopula}}) is generated and conditioned under the data of the neighbouring locations. The 1-dimensional conditional distribution is approximated and evaluated for the given fraction \code{p}. This conditioned fraction is than passed on to the marginal quantile function and evaluated providing an estimate.
}
\value{
According to the \code{predLocs} slot of the provided \code{\linkS4class{neighbourhood}} a spatial data structure extended with an variable holding the predicted values.
}
\author{
Benedikt Graeler
}


\seealso{
\code{\link{spCopPredict}}
}
\examples{
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

# get some prediction data
library("sp")
data("meuse.grid")
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE

data("meuse")
coordinates(meuse) <- ~x+y

meuse$rtZinc <- rank(meuse$zinc)/(length(meuse)+1)

# obtain the prediction neighbourhoods
predMeuseNeigh <- getNeighbours(meuse[1:4,], meuse.grid[c(9:12,16:19,25:28),],
                                5, "rtZinc", prediction=TRUE, min.dist=-1)

qMar <- function(x) {
  qlnorm(x,mean(log(meuse$zinc)),sd(log(meuse$zinc)))
}

# predict using the Gaussian Copula
predMedian <- spGaussCopPredict(meuseGaussCorFun, predMeuseNeigh,
                                meuse[1:4,], meuse.grid[c(9:12,16:19,25:28),],
                                list(q=qMar))

\dontrun{
  spplot(predMedian,"quantile.0.5", 
         sp.layout=list("sp.points", meuse, pch = 19, col = "red"),
         col.regions=bpy.colors())
}
}

\keyword{ prediction }
\keyword{ distribution }