\name{spCopPredict}
\alias{spCopPredict}

\title{
spatial prediction based on a spatial vine copula
}
\description{
A spatial vine copula is used to predict values at unobserved locations conditioned on observations of a local neighbourhood.
}
\usage{
spCopPredict(predNeigh, dataLocs, predLocs, spVine, margin,
             method = "quantile", p = 0.5, ...)
}
\arguments{
  \item{predNeigh}{the \code{\linkS4class{neighbourhood}} used for prediction, its slot \code{prediction} must be \code{TRUE}.}
  \item{dataLocs}{some \code{\linkS4class{Spatial}} class providing the data for the prediction.}
  \item{predLocs}{some \code{\linkS4class{Spatial}} class providing the prediction locations.}
  \item{spVine}{the spatial vine copula describing the spatial dependence}
  \item{margin}{the marginal distribution as a list with entries named "d" for the density function (PDF), "q" for the quantile function and "p" for cumulative distribution function (CDF).}
  \item{method}{one of \code{"quantile"} or \code{"expectation"} denoting the type of predictor.}
  \item{p}{only used for the quantile predictor indicating the desired fraction the quantile should correspond to.}
  \item{...}{Further arguments passed to \code{\link{integrate}}.}
}
\details{
Predictions are done based on \code{\link{condSpVine}} through numerical integration/optimisation.
}
\value{
A \code{\link[sp:Spatial-class]{Spatial}} object of the same type as provided in the slot \code{locations} of the argument \code{predNeigh}.
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{condSpVine}}
}
\examples{
library("sp")
library("VineCopula")
data("meuse.grid")
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE

data("meuse")
coordinates(meuse) <- ~x+y
data("spCopDemo")

calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(), tCopula(),
                                  frankCopula(), normalCopula(), claytonCopula(),
                                  claytonCopula(), claytonCopula(), claytonCopula(),
                                  claytonCopula(), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

spVineCop <- spVineCopula(spCop, vineCopula(4L))

meuse$rtZinc <- rank(meuse$zinc)/(length(meuse)+1)

dataLocs <- meuse[1:4,]
predLocs <- meuse.grid[c(9:12,16:19,25:28),]
predMeuseNeigh <- getNeighbours(dataLocs, predLocs, 
                                5, "rtZinc", prediction=TRUE, min.dist=-1)

qMar <- function(x) {
  qlnorm(x,mean(log(meuse$zinc)),sd(log(meuse$zinc)))
}

predMedian <- spCopPredict(predMeuseNeigh, dataLocs, predLocs,
                           spVineCop, list(q=qMar), "quantile", p=0.5)

\dontrun{
  spplot(predMedian, "quantile.0.5", 
         sp.layout=list("sp.points", meuse, pch = 19, col = "red"),
         col.regions=bpy.colors())
}
}

\keyword{ distribution }
\keyword{ prediction }
\keyword{ spatial }
