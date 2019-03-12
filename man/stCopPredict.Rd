\name{stCopPredict}
\alias{stCopPredict}

\title{
Spatio-Temporal Prediction based on a Spatio-Temporal Vine Copula
}
\description{
A spatio-temporal vine copula is used to predict values at unobserved spatio-temporal locations conditioned on observations of a local spatio-temporal neighbourhood.
}
\usage{
stCopPredict(predNeigh, dataST, predST, stVine, margin, method = "quantile", p = 0.5, ...)
}
\arguments{
  \item{predNeigh}{the \code{\linkS4class{neighbourhood}} used for prediction, its slot \code{prediction} must be \code{TRUE}.}
  \item{dataST}{the same \code{\linkS4class{ST}} object as used in the generation of \code{predNeigh} providing the data for interpolation.}
  \item{predST}{the same \code{\linkS4class{ST}} object as used in the generation of \code{predNeigh} providing the prediction locations for interpolation.}
  \item{stVine}{the spatio-temporal vine copula describing the spatio-temporal dependence}
  \item{margin}{the marginal distribution as a list with entries named "d" for the density function (PDF), "q" for the quantile function and "p" for cumulative distribution function (CDF).}
  \item{method}{one of \code{"quantile"} or \code{"expectation"} denoting the type of predictor.}
  \item{p}{only used for the quantile predictor indicating the desired fraction the quantile should correspond to.}
  \item{...}{Further arguments passed to \code{\link{integrate}}.}
}
\details{
Predictions are done based on \code{\link{condStVine}} through numerical integration/optimisation.
}
\value{
A \code{\link[spacetime:ST-class]{ST}} object of the same type as provided in the slot \code{locations} of the argument \code{predNeigh}.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{condStVine}} and \code{\link{spCopPredict}} for the spatial version.
}
\examples{
library("sp")
library("spacetime")
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

stCop <- stCopula(components=list(spCopT0, spCopT1),
                  tlags=-(0:1))

stVineCop <- stVineCopula(stCop, vineCopula(4L))

sp <- SpatialPoints(matrix(c(181000,181100,333500,333600),2))
time <- as.POSIXct("2014-03-18")+60*60*24*c(0,1,2)
data <- data.frame(var1=runif(6))

stData <- STFDF(sp, time, data)
stQuerry <- STF(SpatialPoints(matrix(c(181000,181200,333600,333600),2)),
                time[2:3])

stNeigh <- getStNeighbours(stData=stData, ST=stQuerry, prediction=TRUE, spSize=3,
                           tlags=-(0:1))

stCopPredict(stNeigh, stData, stQuerry, stVineCop, list(q=qunif), "quantile", 0.5)
}

\keyword{ distribution }
\keyword{ prediction }
\keyword{ spatial }
