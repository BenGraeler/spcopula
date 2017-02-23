\name{getStNeighbours}
\alias{getStNeighbours}

\title{
Creating Local Spatio-Temporal Neighbourhoods
}
\description{
This function calculates local spatio-temporal neighbourhoods to be used for fitting of spatio-temporal vine copulas and for prediction using spatio-temporal vine copulas.
}
\usage{
getStNeighbours(stData, ST, spSize = 4, tlags=-(0:2), var = names(stData@data)[1], 
                coVar=character(), timeSteps=NA, prediction=FALSE, min.dist = 0.01)
}
\arguments{
  \item{stData}{
some spatio-temporal data frame holding the data used for estimation/prediction
}
  \item{ST}{
A spatio-temporal object defining the prediction locations, might be missing if the spatio-temporal neighbourhood is used for fitting.
}
  \item{spSize}{
The spatial size of the neighbourhood including the location of interest (for fitting as well for prediction).
}
  \item{tlags}{The temporal lags to be used in the spatio-temporal neighbourhood.}
  \item{var}{the variable name of interest, by default the first variable is used}
  \item{coVar}{the covariate's name}
  \item{timeSteps}{The number of time instances that should randomly be selected from \code{stData}. The default, \code{NA}, selects all locations.}
  \item{prediction}{whether the neighbourhood should be used for prediction (TRUE) or spatial/Spatio-temporal vine copula fitting.}
  \item{min.dist}{
the minimal distance for a location to be included. Must be larger than 0 for fitting purposes and might be 0 for prediction.
}
}
\value{
An object of \code{\linkS4class{stNeighbourhood}}.
}
\author{
Benedikt Graeler
}

\seealso{
See \code{\link{stNeighbourhood}} for the native constructor of a \code{\linkS4class{stNeighbourhood}} class. The pure spatial version can be found at \code{\link{getNeighbours}}.
}
\examples{
library("sp")
library("spacetime")

sp <- SpatialPoints(matrix(c(181000,181100,333500,333600),2))
time <- as.POSIXct("2014-03-18")+60*60*24*c(0,1,2)
data <- data.frame(measure=runif(6))

stData <- STFDF(sp, time, data)
stQuerry <- STF(SpatialPoints(matrix(c(181000,181200,333600,333600),2)),
                time[2:3])

getStNeighbours(stData=stData, ST=stQuerry, var="measure", spSize=3, 
                tlags=-(0:1), prediction=TRUE)
}
\keyword{ spatial }