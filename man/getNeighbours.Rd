\name{getNeighbours}
\alias{getNeighbours}

\title{
Creating Local Neighbourhoods
}
\description{
This function calculates a local neighbourhood to be used for fitting of spatial/spatio-temporal vine copulas and for prediction using spatial/spatio-temporal vine copulas.
}
\usage{
getNeighbours(dataLocs, predLocs, size = 5, var = names(dataLocs)[1], 
              coVar=character(), prediction=FALSE, min.dist = 0.01)
}
\arguments{
  \item{dataLocs}{
some spatial data frame holding the data used for estimation/prediction
}
  \item{predLocs}{A spatial object defining the prediction locations, might be missing if the neighbourhood is used for fitting}
  \item{size}{The size of the neighbourhood including the location of interest (for fitting as well for prediction)}
  \item{var}{the variable name of interest, by default the first variable is used}
  \item{coVar}{the name of the covariate, if any}
  \item{prediction}{whether the neighbourhood should be used for prediction (TRUE) or spatial/Spatio-temporal vine copula fitting}
  \item{min.dist}{the minimal distance for a location to be included. Must be larger than 0 for fitting purposes and might be 0 for prediction}
}
\value{An object of \code{\linkS4class{neighbourhood}}.}
\author{Benedikt Graeler}

\seealso{
See \code{\link{neighbourhood}} for the native constructor of a \code{\linkS4class{neighbourhood}} class.
}
\examples{
library("sp")
spdf <- data.frame(x=c(112,154,212,289,345),y=c(124,198,85,168,346),measure=rlnorm(5))
coordinates(spdf) <- ~x+y

getNeighbours(spdf,size=4)
}
\keyword{ spatial }