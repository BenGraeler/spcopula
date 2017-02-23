\name{reduceNeighbours}
\alias{reduceNeighbours}

\title{
Selecting the strongest correlated neighbours
}
\description{
A function selecting the strongest correlated neighbours from a larger set of neighbours
}
\usage{
reduceNeighbours(stNeigh, stDepFun, n,
                 prediction=stNeigh@prediction, dropEmpty=!prediction)
}

\arguments{
  \item{stNeigh}{
the proxy neighbourhood to be investigated
}
  \item{stDepFun}{
a spatio-temporal dependence function that return correlation estimates based on a spatial and temporal distance
}
  \item{n}{the number of neighbours to be selected.}
  \item{prediction}{whether the neighbourhood is used for prediction (the data slot does not provide the central location's data)}
  \item{dropEmpty}{whether empty neighbourhoods (i.e. neighbourhoods with \code{NA}'s at the central location) shall be dropped.}
}
\value{A spatio-temporal neighbourhood \code{\linkS4class{stNeighbourhood}} with fewer neighbours.}

\author{
Benedikt Graeler
}


\seealso{
\code{\link{getStNeighbours}}
}

\examples{
library("sp")
library("spacetime")

sp <- SpatialPoints(matrix(c(181000,181100,333500,333600),2))
time <- as.POSIXct("2014-03-18")+60*60*24*c(0,1,2,3,4)
data <- data.frame(var1=runif(10))

stData <- STFDF(sp, time, data)

stNeigh <- getStNeighbours(stData, spSize=2, tlags=-(0:2))

reduceNeighbours(stNeigh, function(h,delta) return(1/h/delta), 2)
}