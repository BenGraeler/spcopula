\name{condCovariate}
\alias{condCovariate}

\title{
Conditioning of a Covariate
}
\description{
Conditions the covariate on the central location in a \code{\linkS4class{stCoVarVineCopula}}.
}
\usage{
condCovariate(stNeigh, coVarCop)
}

\arguments{
  \item{stNeigh}{a spatio-temporal neighbourhood \code{\linkS4class{stNeighbourhood}} with a covariate.}
  \item{coVarCop}{a function taking spatial and temporal indices and returning a \code{\linkS4class{copula}} object describing the dependence between variable of interest and covariate at this location in space and time.}
}

\value{A vector of conditioned data, i.e. covariate|variable of interest}

\author{
Benedikt Graeler
}

\seealso{\code{\linkS4class{stNeighbourhood}}}

\examples{
library("sp")
library("spacetime")

sp <- SpatialPoints(matrix(c(181000,181100,333500,333600),2))
time <- as.POSIXct("2014-03-18")+60*60*24*c(0,1,2)
data <- data.frame(var=runif(6))
data$coVar <- invdduCopula(data$var, gumbelCopula(7), runif(6)) 

stData <- STFDF(sp, time, data)
stQuerry <- STF(SpatialPoints(matrix(c(181000,181200,333600,333600),2)),
                time[2:3])

stNeigh <- getStNeighbours(stData=stData, ST=stQuerry, 
                           spSize=3, tlags=-(0:1),
                           var="var", coVar="coVar", prediction=TRUE)

condCovariate(stNeigh, function(x) gumbelCopula(7))
}
