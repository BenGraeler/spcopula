\name{criticalLevel}
\alias{criticalLevel}
\title{
Calculating the critical level for a given Kendall Return Period
}
\description{
The critical level corresponding to a given Kendall Return Period and a given 
Kendall distribution or its underlying copula.
}
\usage{
criticalLevel(kendallFun, KRP = c(100, 1000), mu = 1, copula)
}
\arguments{
  \item{kendallFun}{
The Kendall distribution function.
}
  \item{KRP}{
The desired Kendall Return Period(s); the default are 100 and 1000.
}
  \item{mu}{
The recurrence of an event per time unit; the default is 1.
}
  \item{copula}{
The copula describing the dependencies between the investigated variables. Only necessary if \code{kendallFun} is not provided.
}
}
\value{
A vector of the same length as \code{KRP} with the corresponding critical levels.
}
\references{
Salvadori, G., De Michele, C., and Durante, F.: On the return period and design in a multivariate framework, Hydrol. Earth Syst. Sci., 15, 3293-3305, doi:10.5194/hess-15-3293-2011, 2011.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{kendallRP}} for the reverse direction: what is the Kendall Return Period for a given critical level?
}
\examples{
criticalLevel(getKendallDistr(frankCopula(.7)), KRP=c(10,100,1000))
}
\keyword{ survival }
\keyword{ multivariate }