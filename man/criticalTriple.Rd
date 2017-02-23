\name{criticalTriple}
\alias{criticalTriple}
\title{
calculate critical triples
}
\description{
For two given arguments, the third argument is calculated for a given critical level (cumulated probability level) and copula.
}
\usage{
criticalTriple(copula, cl, u, ind, tol=sqrt(.Machine$double.eps))
}
\arguments{
  \item{copula}{
The three-dimensional copula.
}
  \item{cl}{
The critical level (cumulative probability level).
}
  \item{u}{
The two given arguments.
}
  \item{ind}{
The index of the given arguments.
}
\item{tol}{
  The tolerance value as used by \code{\link{optimise}}
}
}
\value{
The third argument for which the given critical level (cumulative probability level) is achieved.
}
\author{
Benedikt Graeler
}
\note{
\code{\link{optimise}} is used to find the third argument.
}
\seealso{
\code{\link{criticalPair}}
}
\examples{
w <- criticalTriple(frankCopula(0.7,dim=3), 0.9, c(.97,.97), c(1,2))

# check the triple
pCopula(c(0.97, 0.97, w), frankCopula(0.7, dim=3))

}
\keyword{ multivariate }
\keyword{ distribution }
