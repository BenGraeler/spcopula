\name{criticalPair}
\alias{criticalPair}
\title{
Calculate Critical Pairs
}
\description{
For a given argument, the second argument is calculated for a given critical level (cumulated probability level) and copula.
}
\usage{
criticalPair(copula, cl, u, ind, tol=sqrt(.Machine$double.eps))
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
The second argument for which the given critical level (cumulative probability level) is achieved.
}
\author{
Benedikt Graeler
}
\note{
\code{\link{optimise}} is used to find the third argument.
}

\seealso{
\code{\link{criticalTriple}}
}
\examples{
v <- criticalPair(frankCopula(0.7), 0.9, u=.97, 1)
pCopula(c(0.97, v),frankCopula(0.7))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
