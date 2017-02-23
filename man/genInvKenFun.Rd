\name{genInvKenFun}
\alias{genInvKenFun}
\title{
Generate the inverse Kendall distribution function
}
\description{
The inverse of a (empirical) Kendall distribution function is generated based on numerical inversion using optimise.
}
\usage{
genInvKenFun(kenFun, tol)
}
\arguments{
  \item{kenFun}{
The (empirical) Kendall distribution function to be inverted.
}
  \item{tol}{
Tolerance passed on to \code{\link{optimise}}.
}
}
\value{
The inverse of a Kendall distribution function as a function.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{genEmpKenFun}} and \code{\link{getKendallDistr}} for the Kendall distribution function.
}
\examples{
frankKenDistrFun <- getKendallDistr(frankCopula(.5))
frankInvKenDistrFun <- genInvKenFun(frankKenDistrFun)

frankInvKenDistrFun(.8)
frankKenDistrFun(frankInvKenDistrFun(.8))
}
\keyword{ multivariate }
\keyword{ distribution }
