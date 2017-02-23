\name{asCopula}
\alias{asCopula}

\title{Constructor of an asymmetric copula with cubic and quadratic sections (Nelsen 2006).}
\description{A copula from the copula family of an asymmetric copula with cubic and quadratic sections for a given pair of parameters \code{param} (Nelsen 2006).}
\usage{
asCopula(param)
}
\arguments{
  \item{param}{
The copula is defined by the two-place parameter \code{param}.
}
}

\value{
An \code{\linkS4class{asCopula}} object that behaves (almost) like any other \code{\linkS4class{copula}} from the \code{copula} package.
}
\references{
Example 3.16 in: Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer
}
\author{
Benedikt Graeler
}

\seealso{ 
See also \code{\link{cqsCopula}}.
}
\examples{
cop <- asCopula(c(-2.1,0.8))
cop

smpl <- rCopula(100, cop)

par(mfrow=c(2,2), mai=c(0.4,0.4,0.4,0.4))
plot(smpl,asp=1, main="sample")
persp(cop, dCopula, asp=1, main="density")

contour(cop, pCopula, asp=1, main="CDF")
points(rCopula_y(0.4, cop, 10), col="red")

copFit <- fitCopula(cop, smpl)
tau(cop)
rho(cop)
lambda(cop)

# too few points, but faster for package compilation
kenCop <- getKendallDistr(cop, smpl)
curve(kenCop, main="Kendall function", asp=1)
}

\keyword{asymmetric}
\keyword{copula}
\keyword{cubic quadratic sections}
