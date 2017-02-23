\name{cqsCopula}
\alias{cqsCopula}
\title{
Constructor of a symmetric copula with cubic quadratic sections. 
}
\description{
A copula from the copula family of a symmetric copula with cubic and quadratic sections (\code{\linkS4class{cqsCopula}}) for a given pair of parameters \code{param} (Nelsen 2006).
}
\usage{
cqsCopula(param=c(0,0), fixed="")
}

\arguments{
  \item{param}{
The copula is defined by the two-place parameter \code{param}. The default represents the product copula. 
}
\item{fixed}{a single character referring to "a" or "b" fixing one of the parameters. This allows to estimate the remaining parameter through inversion of Kendall's tau \code{\link{iTau}} or Spearman's rho \code{\link{iRho}}.}
}

\value{
An \code{\linkS4class{cqsCopula}} object that behaves (almost) like any other \code{\linkS4class{copula}} from the \code{\link{copula-package}}.
}
\references{
Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer
}
\author{
Benedikt Graeler
}
\seealso{
See also \code{\link{asCopula}}.
}
\examples{
cop <- cqsCopula(c(-2.1,0.8))
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
\keyword{copula}
\keyword{cubic quadratic sections}
