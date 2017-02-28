\name{trunCopula}
\alias{trunCopula}

\title{
Constructor of the truncated copula class.
}
\description{
Returns a class \code{\linkS4class{trunCopula}} representing a truncated copula where mass of a continuous copula is concentrated in the boubndary curve. 
}
\usage{
trunCopula(copula, contPar, approx.u = 1:1000/1000)
}

\arguments{
  \item{copula}{the underlying continuous copula}
  \item{contPar}{The power 'p' of the polynomial boundary 'x^p'.}
  \item{approx.u}{a discretisation of the u-axis for numerical approximations}
}
\value{An object of \code{\linkS4class{trunCopula}}}
\author{
Benedikt Graeler}

\examples{

cop <- trunCopula(gumbelCopula(3), 1.2)
cop

smpl <- rCopula(100, cop)

par(mfrow=c(2,2), mai=c(0.4,0.4,0.4,0.4))
plot(smpl,asp=1, main="sample")

contour(cop, pCopula, asp=1, main="CDF", n=201)
points(rCopula_y(0.4, cop, 10), col="red")

# too few points, but faster for package compilation
kenCop <- getKendallDistr(cop, smpl)
curve(kenCop, main="Kendall function", asp=1)
}