\name{mixtureCopula}
\alias{mixtureCopula}

\title{
Constructor of a mixture copula}
\description{
An instance of a convex mixture of two copulas is generated.
}
\usage{
mixtureCopula(param = c(0.2, 0.2, 0.5), memberCops = c(normalCopula(0), claytonCopula(1)))
}

\arguments{
  \item{param}{A numeric vector: the set of parameters for the first and second copula appended by the mixture coefficient.}
  \item{memberCops}{A list of two copulas.}
}
\value{
An object of class \linkS4class{mixtureCopula}}
\author{
Benedikt Graeler}

\examples{
cop <- mixtureCopula(c(5,-0.8,2,0.5),
                     memberCops = list(gumbelCopula(), tCopula()))
cop

smpl <- rCopula(100, cop)

par(mfrow=c(2,2), mai=c(0.4,0.4,0.4,0.4))
plot(smpl,asp=1, main="sample")
persp(cop, dCopula, asp=1, main="density")

contour(cop, pCopula, asp=1, main="CDF")
points(rCopula_y(0.4, cop, 10), col="red")

tau(cop)
rho(cop)
lambda(cop)

# too few points, but faster for package compilation
kenCop <- getKendallDistr(cop, smpl)
curve(kenCop, main="Kendall function", asp=1)}
