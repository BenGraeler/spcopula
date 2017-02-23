\name{dduCopula}
\alias{dduCopula}
\alias{ddvCopula}

\title{
partial derivatives of copulas
}
\description{
Similar to the \code{\link{dCopula}} and \code{\link{pCopula}} the function \code{dduCopula} evaluates the partial derivative \eqn{\frac{\partial}{\partial u} C(u,v)} and the function \code{ddvCopula} evaluates the partial derivative \eqn{\frac{\partial}{\partial v} C(u,v)} of the provided copula. 
}
\usage{
dduCopula(u, copula, ...)
ddvCopula(u, copula, ...)
}

\arguments{
  \item{u}{
Pairs of values for which the partial derivative should be evaluated.
}
  \item{copula}{
The copula object representing the family member of interest.
}
  \item{\dots}{additional arguments can be passed on to the underlying functions.}
}
\value{
A vector of the evaluated partial derivatives of the same length as rows in \code{u}.
}
\author{
Benedikt Graeler
}

\examples{
####################################
## Asymmetric vs. Gaussian copula ##
####################################

asCop <- asCopula(c(-2,1))
asCopSmpl <- rCopula(100,asCop)

unitScatter(smpl=asCopSmpl)

# conditional probabilities of an asymmetric copula given u
asGivenU <- dduCopula(asCopSmpl,asCop)

# vs. conditional probabilities of an asymmetric copula given v
asGivenV <- ddvCopula(asCopSmpl[,c(2,1)],asCop)
unitScatter(smpl=cbind(asGivenU, asGivenV))

normalCop <- normalCopula(.6)
normCopSmpl <- rCopula(100,normalCop)

unitScatter(smpl=normCopSmpl)

# conditional probabilities of a Gaussian copula given u
normGivenU <- dduCopula(normCopSmpl,normalCop)

# vs. conditional probabilities of a Gaussian copula given v
normGivenV <- ddvCopula(normCopSmpl[,c(2,1)],normalCop)
unitScatter(smpl=cbind(normGivenU, normGivenV))
}

\keyword{partial derivative}
\keyword{conditional probabilities}
