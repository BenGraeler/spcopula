\name{tawn3pCopula}
\alias{tawn3pCopula}

\title{
Tawn Copula Family constructor using all three parameters
}
\description{
Constructor of the Tawn Copula family with three parameters allowing for asymmetric copula families.
}
\usage{
tawn3pCopula(param = c(0.5, 0.5, 2))
}

\arguments{
  \item{param}{
  The vector defining the three parameters of the Tawn Copula.
}
}
\details{
Pickand's dependece function is given through: 
A(t) = (1-beta)*t + (1-alpha)*(1-t) + ((alpha*(1-t))^theta+(beta*t)^theta)^(1/theta)
}
\value{
An instance of the \code{\linkS4class{tawn3pCopula}} class.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{tawnCopula}} for a symmetric one parameter, \code{\link{tawnT1Copula}} and \code{\link{tawnT2Copula}} for asymmetric two-parameter sub-families of the Tawn copula family. 
}

\examples{
persp(tawn3pCopula(c(0.4,0.9,4)), dCopula)
persp(tawn3pCopula(c(0.4,0.9,4)), pCopula)

}

\keyword{ copula }
\keyword{ distribution }
