\name{empSurCopula}
\alias{empSurCopula}
\title{
Constructor of an empirical survival copula class
}
\description{
Constructs an object of the empirical survival copula class \code{\linkS4class{empSurCopula}}, see \code{\link{genEmpSurCop}} for a simplified version.
}
\usage{
empSurCopula(sample=NULL, copula)
}
\arguments{
  \item{sample}{
A sample from a provided or unknown copula family.
}
  \item{copula}{
The underlying theoretical copula, in case it is known or a sample should be generated.
}
}
\value{
An object of \code{\linkS4class{empSurCopula}}.
}
\author{
Benedikt Graeler
}
\note{
Its implementation of \code{\link{pCopula}} is based on C-code from \code{\link{copula-package}}.
}

\seealso{
\code{\link{genEmpSurCop}} for a simplified constructor with sample length control.
}
\examples{
empCop <- empSurCopula(rCopula(500,frankCopula(0.7)))
str(empCop)

empCop <- empSurCopula(copula=frankCopula(0.7))
str(empCop)

empCop <- empSurCopula(rCopula(500,frankCopula(0.7)), frankCopula(0.7))
str(empCop)

# the empirical value
pCopula(c(0.3, 0.5), empCop)

# the theoretical value
pCopula(c(0.3, 0.5), frankCopula(0.7))
}
\keyword{ multivariate }
