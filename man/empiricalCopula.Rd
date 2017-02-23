\name{empiricalCopula}
\alias{empiricalCopula}
\title{
Constructor of an empirical copula class
}
\description{
Constructs an object of the empirical copula class \code{\linkS4class{empiricalCopula}}, see \code{\link{genEmpCop}} for a simplified version.
}
\usage{
empiricalCopula(sample=NULL, copula)
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
An object of \code{\linkS4class{empiricalCopula}}.
}
\author{
Benedikt Graeler
}
\note{
Its implementation of \code{\link{pCopula}} is based on C-code from \code{\link{copula-package}}.
}

\seealso{
\code{\link{genEmpCop}} for a simplified constructor with sample length control.
}
\examples{
empCop <- empiricalCopula(rCopula(500,frankCopula(0.7)))
str(empCop)

empCop <- empiricalCopula(copula=frankCopula(0.7))
str(empCop)

empCop <- empiricalCopula(rCopula(500,frankCopula(0.7)), frankCopula(0.7))
str(empCop)

# the empirical value
pCopula(c(0.3, 0.5), empCop)

# the theoretical value
pCopula(c(0.3, 0.5), frankCopula(0.7))
}
\keyword{ multivariate }
