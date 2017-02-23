\name{kendallDistribution}
\alias{kendallDistribution}
\alias{kendallDistribution,claytonCopula-method}
\alias{kendallDistribution,frankCopula-method}
\alias{kendallDistribution,gumbelCopula-method}
\alias{kendallDistribution,BB1Copula-method}
\alias{kendallDistribution,BB6Copula-method}
\alias{kendallDistribution,BB7Copula-method}
\alias{kendallDistribution,BB8Copula-method}
\alias{kendallDistribution,joeBiCopula-method}
\title{
The Kendall distribution
}
\description{
The Kendall distribution for a given copula.
}
\usage{
kendallDistribution(copula, t)
}
\arguments{
  \item{copula}{
A copula of which the Kendall Distribution should be calculated.
}
  \item{t}{
The quantile of the copula for which the Kendall distribution should be calculated.
}
}
\value{
A vector of Kendall distribution values the same length as \code{t}.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{getKendallDistr}}}
\examples{
kendallDistribution(frankCopula(.8),.8)
}