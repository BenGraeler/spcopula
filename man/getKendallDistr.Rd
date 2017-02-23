\name{getKendallDistr}
\alias{getKendallDistr}
\alias{getKendallDistr,claytonCopula-method}
\alias{getKendallDistr,frankCopula-method}
\alias{getKendallDistr,gumbelCopula-method}
\alias{getKendallDistr,BB1Copula-method}
\alias{getKendallDistr,BB6Copula-method}
\alias{getKendallDistr,BB7Copula-method}
\alias{getKendallDistr,BB8Copula-method}
\alias{getKendallDistr,joeBiCopula-method}

\title{
Retrieving the Kendall Distribution function for a given copula
}
\description{
Depending on the copula passed, an analytical or empirical Kendall distribution function (\code{\link{kendallDistribution}}) is returned.
}
\usage{
getKendallDistr(copula, sample = NULL)
}
\arguments{
  \item{copula}{
  A copula object.
}
  \item{sample}{
  In case an empirical Kendall Distribution function has to be used, one can pass a sample of the copula that will then be used to derive the empirical Kendall distribution function. Of it is not provided, it will be sampled on the fly. Note, that fixing a sample might speed up computations considerably.
}
}
\value{
A function that takes quantiles of the copula and returns the Kendall Distribution value for this quantile.}
\author{
Benedikt Graeler
}


\seealso{
\code{\link{kendallDistribution}}
}
\examples{
frankKenDistrFun <- getKendallDistr(frankCopula(.5))
frankKenDistrFun(.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
