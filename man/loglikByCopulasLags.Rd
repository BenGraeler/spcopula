\name{loglikByCopulasLags}
\alias{loglikByCopulasLags}

\title{Log-likelihoods by copula family and spatial lag class}

\description{
This function calculates the log-likelihood for a set of provided copula families per lag class. The copulas' parameters are either fitted by a provided distance dependent function, or through the function \code{\link{fitCopula}} per lag and copula family. 
}
\usage{
loglikByCopulasLags(bins, data, families = c(normalCopula(),
                                             tCopula(),
                                             claytonCopula(),
                                             frankCopula(),
                                             gumbelCopula()),
                    calcCor, lagSub=1:length(bins$meanDists))
}

\arguments{
  \item{bins}{Typically, the output of the \code{\link{calcBins}} function providing columns \code{meanDists} to tune the copula's parameter and \code{lagData} to evaluate the log-likelihood.}
  \item{data}{the spatial data set used to derive the bins from for data retrieval.}
  \item{families}{A vector of representatives of the bivariate copula families that should be investigated.}
  \item{calcCor}{The tuning function for the copulas' parameters. Typically, the output of the \code{\link{fitCorFun}} function. Its either based on Kendall's tau, Spearman's rho or the identity. This function needs to return one of "kendall", "spearman" or "id" when passed "NULL". This parameter might be "missing" for static convex spatial copulas.}
  \item{lagSub}{A set of indices indicating which of the bins to use.}
}

\value{A list containing a matrix (\code{loglik}) of spatial lags (rows) and copula family names (columns) holding the calculated log-likelihood value and a list of the corresponding copula fits.}

\author{Benedikt Graeler}

\seealso{
See also \code{\link{calcBins}}, \code{\link{fitCorFun}} or \code{\link{spCopula}}.
}

\examples{
# load spatial data
library("sp")
data("meuse")
coordinates(meuse) <- ~x+y

# drop margins
meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))

# load parts of the workspace from demo(spCopula)
data("spCopDemo") 

calcKTauPol <- fitCorFun(bins, degree=3)

loglikTau <- loglikByCopulasLags(bins, meuse, calcCor=calcKTauPol)
loglikTau$loglik
}

\keyword{spcopula}
\keyword{log-likelihood}