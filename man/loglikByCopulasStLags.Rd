\name{loglikByCopulasStLags}
\alias{loglikByCopulasStLags}

\title{Log-likelihoods by copula family and spatio-temporal lag class}

\description{This function works through a set of copula families and evaluates the best fitting one for each spatio-temporal lag.}

\usage{
loglikByCopulasStLags(stBins, data, families = c(normalCopula(),
                                             tCopula(),
                                             claytonCopula(),
                                             frankCopula(),
                                             gumbelCopula()),
                      calcCor, lagSub = 1:length(stBins$meanDists))
}

\arguments{
  \item{stBins}{a spatio-temporal bining typically prodsuced with \code{\link{calcBins}}.}
  \item{data}{the spatial data set used to derive the bins from for data retrieval.}
  \item{families}{A vector of representatives of the bivariate copula families that should be investigated.}
  \item{calcCor}{a two place function tuning the copulas' parameters. Typically, the output of the \code{\link{fitCorFun}} function. Its either based on Kendall's tau, Spearman's rho or the identity. This parameter might be "missing" for static convex spatial copulas.}
  \item{lagSub}{A set of indices indicating which of the bins to use.}
}

\value{A list containing the result of \code{\link{loglikByCopulasLags}} per temporal lag: i.e. a list containing a matrix (\code{loglik}) of spatial lags (rows) and copula family names (columns) holding the calculated log-likelihood value and a list of the corresponding copula fits.}

\author{Benedikt Graeler}

\seealso{\code{\link{calcBins}}, \code{\link{loglikByCopulasLags}}}

\examples{
# the spatial case
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
