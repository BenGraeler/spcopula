\name{fitSpCopula}
\alias{fitSpCopula}

\title{
Spatial Copula Fitting
}

\description{
A bivariate spatial copula is composed out of a set of bivariate copulas. These are combined using a convex linear combination with weights based on distances where for copulas with a 1-1 correspondence of Kendall's tau or Spearman's rho a dependence function providing measures of association based on distances might be used. This function estimates a spatial dependence function, evaluates the log-likelihood per family and lag class, selects the best fits and composes a spatial bivariate copula.
}
\usage{
fitSpCopula(bins, data, cutoff = NA, families = c(normalCopula(), tCopula(),
            claytonCopula(), frankCopula(), gumbelCopula()), ...)
}
\arguments{
  \item{bins}{the bins to be used, typically output from \code{\link{calcBins}}.}
  \item{data}{the spatial dataset that ahs been used to generate \code{bins}.}
  \item{cutoff}{The maximal distance to be used in the fit.}
  \item{families}{The set of families to be investigated.}
  \item{\dots}{Passed on to the function \code{\link{fitCorFun}}.}
}

\value{
A \code{\linkS4class{spCopula}} object.
}

\references{
Graeler, B. & E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. 
Poster at: Spatial Statistics 2011 - Mapping global change. Enschede, The Netherlands, 23-25 March 2011. 
DOI: 10.1016/j.proenv.2011.07.036
}
\author{
Benedikt Graeler
}

\seealso{
Take a look at \code{\link{fitCorFun}}, \code{\link{loglikByCopulasLags}}, \code{\link{composeSpCopula}} and \code{\linkS4class{spCopula}} to gain more control over the fitting procedure.
}
\examples{
# reload some spatial data
library("sp")
data("meuse")
coordinates(meuse) <- ~x+y

# drop margins
meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))

# load data from a provided binning
data("spCopDemo")

fitSpCopula(bins, meuse, 600)
}

\keyword{spatial}
\keyword{multivariate}
\keyword{distribution}
