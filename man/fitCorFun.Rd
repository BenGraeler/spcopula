\name{fitCorFun}
\alias{fitCorFun}

\title{
Automated fitting of a correlation function to the correlogram
}
\description{
Polynomials of different degrees can be fitted to the spatial/spatio-temporal correlogram calculated using \code{\link{calcBins}}. This function will be used to adjust the copula parameter in the spatial/spatio-temporal copula.
}
\usage{
fitCorFun(bins, degree = 3, cutoff = NA, tlags, bounds = c(0, 1), 
          cor.method = NULL, weighted = FALSE)
}
\arguments{
  \item{bins}{Typically the output of \code{\link{calcBins}}. Any \code{data.frame} with a columns \code{lagCor} and \code{meanDists} in the first two columns will do.}
  \item{degree}{The degree of polynomial to be fitted - recycled if needed.}
  \item{cutoff}{Maximal distance to which lags should be included in the polynomial fit.}
  \item{tlags}{The temporal lags used for the genration of \code{bins}.}
  \item{bounds}{Bounds of the correlation values. The default is set [0,1] not allowing any negative relationship but perfect positive dependence.}
  \item{cor.method}{The output of \code{\link{calcBins}} has an attribute \code{cor.method}, in case this is not present the parameter \code{cor.method} will be used. In case the parameter \code{cor.method} is not \code{NULL} and the attribute \code{cor.method} is present, they will be compared.}
  \item{weighted}{shall the residuals be weighted by the number of points in the lag class?}
}

\value{Returns a one/two-place function that provides correlation estimates for every separating spatial/spatio-temporal distance.}

\author{Benedikt Graeler}

\seealso{See also \code{\link{calcBins}} and \code{\link{spCopula}}.}

\examples{
# a simplified bins object (from demo(spcopula))
bins <- list(meanDists=c(64, 128, 203, 281, 361, 442, 522, 602, 681, 760), 
             lagCor=c(0.57,  0.49, 0.32, 0.29, 0.15, 0.14, 0.10, -0.00, 0.03, -0.01))
attr(bins,"cor.method") <- "kendall"

# plot the correlogram
plot(lagCor~meanDists,bins)

# fit and plot a linear model
calcKTauLin <- fitCorFun(bins, degree=1, cutoff=600)
curve(calcKTauLin,0, 1000, col="red",add=TRUE)

# fit and plot a polynomial model
calcKTauPol <- fitCorFun(bins, degree=5)
curve(calcKTauPol,0, 1000, col="purple",add=TRUE)
}
\keyword{correlogram}
\keyword{spcopula}
