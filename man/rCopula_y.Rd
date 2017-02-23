\name{rCopula_y}
\alias{rCopula_y}
\alias{rCopula_y-methods}
\alias{rCopula_y,ANY-method}
\alias{rCopula_y,numeric,copula-method}

\title{
Sampling from a given contour level
}
\description{
This method draws a sample from a given copula for a given contour level 
respecting the copula density along the contour.
}
\usage{
rCopula_y(y, copula, n = 1, n.disc = 100)
}

\arguments{
  \item{y}{
A vector of contour levels, i.e. y = C(u,v).
}
  \item{copula}{
a copula object that is used
}
  \item{n}{
The number of samples; can only be larger 1 if a single contour level \code{"y"} is provided.
}
  \item{n.disc}{
An integer giving the approximate number of evaluation points of the contour curve(s) per axis. 
}
}
\details{
A) The contour curve is evaluated for each value of \code{y}. The length of this contour curve is approximated as piece-wise linear and the copula density along this curve is integrated. 
B) Taking am appropriate uniform random number and applying the inverse of the (almost) CDF yields a distance d on the contour curve. Based on the piece-wise linear function, a pair (u,v) is found that represents the above distance d. Depending if the pair lies in the upper or lower triangle, u or v is fixed and \code{\link{qCopula_u}} respectively \code{\link{qCopula_v}} is applied to find the "exact" random pair on the contour curve. The procedure starting in A is repeated for each value in \code{y} or the procedure starting in B is repeated \code{n}-times if only a single \code{y} is provided.
}
\value{
A matrix of random pairs with as many rows as values in \code{y} or \code{n} if only a single value \code{y} is given.
}
\author{
Benedikt Graeler
}


\seealso{
\code{\link{rCopula}}
}
\examples{
contour(frankCopula(4), pCopula, asp=1)
points(rCopula_y(0.4, frankCopula(4), 100), asp=1, xlim=c(0.8,1), col="red")
}
\keyword{ multivariate }
\keyword{ distribution }
