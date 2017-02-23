\name{stVineCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{stVineCopula-class}
\alias{fitCopula,stVineCopula-method}

\title{Class \code{"stVineCopula"}}
\description{
An object representing a spatio-temporal vine copula.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stVineCopula", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{stCop}:}{Object of class \code{"stCopula"} the bivariate spatio-temporal copula in the first tree. }
    \item{\code{topCop}:}{Object of class \code{"copula"} the remaining vine copula. }
    \item{\code{dimension}:}{Object of class \code{"integer"} giving the dimension of the spatio-temporal copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.names}:}{Object of class \code{"character"} empty }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{fullname}:}{Object of class \code{"character"} describing the copula family. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{fitCopula}{\code{signature(copula = "stVineCopula")}: to fit a spatio-temporal vine based on a spatio-temporal neighbourhood (\code{\linkS4class{stNeighbourhood}}) and a spatio-temporal bivariate copula (\code{\linkS4class{stCopula}}).}
	 }
}
\author{
Benedikt Graeler
}

\seealso{\code{\linkS4class{spVineCopula}}}
\examples{
showClass("stVineCopula")
}
\keyword{classes}