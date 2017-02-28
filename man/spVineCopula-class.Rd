\name{spVineCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{spVineCopula-class}
\alias{mixedSpVineCopula-class}
\alias{pureSpVineCopula-class}
\alias{fitCopula,spVineCopula,ANY-method}

\title{Class \code{"spVineCopula"}}
\description{
A S4-class representing a spatial vine copula structure to model spatial processes. This spatial vine might be composed out of purely spatial trees or being mixed as lower spatial trees with an upper non-spatially varying joined copula.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("spVineCopula", ...)} or to the constructor \code{\link{spVineCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{spCop}:}{Object of class \code{\linkS4class{spCopula}} A bivariate spatially varying copula. }
    \item{\code{topCop}:}{Object of class \code{\linkS4class{vineCopula}}, \code{\linkS4class{copula}} or \code{NULL}. The top copula being constant over space. It is either a R-Vine (\code{\linkS4class{vineCopula}}) or a bivariate copula (\code{\linkS4class{copula}}) joining the spatial components (neighbours) or all trees of the vine are spatial bivariate copulas (\code{NULL}).}
    \item{\code{dimension}:}{Object of class \code{"integer"} giving the dimension of the spVineCopula (neighbourhood). }
    \item{\code{parameters}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.names}:}{Object of class \code{"character"} empty }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{fullname}:}{Object of class \code{"character"} A descriptive name of the this copula family.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
No methods defined with class \code{spVineCopula} in the signature.
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{spCopula}}, \code{\linkS4class{vineCopula}}
}
\examples{
showClass("spVineCopula")
}
\keyword{classes}
