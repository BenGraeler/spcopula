\name{trunCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{trunCopula-class}
\alias{fitCopula,trunCopula,matrix-method}
\alias{qCopula_v,trunCopula-method}
\alias{rCopula_y,numeric,trunCopula-method}

\title{Class \code{"trunCopula"}}
\description{
A representation of a truncated copula where mass below a truncation boundary is pushed on to the the boundary.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("trunCopula", ...)}. Or via the simplified constructor \code{\link{trunCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{trunFamily}:}{Object of class \code{"copula"} of the underlying continuous copula.}
    \item{\code{contPar}:}{Object of class \code{"numeric"} between 1 and infintiy defining the power 'p' of the truncation boundary 'u^p'.}
    \item{\code{.tools}:}{Object of class \code{"list"} containing a couple of helper functions that are assigned via the constructor.}
    \item{\code{dimension}:}{Object of class \code{"integer"} giving the dimesnion (currently only 2). }
    \item{\code{parameters}:}{Object of class \code{"numeric"} representing the continuous copula parameters and the power 'p' as of 'contPar'.}
    \item{\code{param.names}:}{Object of class \code{"character"} giving the parameter names. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} giving the parameter lower bounds.}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} giving the parameter upper bounds.}
    \item{\code{fullname}:}{Object of class \code{"character"} <deprecated>. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{parCopula}"}, by class "copula", distance 2.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 3.
}
\section{Methods}{
  \describe{
    \item{fitCopula}{\code{signature(copula = "trunCopula", data = "matrix")}: ... }
    \item{qCopula_v}{\code{signature(copula = "trunCopula", p = "numeric", v = "numeric")}: ... }
    \item{rCopula_y}{\code{signature(y = "numeric", copula = "trunCopula")}: ... }
	 }
}
\author{
Benedikt Graeler
}

\examples{
showClass("trunCopula")
}

\keyword{classes}
