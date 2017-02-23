\name{stCoVarVineCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{stCoVarVineCopula-class}

\title{Class \code{"stCoVarVineCopula"}}
\description{
A S4-class  representing a spatio-temporal covariate vine copula
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stCoVarVineCopula", ...)} or by calls to the constructor \code{\link{stCoVarVineCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{coVarCop}:}{Object of class \code{"function"} returning a bivariate copula object for absolute pairs of spatio-temporal indices.}
    \item{\code{stCop}:}{Object of class \code{"stCopula"} modelling the spatio-temporal dependence on the first tree.}
    \item{\code{topCop}:}{Object of class \code{"copula"} joining the spatio-temporal tree with its additional covariate copula to a full vine copula.}
    \item{\code{dimension}:}{Object of class \code{"integer"} denoting the dimension of the overall copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"}, not used.}
    \item{\code{param.names}:}{Object of class \code{"character"}, not used.}
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}, not used.}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}, not used.}
    \item{\code{fullname}:}{Object of class \code{"character"} providing a textual short description of this class.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
No methods defined with class "stCoVarVineCopula" in the signature.
}
\author{
Benedikt Graeler
}

\seealso{\code{\linkS4class{stVineCopula}}}
\examples{
showClass("stCoVarVineCopula")
}
\keyword{classes}
