\name{spGeomCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{spGeomCopula-class}

\title{Class \code{"spGeomCopula"}}
\description{
A spatial copula class based on geometric means
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("spGeomCopula", ...)} or \code{spGeomCopula(...)}
}
\section{Slots}{
  \describe{
    \item{\code{components}:}{Object of class \code{"list"} ~~ }
    \item{\code{distances}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{unit}:}{Object of class \code{"character"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{param.names}:}{Object of class \code{"character"} ~~ }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fullname}:}{Object of class \code{"character"} ~~ }
    \item{\code{dimension}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{dimCopula}"}, by class "copula", distance 2.
Class \code{"\linkS4class{parCopula}"}, by class "copula", distance 2.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 3.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "matrix", copula = "spGeomCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "spGeomCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "spGeomCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "spGeomCopula")}: ... }
    \item{invdduCopula}{\code{signature(u = "numeric", copula = "spGeomCopula", y = "ANY")}: ... }
    \item{invddvCopula}{\code{signature(v = "numeric", copula = "spGeomCopula", y = "ANY")}: ... }
	 }
}
\author{
Benedikt Graeler
}
\examples{
showClass("spGeomCopula")
}
\keyword{classes}
