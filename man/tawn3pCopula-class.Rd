\name{tawn3pCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{tawn3pCopula-class}
\alias{A,tawn3pCopula-method}
\alias{dAdu,tawn3pCopula-method}
\alias{dCopula,matrix,tawn3pCopula-method}
\alias{dCopula,numeric,tawn3pCopula-method}
\alias{fitCopula,tawn3pCopula-method}
\alias{pCopula,matrix,tawn3pCopula-method}
\alias{pCopula,numeric,tawn3pCopula-method}
\alias{rCopula,numeric,tawn3pCopula-method}
\alias{dduCopula,ANY,tawn3pCopula-method}
\alias{ddvCopula,ANY,tawn3pCopula-method}

\title{Class \code{"tawn3pCopula"}}
\description{
S4-Class representing the three parameter Tawn Copula family. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tawn3pCopula", ...)} or more likely with the help of the constructor \code{\link{tawn3pCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"} defining the cdf and pdf of the copula.}
    \item{\code{dimension}:}{Object of class \code{"integer"} defining the dimension (currently fixed at 2.) }
    \item{\code{parameters}:}{Object of class \code{"numeric"} providing the three parameters. }
    \item{\code{param.names}:}{Object of class \code{"character"} giving the names of the parameters. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} defining the lower bounds of the three parameters. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} defining the upper bounds of the three parameters. }
    \item{\code{fullname}:}{Object of class \code{"character"} providing a descriptive name for the copula family.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{evCopula}"}, directly.
Class \code{"\linkS4class{copula}"}, by class "evCopula", distance 2.
Class \code{"\linkS4class{Copula}"}, by class "evCopula", distance 3.
}
\section{Methods}{
  \describe{
    \item{A}{\code{signature(copula = "tawn3pCopula")}: ... }
    \item{dAdu}{\code{signature(copula = "tawn3pCopula")}: ... }
    \item{dCopula}{\code{signature(u = "matrix", copula = "tawn3pCopula")}: ... }
    \item{dCopula}{\code{signature(u = "numeric", copula = "tawn3pCopula")}: ... }
    \item{fitCopula}{\code{signature(copula = "tawn3pCopula")}: ... }
    \item{pCopula}{\code{signature(u = "matrix", copula = "tawn3pCopula")}: ... }
    \item{pCopula}{\code{signature(u = "numeric", copula = "tawn3pCopula")}: ... }
    \item{rCopula}{\code{signature(n = "numeric", copula = "tawn3pCopula")}: ... }
	 }
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{tawnCopula}} for a symmetric one parameter, \code{\linkS4class{tawnT1Copula}} and \code{\linkS4class{tawnT2Copula}} for asymmetric two-parameter sub-families of the Tawn copula family. 
}
\examples{
showClass("tawn3pCopula")
}
\keyword{classes}
