\name{asCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{asCopula-class}
\alias{dduCopula,matrix,asCopula-method}
\alias{dduCopula,numeric,asCopula-method}
\alias{ddvCopula,matrix,asCopula-method}
\alias{ddvCopula,numeric,asCopula-method}
\alias{fitCopula,asCopula,ANY-method}
\alias{invdduCopula,numeric,asCopula,numeric-method}
\alias{invddvCopula,numeric,asCopula,numeric-method}

\title{Class \code{"asCopula"}}
\description{A class representing an asymmetric copula family with cubic and quadratic sections (Nelsen 2006).}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("asCopula", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{Object of class \code{"integer"}: defining the dimension of the copula.}
    \item{\code{parameters}:}{Object of class \code{"numeric"}: a vector of length 2 providing parameters a and b.}
    \item{\code{param.names}:}{Object of class \code{"character"}: the names of the parameters ("a", "b") }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} giving the lower bounds of the parameters. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} giving the upper bounds of the parameters. }
    \item{\code{fullname}:}{Object of class \code{"character"} the descriptive full name of the copula }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "matrix", copula = "asCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "asCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "asCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "asCopula")}: ... }
    \item{fitCopula}{\code{signature(copula = "asCopula")}: ... }
    \item{invdduCopula}{\code{signature(u = "numeric", copula = "asCopula", y = "numeric")}: ... }
    \item{invddvCopula}{\code{signature(v = "numeric", copula = "asCopula", y = "numeric")}: ... }
	 }
}
\references{
Example 3.16 in: Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer
}
\author{
Benedikt Graeler
}

\seealso{ See also \code{\linkS4class{cqsCopula}} and \code{\linkS4class{copula}}.
}
\examples{
showClass("asCopula")
}
\keyword{classes}
\keyword{asymmetric copula}
\keyword{copula}