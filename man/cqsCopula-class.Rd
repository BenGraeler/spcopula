\name{cqsCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{cqsCopula-class}
\alias{dduCopula,matrix,cqsCopula-method}
\alias{dduCopula,numeric,cqsCopula-method}
\alias{ddvCopula,matrix,cqsCopula-method}
\alias{ddvCopula,numeric,cqsCopula-method}
\alias{fitCopula,cqsCopula-method}
\alias{invdduCopula,numeric,cqsCopula,numeric-method}
\alias{invddvCopula,numeric,cqsCopula,numeric-method}

\title{Class \code{"cqsCopula"}}
\description{A class representing a symmetric copula family with cubic and quadratic sections (Nelsen 2006).}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cqsCopula", ...)} or by the function \code{\link{cqsCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{Object of class \code{"integer"}: defining the dimension of the copula.}
    \item{\code{parameters}:}{Object of class \code{"numeric"}: a vector of length 2 providing parameters a and b.}
    \item{\code{param.names}:}{Object of class \code{"character"}: the names of the parameters ("a", "b") }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} giving the lower bounds of the parameters. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} giving the upper bounds of the parameters. }
    \item{\code{fullname}:}{Object of class \code{"character"} the descriptive full name of the copula }
    \item{\code{fixed}:}{Object of class \code{"character"}: denoting whether parameter "a" or "b" shall be kept fixed an empty character otherwise.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "matrix", copula = "cqsCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "cqsCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "cqsCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "cqsCopula")}: ... }
    \item{fitCopula}{\code{signature(copula = "cqsCopula")}: ... }
    \item{invdduCopula}{\code{signature(u = "numeric", copula = "cqsCopula", y = "numeric")}: ... }
    \item{invddvCopula}{\code{signature(v = "numeric", copula = "cqsCopula", y = "numeric")}: ... }
   }
}
\references{
Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer
}
\author{
Benedikt Graeler
}

\seealso{ See also \code{\linkS4class{asCopula}} and \code{\linkS4class{copula}}.
}
\examples{
showClass("cqsCopula")
}
\keyword{classes}
\keyword{copula}