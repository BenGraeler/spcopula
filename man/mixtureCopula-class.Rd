\name{mixtureCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{mixtureCopula-class}
\alias{dduCopula,ANY,mixtureCopula-method}
\alias{ddvCopula,ANY,mixtureCopula-method}
\alias{fitCopula,mixtureCopula,ANY-method}
\alias{invdduCopula,numeric,mixtureCopula,numeric-method}
\alias{invddvCopula,numeric,mixtureCopula,numeric-method}

\title{Class \code{"mixtureCopula"}}
\description{
The \code{mixtureCopula} contains two copulas that are combined in a convex manner to a new copula.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mixtureCopula", ...)} or by calling the constructor \code{\link{mixtureCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{memberCops}:}{Object of class \code{"list"}: holding the copulas used in the convex combination. }
    \item{\code{dimension}:}{Object of class \code{"integer"}: the dimension of both input and the resulting copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"}: the set of parameters for the first and second copula as well as the mixing coefficient.}
    \item{\code{param.names}:}{Object of class \code{"character"}: names of all parameters. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}: lower bounds of all parameters. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}: upper bounds of all parameters. }
    \item{\code{fullname}:}{Object of class \code{"character"}: the name of the copula. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{parCopula}"}, by class "copula", distance 2.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 3.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "ANY", copula = "mixtureCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "ANY", copula = "mixtureCopula")}: ... }
    \item{fitCopula}{\code{signature(copula = "mixtureCopula")}: ... }
    \item{invdduCopula}{\code{signature(u = "numeric", copula = "mixtureCopula", y = "numeric")}: ... }
    \item{invddvCopula}{\code{signature(v = "numeric", copula = "mixtureCopula", y = "numeric")}: ... }
	 }
}

\author{
Benedikt Graeler}

\examples{
showClass("mixtureCopula")
}
\keyword{classes}
