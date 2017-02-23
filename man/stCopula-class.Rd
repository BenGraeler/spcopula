\name{stCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{stCopula-class}
\alias{show,stCopula-method}
\alias{dduCopula,numeric,stCopula-method}
\alias{dduCopula,matrix,stCopula-method}
\alias{ddvCopula,numeric,stCopula-method}
\alias{ddvCopula,matrix,stCopula-method}

\title{Class \code{"stCopula"}}
\description{
A class representing bivariate spatio-temporal copulas.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stCopula", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{spCopList}:}{Object of class \code{"list"} the list of bivariate spatial copulas \code{\linkS4class{spCopula}}, one per time lag }
    \item{\code{tlags}:}{Object of class \code{"numeric"} the time lags }
    \item{\code{tres}:}{Object of class \code{"character"} the temporal resolution }
    \item{\code{dimension}:}{Object of class \code{"integer"} the dimension: 2 }
    \item{\code{parameters}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.names}:}{Object of class \code{"character"} empty }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{fullname}:}{Object of class \code{"character"} a descriptive name }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "stCopula")}: ... }
	 }
}
\references{
Graeler, B. & E. Pebesma (2012): Modelling Dependence in Space and Time with Vine Copulas. 
Presentation at: Geostats 2012, Oslo, Norway, 11-15 June 2012
}
\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{spCopula}}
}
\examples{
showClass("stCopula")
}
\keyword{classes}
