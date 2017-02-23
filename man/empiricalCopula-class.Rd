\name{empiricalCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{empiricalCopula-class}

\title{Class \code{"empiricalCopula"}}
\description{
A class representing an empirical version of a given copula.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("empiricalCopula", ...)}, 
through the constructor \code{\link{empiricalCopula}} or the simplified constructor \code{\link{genEmpCop}}.

}
\section{Slots}{
  \describe{
    \item{\code{sample}:}{Object of class \code{"matrix"} representing the sample. }
    \item{\code{dimension}:}{Object of class \code{"integer"} giving the dimension. }
    \item{\code{parameters}:}{Object of class \code{"numeric"} giving the parameters of the underlying copula family if present, \"NA\" otherwise.}
    \item{\code{param.names}:}{Object of class \code{"character"} giving the parameter names of the underlying copula family if present, \"unknown\" otherwise. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}  giving the lower bound of the underlying copula family if present, \"NA\" otherwise.}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}  giving the upper bound of the underlying copula family if present, \"NA\" otherwise.}
    \item{\code{fullname}:}{Object of class \code{"character"}  giving a descriptive name including the underlying copula family's name if present. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
No additional methods defined with class \code{empiricalCopula} in the signature, but \code{\link{pCopula}} works as expected.
}
\author{
Benedikt Graeler}
\note{
Its implementation of \code{\link{pCopula}} is based on C-code from \code{\link{copula-package}}.
}


\seealso{\code{\link{genEmpCop}}}
\examples{
showClass("empiricalCopula")
}
\keyword{classes}
