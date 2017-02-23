\name{spCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{spCopula-class}
\alias{dduCopula,matrix,spCopula-method}
\alias{ddvCopula,matrix,spCopula-method}
\alias{dduCopula,numeric,spCopula-method}
\alias{ddvCopula,numeric,spCopula-method}

\title{Class \code{"spCopula"}}
\description{
A class representing a bivariate spatial copula.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("spCopula", ...)} or by the function \code{\link{spCopula}}.
}

\section{Slots}{
  \describe{
    \item{\code{components}:}{Object of class \code{"list"} A list of copula components per lag.}
    \item{\code{distances}:}{Object of class \code{"numeric"} The mean distances of the lag classes. }
    \item{\code{calibMoa}:}{Object of class \code{"function"} A function providing the copulas' parameter based on a measure of association.}
    \item{\code{unit}:}{Object of class \code{"character"} The spatial unit. }
    \item{\code{dimension}:}{Object of class \code{"integer"} The dimension of the copula (2). }
    \item{\code{parameters}:}{Object of class \code{"numeric"} The parameters of the components. }
    \item{\code{param.names}:}{Object of class \code{"character"} Names of the above parameters. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} The lower bounds of the above parameters. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} The upper bounds of the above parameters. }
    \item{\code{fullname}:}{Object of class \code{"character"}A phrase describing the spatial copula.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "numeric", copula = "spCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "spCopula")}: ... }
	 }
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}

\examples{
# data from demo(spcopula_estimation)
data("spCopDemo") 

calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")
                  
dCopula(u=matrix(c(.3,.3,.7,.7),ncol=2),spCop,h=c(200,400))
pCopula(u=matrix(c(.3,.3,.7,.7),ncol=2),spCop,h=c(200,400))
}
\keyword{classes}
\keyword{spcopula}
\keyword{copula}
