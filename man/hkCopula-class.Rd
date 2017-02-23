\name{hkCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{hkCopula-class}

\title{Class \code{"hkCopula"}}
\description{
An object representation of Hierarchical Kendall Copulas
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hkCopula", ...)} and via the simplified constructor \code{hkCopula}.
}
\section{Slots}{
  \describe{
    \item{\code{nestingCop}:}{Object of class \code{"copula"} describing the nesting copula at the top.}
    \item{\code{clusterCops}:}{Object of class \code{"list"} with \code{"list"} entries of clustering copulas and their corresponding indices.}
    \item{\code{kenFuns}:}{Object of class \code{"list"} providing the Kendall distribution functions corresponding to the copulas in \code{"clusterCops"}.}
    \item{\code{dimension}:}{Object of class \code{"integer"} giving the dimension of the entire Copula.}
    \item{\code{parameters}:}{Object of class \code{"numeric"} that is empty for compatibility. }
    \item{\code{param.names}:}{Object of class \code{"character"} that is empty for compatibility. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} that is empty for compatibility. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} that is empty for compatibility. }
    \item{\code{fullname}:}{Object of class \code{"character"} describing the copula family. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{parCopula}"}, by class "copula", distance 2.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 3.
}
\section{Methods}{
No methods defined with class "hkCopula" in the signature.
}
\references{
Brechmann, Eike Christian. "Hierarchical Kendall copulas: Properties and inference." Canadian Journal of Statistics 42.1 (2014): 78-108.
}
\author{
Benedikt Graeler
}

\seealso{Vine copulas for a tree structure \code{\linkS4class{vineCopula}} and \code{\link{kendallDistribution}} for details on the Kendall distribution.
}

\examples{
showClass("hkCopula")
}
\keyword{classes}
