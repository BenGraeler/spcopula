\name{hkCopula}
\alias{hkCopula}

\title{Constructor of a hierarchical Kendall copula}
\description{A copula object from the copula family of hierarchical Kendall copulas. Cluster copulas are linked via their Kendall distribution functions in a nesting copula.}

\usage{
hkCopula(nestingCop, clusterCops, kenFuns=NULL)
}

\arguments{
  \item{nestingCop}{The nesting copula that combines the clusters. Any copula is allowed.}
  \item{clusterCops}{A list of lists defining the cluster copulas and the indices used. Any copula is allowed, but the evaluation of the Kendall distribution function is time consuming for many families. Archimedean copulas are preferred from a computation time point of view.}
  \item{kenFuns}{A list of corresponding Kendall distribution functions of the same length as \code{clusterCops}.}
  
}

\value{
A \code{\linkS4class{hkCopula}} object that behaves (almost) like any other \code{\linkS4class{copula}} from the \code{copula} package.
}
\references{
Brechmann, Eike Christian. "Hierarchical Kendall copulas: Properties and inference." Canadian Journal of Statistics 42.1 (2014): 78-108.
}
\author{
Benedikt Graeler
}

\seealso{ 
See also \code{\linkS4class{hkCopula}}, \code{\linkS4class{vineCopula}}.
}
\examples{
hkCop <- hkCopula(gumbelCopula(5),list(list(frankCopula(3), c(1,2)),
                                       list(claytonCopula(3), c(3,4))))
hkCop

dCopula(matrix(runif(400), 100), hkCop)
}

\keyword{asymmetric}
\keyword{copula}
\keyword{cubic quadratic sections}
