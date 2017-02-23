\name{calcSpTreeDists}
\alias{calcSpTreeDists}

\title{
Calculating the conditional spatial distances
}
\description{
All spatial distances of all involved neighbours are calculated according to multiple spatial copulas forming a C-Vine.
}
\usage{calcSpTreeDists(neigh, dataLocs, n.trees)}

\arguments{
  \item{neigh}{the initial \code{\linkS4class{neighbourhood}}}
  \item{dataLocs}{the initial \code{\linkS4class{Spatial}} data set.}
  \item{n.trees}{the number of spatial trees}
}

\value{
A \code{list} of the spatial distance of the pairs involved in the conditional data pairs.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{neighbourhood}}
}

\keyword{ spatial }
\keyword{ multivariate }
\keyword{distribution}
