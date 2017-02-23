\name{dropStTree}
\alias{dropStTree}

\title{
Drops a spatio-temporal tree
}
\description{
A spatio-temporal tree is dropped in order to fit/evaluate the following copula.
}
\usage{dropStTree(stNeigh, dataLocs, stCop)}

\arguments{
  \item{stNeigh}{the current spatio-temporal \code{\linkS4class{stNeighbourhood}}}
  \item{dataLocs}{the data locations (the same as used for the generation of the spatio-temporal neighbourhood).}
  \item{stCop}{the current spatio-temporal copula performing the conditioning}
}

\value{
A conditioned spatio-temporal \code{\linkS4class{stNeighbourhood}} of dimension 1 less.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{stNeighbourhood}} and \code{\link{dropSpTree}}
}

\keyword{ spatio-temporal }