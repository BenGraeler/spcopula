\name{unitScatter}
\alias{unitScatter}
\alias{univScatter}
\title{
A scatter plot on the unit-square
}
\description{
Plots a sample on the unit-square. If needed, it will transform the sample using the rank order transformation as in \code{\link{rankTransform}}.
}
\usage{
unitScatter(var = NULL, smpl, ...)
}
\arguments{
  \item{var}{
Column IDs or variable names to be used. If not provided, the first two columns will be used.
}
  \item{smpl}{
The data set to be used.
}
\item{...}{
Passed on to \code{\link{xyplot}}.
}
}
\value{
A lattice plot.
}
\author{
Benedikt Graeler
}
\examples{
data("loss")
rt_loss <- rankTransform(loss[,1:2])
unitScatter(smpl=rt_loss)
}
\keyword{ hplot }

