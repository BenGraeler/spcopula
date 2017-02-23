\name{dependencePlot}
\alias{dependencePlot}
\title{
Kernel smoothed scatter plot
}
\description{
Plots a kernel smoothed scatter plot of the provided rank-transformed sample. The work is done by the function \code{\link{panel.smoothScatter}}.
}
\usage{
dependencePlot(var = NULL, smpl, bandwidth = 0.075,  main="Strength of dependence",
               transformation = function(x) x, margin=NULL, ...)
}

\arguments{
  \item{var}{Column IDs or variable names to be used. If not provided, the first two columns will be used.}
  \item{smpl}{a matrix (two-columns at least) holding the data}
  \item{bandwidth}{the bandwidth passed to the smoothing kernel}
  \item{main}{the title of the plot}
  \item{transformation}{a transformation passed to the kernel}
  \item{margin}{a quantile function to back transform the uniform margins to any other desired marginal distribution (typically \code{\link{qnorm}}). Note that the \code{bandwidth} parameter might need adjustment.}
  \item{\dots}{passed on to the function \code{\link{panel.smoothScatter}}}
}
\details{
see \code{\link{panel.smoothScatter}}
}
\value{
A plot.
}
\author{
Benedikt Graeler
}
\note{
This is simple wrapper to \code{\link{smoothScatter}}.
}

\seealso{
\code{\link{smoothScatter}}
}
\examples{
\dontrun{dependencePlot(smpl=rCopula(500,asCopula(c(-1,1))))}
}
\keyword{plot}
