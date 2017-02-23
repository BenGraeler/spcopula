\name{stNeighbourhood-class}
\Rdversion{1.1}
\docType{class}
\alias{stNeighbourhood-class}
\alias{[,stNeighbourhood,numeric,ANY,ANY-method}
\alias{names,stNeighbourhood-method}
\alias{show,stNeighbourhood-method}


\title{Class \code{"stNeighbourhood"}}
\description{
An object representing a set of spatio-temporal neighbourhoods including data, spatio-temporal indices and spatio-temporal distances.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stNeighbourhood", ...)} or through the simplified call to \code{\link{getStNeighbours}}.
}

\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"data.frame"} holding the data of spatio-temporal neighbourhoods. }
    \item{\code{distances}:}{Object of class \code{"array"} with the following three dimensions: number of spatio-temporal neighbourhoods, size of each spatio-temporal neighbourhood, 2 (spatial and temporal distance).}
    \item{\code{index}:}{Object of class \code{"array"} with the following three dimensions: number of spatio-temporal neighbourhoods, size of each spatio-temporal neighbourhood, 2 (spatial and temporal index). }
    \item{\code{var}:}{Object of class \code{"character"}; the name of the variable that is estimated.}
    \item{\code{coVar}:}{Object of class \code{"character"}; the name of the covariate.}
    \item{\code{prediction}:}{Object of class \code{"logical"} whether prediction or fitting takes place. }
  }
}

\section{Methods}{
  \describe{
    \item{[[}{\code{signature(x = "neighbourhood", i = "numeric", j = "missing")}: select distinct rows from the neighbourhood slots.}
    \item{names}{\code{signature(x = "neighbourhood")}: retrieve the names of the data slot.}
    \item{show}{\code{signature(object = "neighbourhood")}: print some useful information. }
   }
}

\author{
Benedikt Graeler
}

\seealso{\code{\linkS4class{neighbourhood}}, \code{\link{getStNeighbours}}}

\examples{
showClass("stNeighbourhood")
}
\keyword{classes}
