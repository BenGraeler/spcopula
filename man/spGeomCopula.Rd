\name{spGeomCopula}
\alias{spGeomCopula}
\title{
Constructor of a the spatial copula class based on geometric means
}
\description{
Generates an instance of a spatial copula class based on geometric means.
}
\usage{
spGeomCopula(components, distances, unit = NULL)
}
\arguments{
  \item{components}{
the list of components to be used in the spatial geometric copula
}
  \item{distances}{
the corresponding distances associated with the components
}
  \item{unit}{
the unit of the distance measure
}
}
\value{
An instance of the spatial copula class based on geometric means.
}
\author{
Benedikt Graeler
}
\examples{
spGeomCopula(list(normalCopula(0.2), gumbelCopula(1.2), claytonCopula(2)),
              c(0, 120, 300),
             "m")
}
