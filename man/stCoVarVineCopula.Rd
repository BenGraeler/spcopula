\name{stCoVarVineCopula}
\alias{stCoVarVineCopula}

\title{
Constructor for \code{\linkS4class{stCoVarVineCopula}}
}
\description{
This function provides a more comfortable way of defining a \code{\linkS4class{stCoVarVineCopula}}.
}
\usage{
stCoVarVineCopula(coVarCop, stCop, topCop)
}

\arguments{
  \item{coVarCop}{
A function returning a returning a bivariate copula object for absolute pairs of spatio-temporal indices
}
  \item{stCop}{
Object of class \code{"stCopula"} modelling the spatio-temporal dependence on the first tree
}
  \item{topCop}{
Object of class \code{"copula"} joining the spatio-temporal tree with its additional covariate copula to a full vine copula.
}
}
\details{
For a spatio-temporal random field Z with covariate Y a c-vine is assumed with data sorted as (z_0, y_0, z_1, .., z_n).
}

\value{An object of class \code{\linkS4class{stCoVarVineCopula}}.}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{stVineCopula}}, \code{\linkS4class{stCoVarVineCopula}}
}

\examples{
# a spatio-temporal C-vine copula (with independent copulas in the upper vine)
library("VineCopula")
spCopT0 <- spCopula(components=list(claytonCopula(8), claytonCopula(4), 
                                    claytonCopula(2), claytonCopula(1),
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300,400,500),
                    unit="km")
spCopT1 <- spCopula(components=list(claytonCopula(4), claytonCopula(2), 
                                    claytonCopula(1), claytonCopula(0.5),
                                    indepCopula()),
                    distances=c(0,100,200,300,400),
                    unit="km")
spCopT2 <- spCopula(components=list(claytonCopula(2), claytonCopula(1), 
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300),
                    unit="km")

stCop <- stCopula(components=list(spCopT0, spCopT1, spCopT2),
                  tlags=-(0:2))

# only a constant copula ius used for the covariate
stCoVarVineCopula(function(x) gumbelCopula(2), stCop, vineCopula(9L))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
