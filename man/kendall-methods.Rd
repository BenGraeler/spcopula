\name{kendall-methods}
\docType{methods}
\alias{kendall}
\alias{kendall-methods}
\alias{kendall,ANY,ANY-method}
\alias{kendall,numeric,claytonCopula-method}
\alias{kendall,numeric,frankCopula-method}
\alias{kendall,numeric,gumbelCopula-method}
\alias{kendall,numeric,BB1Copula-method}
\alias{kendall,numeric,BB6Copula-method}
\alias{kendall,numeric,BB7Copula-method}
\alias{kendall,numeric,BB8Copula-method}
\alias{kendall,numeric,joeBiCopula-method}
\alias{kendall,numeric,joeCopula-method}

\title{ Methods for Function \code{kendall} }
\description{
  Methods for function \code{kendall} in case of Archimedian copulas up to 4 dimensions
}
\section{Methods}{
\describe{

\item{\code{signature(t = "ANY", copula = "ANY")}}{nothing}

\item{\code{signature(t = "numeric", copula = "claytonCopula")}}{Kendall distribtuion for Clayton family.}

\item{\code{signature(t = "numeric", copula = "frankCopula")}}{Kendall distribtuion for Frank family.}

\item{\code{signature(t = "numeric", copula = "gumbelCopula")}}{Kendall distribtuion for Gumbel family.}
}}
\keyword{methods}
