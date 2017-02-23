########################################
##                                    ##
## an empirical copula representation ##
##                                    ##
########################################

# validity
validEmpCopula <- function(object) {
  if(ncol(object@sample) != object@dimension)
    return("Dimension of the copula and the sample do not match.")
  else
    return(TRUE)
}

# class definition
setClass("empiricalCopula",
         representation = representation("copula", sample="matrix"),
         validity = validEmpCopula,
         contains = list("copula")
)

# constructor
empiricalCopula <- function (sample=NULL, copula) {
  if(is.null(sample) && missing(copula))
    stop("At least one parameter of copula or sample must be provided.")
  
  if(is.null(sample))
    return(genEmpCop(copula))
  
  if(missing(copula))
    return(new("empiricalCopula", dimension = as.integer(ncol(sample)), 
               parameters = as.numeric(NA), param.names = "unknown", 
               param.lowbnd = as.numeric(NA), param.upbnd = as.numeric(NA), 
               fullname = "Unkown empirical copula based on a sample.",
               sample=sample))
  
  new("empiricalCopula", dimension = copula@dimension, 
      parameters = copula@parameters, param.names = copula@param.names, 
      param.lowbnd = copula@param.lowbnd, param.upbnd = copula@param.upbnd, 
      fullname = paste("Empirical copula derived from", describeCop(copula, "very short")),
      sample=sample)
}

# simplified constructor
genEmpCop <- function(copula, sample.size=1e5) {
  cat("Note: the copula will be empirically represented by a sample of size:",
      sample.size, "\n")
  empiricalCopula(rCopula(sample.size, copula), copula)
}

# printing
setMethod("describeCop", c("empiricalCopula", "character"),
          function(x, kind = c("short", "very short", "long"), prefix = "", ...) {
            kind <- match.arg(kind)
            name <- "empirical"
            if(kind == "very short") # e.g. for show() which has more parts
              return(paste0(prefix, name, " copula"))
            ## else
            d <- dim(x)
            ch <- paste0(prefix, name, " copula, dim. d = ", d)
            switch(kind <- match.arg(kind),
                   short = ch,
                   long = paste0(ch, "\n", prefix, " param.: ",
                                 capture.output(str(x@parameters,
                                                    give.head=FALSE))),
                   stop("invalid 'kind': ", kind))
          })

## density, not yet needed and hence not implemented ##

## jcdf ##
# from package copula
pempCop.C <- function(u, copula) {
  F.n(u, copula@sample)
}

setMethod("pCopula", signature("numeric", "empiricalCopula"),
          function(u, copula, ...) {
            pempCop.C(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix", "empiricalCopula"), pempCop.C)


tauempCop <- function(copula){
  TauMatrix(copula@sample)[1,2]
}

setMethod("tau",signature("empiricalCopula"), tauempCop)


rhoempCop <- function(copula){
  cor(copula@sample,method="spearman")
}

setMethod("rho",signature("empiricalCopula"), rhoempCop)

setMethod("lambda", signature("empiricalCopula"), 
          function(copula, ...) stop("No evaluation possible, try to plot 'empBivJointDepFun' for a visual assessment."))

##################################
##                              ##
## an empirical survival copula ##
##                              ##
##################################

# constructor
empSurCopula <- function (sample=NULL, copula) {
  if(is.null(sample) && missing(copula))
    stop("At least one parameter of copula or sample must be provided.")
  
  if(is.null(sample))
    return(genEmpSurCop(copula))
  
  if(missing(copula))
    return(new("empSurCopula", dimension = as.integer(ncol(sample)), 
               parameters = as.numeric(NA), param.names = "unknown", 
               param.lowbnd = as.numeric(NA), param.upbnd = as.numeric(NA), 
               fullname = "Unkown empirical survival copula based on a sample.",
               sample=sample))
  
  new("empSurCopula", dimension = copula@dimension, 
      parameters = copula@parameters, param.names = copula@param.names, 
      param.lowbnd = copula@param.lowbnd, param.upbnd = copula@param.upbnd, 
      fullname = paste("Empirical survival copula derived from", describeCop(copula, "very short")),
      sample=sample)
}

# simplified constructor
genEmpSurCop <- function(copula, sample.size=1e5) {
  cat("Note: the survival copula will be empirically represented by a sample of size:",
      sample.size, "\n")
  empSurCopula(1 - rCopula(sample.size, copula), copula)
}

## jcdf ##
# from package copula # 3D
pempSurCop.C <- function(u, copula) {
  if (copula@dimension==2)
    return(apply(u,1,sum) - 1 + F.n(1-u, copula@sample))
  
  if (copula@dimension==3) 
    return(apply(u,1,sum) - 2 + F.n(cbind(1-u[, 1:2, drop=F],1),
                                    copula@sample) + F.n(cbind(1-u[,1, drop=F], 1, 1-u[,3, drop=F]),
                                                         copula@sample) + F.n(cbind(1-u[,2:3, drop=F], 1),
                                                                              copula@sample) - F.n(1-u, copula@sample))
  
  stop("The empirical survival copula is only implemented for 2 and 3 dimensions.")
}

setMethod("pCopula", signature("numeric", "empSurCopula"),
          function(u, copula, ...) {
            pempSurCop.C(matrix(u, ncol=copula@dimension), copula)
          })
setMethod("pCopula", signature("matrix", "empSurCopula"), pempSurCop.C)


# tauempCop <- function(copula){*-
#   TauMatrix(copula@sample)[1,2]
# }
# 
# setMethod("tau",signature("empiricalCopula"), tauempCop)
# 
# 
# rhoempCop <- function(copula){
#   cor(copula@sample,method="spearman")
# }
# 
# setMethod("rho",signature("empiricalCopula"), rhoempCop)
# 
# setMethod("lambda", signature("empiricalCopula"), 
#           function(copula, ...) stop("No evaluation possible, try to plot 'empBivJointDepFun' for a visual assessment."))
# 


# Vine Copula - empirical evaluation
## jcdf ##
pvineCopula <- function(u, copula) {
  empCop <- genEmpCop(copula, 1e5)
  
  return(pCopula(u, empCop))
}

setMethod("pCopula", signature("numeric","vineCopula"), 
          function(u, copula) {
            pvineCopula(matrix(u, ncol=copula@dimension), copula)
          })
setMethod("pCopula", signature("data.frame","vineCopula"), 
          function(u, copula) pvineCopula(as.matrix(u), copula))
setMethod("pCopula", signature("matrix","vineCopula"), pvineCopula)