######################################################
##                                                  ##
## a symmetric copula with cubic quadratic sections ##
##                                                  ##
######################################################
# (see Example 3.16 in: Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer)

validCqsCopula <- function(object) {
  if (object@dimension != 2)
    return("Only copulas with cubic quadratic sections of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param > upper | param < lower))
    return("Parameter value out of bound")
  if (object@fixed != ""){
    if(!("a" %in% object@fixed | "b" %in% object@fixed))
      return("The slot fixed may only refer to \"a\" or \"b\".")
    if ("a" %in% object@fixed & "b" %in% object@fixed)
      return("Only one of the parameters may be kept fixed.")
  }
  else return (TRUE)
}

setClass("cqsCopula",
         representation = representation("copula",fixed="character"),
         validity = validCqsCopula,
         contains = list("copula")
)

# constructor 
cqsCopula <- function (param=c(0,0), fixed="") {
  new("cqsCopula", dimension = as.integer(2), parameters = param, 
      param.names = c("a", "b"), param.lowbnd = c(limA(param[2]),-1),
      param.upbnd = c(1, 1), 
      fullname = "copula family with cubic quadratic sections", fixed=fixed)
}

## printing
setMethod("describeCop", c("cqsCopula", "character"),
          function(x, kind = c("short", "very short", "long"), prefix = "", ...) {
            kind <- match.arg(kind)
            if(kind == "very short") # e.g. for show() which has more parts
              return(paste0(prefix, "CQS copula"))
            
            name <- "cubic-quadratic sections"
            d <- dim(x)
            ch <- paste0(prefix, name, " copula, dim. d = ", d)
            switch(kind <- match.arg(kind),
                   short = ch,
                   long = paste0(ch, "\n", prefix, " param.: ",
                                 capture.output(str(x@parameters,
                                                    give.head=FALSE))),
                   stop("invalid 'kind': ", kind))
          })

## density ##
dCQSec <- function (u, copula, log=F) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  
  u1 <- u[, 1]
  u2 <- u[, 2]
  
  if (log)
    return(log(pmax(1-b*(1-2*u2)*(1-2*u1)+(b-a)*(1-u2)*(1-3*u2)*(1-u1)*(1-3*u1),0)))
  else
    return(pmax(1-b*(1-2*u2)*(1-2*u1)+(b-a)*(1-u2)*(1-3*u2)*(1-u1)*(1-3*u1),0))
}

setMethod("dCopula", signature("numeric", "cqsCopula"),
          function(u, copula, log) {
            dCQSec(matrix(u,ncol=copula@dimension), copula, log)
          })
setMethod("dCopula", signature("matrix", "cqsCopula"), dCQSec)

## jcdf ##
pCQSec <- function (u, copula) {
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
return(u1*u2*(1- b*(1-u1)*(1-u2) + (b-a)*(1-u2)^2*(1-u1)^2))
}

setMethod("pCopula", signature("numeric", "cqsCopula"),
          function(u, copula, ...) {
            pCQSec(matrix(u,ncol=copula@dimension), copula)
          })

setMethod("pCopula", signature("matrix","cqsCopula"), pCQSec)

## partial derivatives ##

# solves a*x^3 + b*x^2 + c*x + d = 0
solveCubicEq <- function(a,b,c,d){
eps <- .Machine$double.eps

  # using the reduced equation z^3 + 3 * p * z + q = 0 with:
  p <- 3*a*c-b^2
  q <- 2*b^3-9*a*b*c+27*a^2*d
  D <- q^2+4*p^3

  z <- matrix(NA,nrow=length(D),ncol=3)

  ind <- abs(D) <= eps
  if(any(ind)){
    z[ind,1] <- 0.5*(-4*q[ind])^(1/3)
    z[ind,2] <- -z[ind,1]
  }

  ind <- D > eps
  if(any(ind)){
    cubeRad <- -4*q[ind]+4*sqrt(D[ind])
    r1 <- sign(cubeRad)*abs(cubeRad)^(1/3)
    cubeRad <- -4*q[ind]-4*sqrt(D[ind])
    r2 <- sign(cubeRad)*abs(cubeRad)^(1/3)
    z[ind,1] <- 0.5*(r1+r2)
  }

  ind <- D < eps
  if(any(ind)){
    phi <- acos(-q[ind]/(2*sqrt(-p[ind]^3)))
    triple <- NULL
    triple <- sqrt(-p[ind])*2*cos(phi/3)
    triple <- cbind(triple,sqrt(-p[ind])*2*cos((phi+2*pi)/3))
    triple <- cbind(triple,sqrt(-p[ind])*2*cos((phi+4*pi)/3))
    z[ind,] <- triple
  }

return((z-b)/(3*a))
}

## partial derivative ddu pCQSec

dduCQSec <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  u1 <- u[, 1]
  u2 <- u[, 2]

  return(u2-b*(u2-u2^2-2*u1*u2+2*u1*u2^2)+(b-a)*(u2-4*u1*u2+3*u1^2*u2-2*u2^2+8*u1*u2^2-6*u1^2*u2^2+u2^3-4*u1*u2^3+3*u1^2*u2^3))
}

setMethod("dduCopula", signature("numeric","cqsCopula"),
          function(u, copula, ...) {
            dduCQSec(matrix(u,ncol=copula@dimension), copula)
          }) 
setMethod("dduCopula", signature("matrix","cqsCopula"), dduCQSec)

## inverse partial derivative ddu

## inverse partial derivative ddu
# seems to be accurate (1.4e-05 is the max out of 1000 random CQSec-copulas for 1000 random pairs (u,v) each.)
invdduCQSec <- function (u, copula, y) {
  stopifnot(length(u) == length(y))
  
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  if (a != b) { # the cubic case
    # solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    usq <- u^2
    c3 <- (b-a)*(1-4*u+3*usq)
    c2 <- (b-a)*(-2+8*u-6*u^2)-b*(-1+2*u)
    c1 <- (b-a)*(1-4*u+3*u^2)-b*(1-2*u)+1
    c0 <- -y
  
    v <- solveCubicEq(c3,c2,c1,c0)
    
    filter <- function(vec){
      vec <- vec[!is.na(vec)]
      return(vec[vec >= 0 & vec <= 1])
    }
  
    return(apply(v,1,filter))
  }
  if(a==0) # and b==0 obvioulsy as well: the independent case
    return(y)
  
  # the qudratic cases remain
  v <- y
  uR <- u[u != 0.5]
  v[u != 0.5] <- (-sqrt((-2*b*uR+b-1)^2-4*y[u != 0.5]*(2*b*uR-b))+2*b*uR-b+1)/(2*b*(2*uR-1))
  
  return(v)
}

setMethod("invdduCopula", signature("numeric","cqsCopula","numeric"), invdduCQSec)

## partial derivative ddv

ddvCQSec <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)

  u1 <- u[, 1]
  u2 <- u[, 2]

  return(u1-b*(u1-2*u1*u2-u1^2+2*u1^2*u2)+(b-a)*(u1-2*u1^2+u1^3-4*u1*u2+8*u1^2*u2-4*u1^3*u2+3*u1*u2^2-6*u1^2*u2^2+3*u1^3*u2^2))
}

setMethod("ddvCopula", signature("numeric","cqsCopula"),
          function(u, copula, ...) {
            ddvCQSec(matrix(u,ncol=copula@dimension), copula)
          })
setMethod("ddvCopula", signature("matrix","cqsCopula"), ddvCQSec)

## inverse partial derivative ddv
invddvCQSec <- function (v, copula, y) {
  stopifnot(length(v)==length(y)) 

  a <- copula@parameters[1]
  b <- copula@parameters[2]

  if (a != b) { # the cubic case
    # solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    vsq <- v^2
    c3 <- (b-a)*(1-4*v+3*vsq)
    c2 <- (b-a)*(-2+8*v-6*vsq)-b*(-1+2*v)
    c1 <- (b-a)*(1-4*v+3*vsq)-b*(1-2*v)+1
    c0 <- -y
  
    u <- solveCubicEq(c3,c2,c1,c0)
    
    filter <- function(vec){
      vec <- vec[!is.na(vec)]
      return(vec[vec >= 0 & vec <= 1])
    }
  
    return(apply(u,1,filter))
  }
  if(a==0) # and b==0 obvioulsy as well: the independent case
    return(y)
  
  # the qudratic cases remain
  u <- y
  vR <- v[v != 0.5]
  u[v != 0.5] <- (-sqrt((-2*b*vR+b-1)^2-4*y[v != 0.5]*(2*b*vR-b))+2*b*vR-b+1)/(2*b*(2*vR-1))
  
  return(u)
}

setMethod("invddvCopula", signature("numeric","cqsCopula","numeric"), invddvCQSec)

## random number generator

rCQSec <- function (n, copula) {
  u <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
    
  res <- cbind(u, invdduCQSec(u, copula, y))
  colnames(res) <- c("u","v")
    
  return(res)
}

setMethod("rCopula", signature("numeric","cqsCopula"), rCQSec)
## fitment

fitCopula.cqs <- function (copula, data, method = "ml", start=c(0,0), 
                           lower=c(-3,-1), upper=c(1,1), 
                           optim.method="L-BFGS-B", optim.control=list(),
                           estimate.variance = FALSE) {
  fit <- switch(method,
                ml=fitCQSec.ml(copula, data, start, lower, upper, optim.control, optim.method),
                itau=fitCQSec.itau(copula, data, estimate.variance),
                irho=fitCQSec.irho(copula, data, estimate.variance),
                stop("Implemented methods for copulas in the spcopula package are: ml, itau, and irho."))
  return(fit)
}

setMethod("fitCopula", signature("cqsCopula"), fitCopula.cqs)

## Fits the copula with cubic and quadratic sections acoording to a measure of association.
## It performs a maximum likelihood evaluation over all possible pairs of 
## parameters a and b generating a copula with the given mesaure of 
## association.
#
# moa
#  measure of association, according to method
# data
#  the bivariate data set as a 2-column matrix within the unitsquare
# method
#  one of kendall or spearman according to the calculation of moa

fitCQSec.itau <- function(copula, data, estimate.variance=FALSE, tau=NULL) {
  if(is.null(tau))
    tau <- TauMatrix(data)[1,2]
  if(copula@fixed=="a")
    esti <- c(copula@parameters[1], iTauCQSec.a(copula@parameters[1], tau))
  if(copula@fixed=="b")
    esti <- c(iTauCQSec.b(copula@parameters[2], tau),copula@parameters[2])
  else
    esti <- fitCQSec.moa(tau, data, method="itau")
  
  copula <- cqsCopula(esti, fixed=copula@fixed)

  new("fitCopula", estimate = esti, var.est = matrix(NA), 
      method = "Inversion of Kendall's tau and MLE", 
      loglik = sum(dCopula(data, copula, log=T)),
      fitting.stats=list(convergence = as.integer(NA)), nsample = nrow(data),
      copula=copula)
}


fitCQSec.irho <- function(copula, data, estimate.variance=FALSE, rho=NULL){
  if(is.null(rho))
    rho <- cor(data,method="spearman")[1,2]
  if(copula@fixed=="a")
    esti <- c(copula@parameters[1], iRhoCQSec.a(copula@parameters[1],rho))
  if(copula@fixed=="b")
    esti <- c(iRhoCQSec.b(copula@parameters[2],rho),copula@parameters[2])
  else
    esti <- fitCQSec.moa(rho, data, method="irho")
  
  copula <- cqsCopula(esti, fixed=copula@fixed)

  new("fitCopula", estimate = esti, var.est = matrix(NA), 
      method = "Inversion of Spearman's rho and MLE", 
      loglik = sum(dCopula(data, copula, log=T)),
      fitting.stats=list(convergence = as.integer(NA)), nsample = nrow(data),
      copula=copula)
}

fitCQSec.moa <- function(moa, data, method="itau", tol=.Machine$double.eps^.5) {
  smpl <- as.matrix(data)

  iRho.CQS <- function(p) {
    iRhoCQSec.b(p,moa)
  }
  
  iTau.CQS <- function(p) {
    iTauCQSec.b(p,moa)
  }
  
  iFun <- switch(method, itau=iTau.CQS, irho=iRho.CQS)

  sec <- function (parameters) {
    res <- NULL
    for(param in parameters) {
      res <- rbind(res, -sum(dCQSec(smpl, cqsCopula(c(iFun(param),param)),log=T)))
    }
    
    return(res)
  }

  b <- optimize(sec,c(-1,1), tol=tol)$minimum

  return(c(iFun(b),b))
}

# maximum log-likelihood estimation of a and b using optim

fitCQSec.ml <- function(copula, data, start, lower, upper, optim.control, optim.method) { 
  if(length(start)!=2) stop("Start values need to have same length as parameters.")
  
  if (copula@fixed=="") {
    optFun <- function(param=c(0,0)) {
      if(any(param > 1) | param[2] < -1 | param[1] < limA(param[2])) 
        return(100)
      return(-sum(log( dCQSec(data, cqsCopula(param)) )))
    }
  
    optimized <- optim(par=start, fn=optFun, method = optim.method, 
                       lower=lower, upper=upper, control = optim.control)
  
    return(new("fitCopula", estimate = optimized$par, var.est = matrix(NA),
               method = "Numerical MLE over the full range.",
               loglik = -optimized$value, fitting.stats= optimized,
               nsample = nrow(data), copula=cqsCopula(optimized$par)))
  } else {  
    optFunFixed <- function(p) {
      param <- switch(copula@fixed, a=c(copula@parameters[1],p),
                      b=c(p,copula@parameters[2]))
      if(any(param > 1) | param[2] < -1 | param[1] < limA(param[2])) 
        return(100)
      return(-sum(log( dCQSec(data, cqsCopula(param)) )))
    }
  
    optimized <- optimise(optFunFixed, c(-3,1))
  
    optPar <- switch(copula@fixed, a=c(copula@parameters[1],optimized$minimum),
                     b=c(optimized$minimum,copula@parameters[2]))
  
    return(new("fitCopula", estimate = optimized$minimum, var.est = matrix(NA),
               method = "Numerical MLE over the full range.",
               loglik = -optimized$objective, fitting.stats=list(),
               nsample = nrow(data), copula=cqsCopula(optPar)))
  }
}

####

iTauCQSec.a <- function(a, tau=0) {
  limits <- limB(a)
  min(max(limits[1],0.5*(sqrt(a^2-250*a-1800*tau+5626)+a-75)),limits[2])
}

iTauCQSec.b <- function(b,tau=0) {
  min(max(limA(b),(b^2 + 75*b + 450*tau)/(b - 25)),1)
}

setMethod("iTau",signature="cqsCopula",
          function(copula, tau) {
            switch(copula@fixed,
                   a=c(copula@parameters[1],iTauCQSec.a(copula@parameters[1],tau)), 
                   b=c(iTauCQSec.b(copula@parameters[2],tau),copula@parameters[2]),
                   stop("iTau may only be used for cqsCopula with one parameter fixed."))
            })

####

tauCQSec <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return( (a*b - 25*a - b^2 - 75*b)/450 )
}

setMethod("tau",signature("cqsCopula"),tauCQSec)

####
# find parameter "a" for parameter "b" under a given measure of association "rho" 
# it may return a value exceeding the limit of "a" which may result in an invalid copula.

iRhoCQSec.a <- function(a, rho=0) {
  limits <- limB(a)
  min(max(limits[1],-a/3 - 4*rho),limits[2])
}

iRhoCQSec.b <- function(b, rho=0) {
  min(max(limA(b),-3*b - 12*rho),1)
}

setMethod("iRho",signature="cqsCopula",
          function(copula, rho) {
            switch(copula@fixed,
                   a=function(copula, rho) c(copula@parameters[1],iRhoCQSec.a(copula@parameters[1],rho)), 
                   b=function(copula, rho) c(iRhoCQSec.b(copula@parameters[2],rho),copula@parameters[2]),
                   stop("iRho may only be used for cqsCopula with one parameter fixed."))
            })

####

rhoCQSec <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return( -(a+3*b)/12 )
}

setMethod("rho",signature("cqsCopula"),rhoCQSec)

setMethod("lambda", signature("cqsCopula"), 
          function(copula, ...) c(lower = 0, upper = 0))