##########################
##                      ##
## an asymmetric copula ##
##                      ##
##########################
# (see Example 3.16 in: Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer)

validAsCopula = function(object) {
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
  else return (TRUE)
}

# the lower bound of the parameter a dependening on the parameter b
limA <- function (b) {
  stopifnot(abs(b) <= 1)
  0.5*(-sqrt(-3*b^2+6*b+9)+b-3)
}

# the lower and upper bound of the parameter b dependening on the parameter a
limB <- function (a) {
  stopifnot(a <=1 & a >= -3)
  if(a>-2)
    return(c(-1,1))
  pmax(pmin(0.5*(c(-1,1)*(sqrt(3)*sqrt(-a^2-2*a+3))+a+3),1),-1)
}

setClass("asCopula",
         representation = representation("copula"),
         validity = validAsCopula,
         contains = list("copula")
)

# constructor
asCopula <- function (param=c(0,0)) {
  val <- new("asCopula", dimension = as.integer(2), parameters = param, 
             param.names = c("a", "b"), param.lowbnd = c(limA(param[2]), -1),
             param.upbnd = c(1, 1), fullname = "asymmetric copula family with cubic and quadratic sections")
  return(val)
}

## printing
setMethod("describeCop", c("asCopula", "character"),
          function(x, kind = c("short", "very short", "long"), prefix = "", ...) {
            kind <- match.arg(kind)
            if(kind == "very short") # e.g. for show() which has more parts
              return(paste0(prefix, "AS-CQS copula"))
            
            name <- "asymmetric"
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

dASC2 <- function (u, copula, log=FALSE) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  
  u1 <- u[, 1]
  u2 <- u[, 2]
  
  if(log)
    return(log(pmax(a * u2 * (((12 - 9 * u1) * u1 - 3) * u2 + u1 * (6 * u1 - 8) + 2) + b * (u2 * ((u1 * (9 * u1 - 12) + 3) * u2 + (12 - 6 * u1) * u1 - 4) - 2 * u1 + 1) + 1,0)))
  else
    return(pmax(a * u2 * (((12 - 9 * u1) * u1 - 3) * u2 + u1 * (6 * u1 - 8) + 2) + b * (u2 * ((u1 * (9 * u1 - 12) + 3) * u2 + (12 - 6 * u1) * u1 - 4) - 2 * u1 + 1) + 1,0))
}

setMethod("dCopula", signature("numeric","asCopula"), 
          function(u, copula, log) {
            dASC2(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","asCopula"), dASC2)

## jcdf ##
pASC2 <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  u1 <- u[, 1]
  u2 <- u[, 2]
  return( u1 * u2 + u1 * u2 * (1 - u1) * (1 - u2) * ((a - b) * u2 * (1 - u1) + b) )
}

setMethod("pCopula", signature("numeric", "asCopula"),
          function(u, copula, ...) {
            pASC2(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix", "asCopula"), pASC2)

## partial derivatives ##
## ddu

dduASC2 <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  u1 <- u[, 1]
  u2 <- u[, 2]

  return(u2*(1 + b*(-1 + 2*u1)*(-1 + u2) - (a - b)*(1 - 4*u1 + 3*u1^2)*(-1 + u2)*u2))
}

setMethod("dduCopula", signature("numeric", "asCopula"),
          function(u, copula, ...) {
            dduASC2(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix", "asCopula"), dduASC2)

## ddv
ddvASC2 <- function (u, copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  u1 <- u[, 1]
  u2 <- u[, 2]

  return( u1 + b*(u1-1)*u1*(2*u2-1) - (a - b)*(-1 + u1)^2*u1*u2*(-2 + 3*u2))
}

setMethod("ddvCopula", signature("numeric", "asCopula"),
          function(u, copula, ...) {
            ddvASC2(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix", "asCopula"),ddvASC2)

## random number generator
# incorporating the inverse of the partial derivative that is solved numerically using optimize

## inverse partial derivative 

invdduASC2 <- function (u, copula, y) {
  stopifnot(length(u) == length(y)) 
  
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  if (a != b) { # the cubic case
    # solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    usq <- u^2
    c3 <- (a-b)*(-3*usq+4*u-1)
    c2 <- (a-b)*(1-4*u+3*usq)+b*(- 1 + 2*u)
    c1 <- 1+b*(1-2*u)
    c0 <- -y

    v <- solveCubicEq(c3,c2,c1,c0) # from cqsCopula.R
    
    filter <- function(vec) {
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
  v[u != 0.5] <- (sqrt(4*y[u != 0.5]*(2*b*uR-b)+(-2*b*uR+b+1)^2)+2*b*uR-b-1)/(2*b*(2*uR-1))
    
  return(v)
}

setMethod("invdduCopula", signature("numeric","asCopula","numeric"),invdduASC2)

## inverse partial derivative ddv
invddvASC2 <- function (v, copula, y) {
  stopifnot(length(v) == length(y)) 

  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  if (a != b) { # the cubic case
    # solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    vsq <- v^2
    c3 <- (a-b)*(2*v-3*vsq)
    c2 <- (a-b)*(-4*v+6*vsq)+b*(-1+2*v)
    c1 <- 1+(a-b)*(2*v - 3*vsq)+b*(1-2*v)
    c0 <- -y

    u <- solveCubicEq(c3,c2,c1,c0) # from cqsCopula.R

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
  u[v != 0.5] <- (sqrt(4*y[v != 0.5]*(2*b*vR-b)+(-2*b*vR+b+1)^2)+2*b*vR-b-1)/(2*b*(2*vR-1))
  
  return(u)
}

setMethod("invddvCopula", signature("numeric","asCopula","numeric"),invddvASC2)

## random number generator
rASC2 <- function (n, copula) {
  u <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
    
  res <- cbind(u, invdduASC2(u, copula, y))
  colnames(res) <- c("u","v")
  
  return(res)
}

setMethod("rCopula", signature("numeric", "asCopula"), rASC2)

## fitment

fitCopulaASC2 <- function (copula, data, method = "ml", start=c(0,0),
                           lower=c(-3,-1), upper=c(1,1), 
                           optim.method="L-BFGS-B", optim.control=list(),
                           estimate.variance = FALSE, call) {
  if(missing(call)) 
    call <- match.call()

  fit <- switch(method, 
                ml=fitASC2.ml(copula, data, start, lower, upper, optim.control, optim.method, call),
                itau=fitASC2.itau(copula, data, estimate.variance, call),
                irho=fitASC2.irho(copula, data, estimate.variance, call),
                stop("Implemented methods for copulas in the spCopula package are: ml, itau, and irho."))
  return(fit)
}

setMethod("fitCopula", signature("asCopula"), fitCopulaASC2)

## Fits the type 2 asymmetric copula acoording to a measure of association.
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

fitASC2.itau <- function(copula, data, estimate.variance, tau=NULL, call) {
  if(missing(call)) 
    call <- match.call()
  if(is.null(tau))
    tau <- TauMatrix(data)[1,2]
  
  esti <- fitASC2.moa(tau, data, method="itau")
  copula <- asCopula(esti)
  
  new("fitCopula", estimate = esti, var.est = matrix(NA), 
      loglik = sum(log(dCopula(data, copula))), nsample = nrow(data),
      method = "Inversion of Kendall's tau and MLE", call = call,
      fitting.stats = list(convergence=as.integer(NA)), copula = copula)
}

fitASC2.irho <- function(copula, data, estimate.variance, rho=NULL){
  if(missing(call)) 
    call <- match.call()
  if(is.null(rho))
    rho <- cor(data,method="spearman")[1,2]
  
  esti <- fitASC2.moa(rho, data, method="irho")
  copula <- asCopula(esti)
  
  new("fitCopula", estimate = esti, var.est = matrix(NA), 
      loglik = sum(log(dCopula(data, copula))), nsample = nrow(data),
      method = "Inversion of Spearman's rho and MLE", call=call,
      fitting.stats = list(convergence=as.integer(NA)), copula = copula)
}

fitASC2.moa <- function(moa, data, method="itau", tol=.Machine$double.eps^.5) {
  smpl <- as.matrix(data)
  iFun <- switch(method, 
                 itau=function(p) iTauASC2(p,moa),
                 irho=function(p) iRhoASC2(p,moa))

  sec <- function (parameters) {
    res <- NULL
    for(param in parameters) {
      res <- rbind(res, -sum(log( dASC2(asCopula(c(iFun(param),param)),u=smpl) )))
    }
    return(res)
  }

  b <- optimize(sec,c(-1,1), tol=tol)$minimum

  return(c(iFun(b),b))
}

# maximum log-likelihood estimation of a and b using optim

fitASC2.ml <- function(copula, data, start, lower, upper, optim.control, optim.method, call) { 
  if(missing(call)) 
    call <- match.call()
  if(length(start)!=2) 
    stop("Start values need to have same length as parameters:")
  
  optFun <- function(param=c(0,0)) {
    if(any(param > 1) | param[2] < -1 | param[1] < limA(param[2])) return(1)
    return(-sum(log( dASC2(data, asCopula(param)))))
  }
  
  optimized <- optim(par=start, fn=optFun, method = optim.method, 
                     lower=lower, upper=upper, control = optim.control)
  
  return(new("fitCopula", 
             estimate = optimized$par, 
             var.est = matrix(NA), 
             loglik = -optimized$value,
             nsample = nrow(data),
             method = "Numerical MLE over the full range.", 
             call = call,
             fitting.stats = optimized,
             copula = asCopula(optimized$par)))
}

####

iTauASC2 <- function(b,tau=0) {
return(min(max(limA(b),(450*tau-75*b+b^2)/(25-b)),1))
}

####

tauASC2 <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return((75*b-b^2+a*(25-b))/450)
}

setMethod("tau",signature("asCopula"),tauASC2)

####
# find parameter "a" for parameter "b" under a given measure of association "rho" 
# it may return a value exceeding the limit of "a" which may result in an invalid copula.

iRhoASC2 <- function(b,rho=0) {
  return(min(max(limA(b),12*rho-3*b),1))
}

####

rhoASC2 <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return((a+3*b)/12)
}

setMethod("rho", signature("asCopula"), rhoASC2)

setMethod("lambda", signature("asCopula"), 
          function(copula, ...) c(lower = 0, upper = 0))