##############################
##                          ##
## a general mixture copula ##
##                          ##
##############################

# class
setClass("mixtureCopula", contains = "copula", slots = list(memberCops= "list"))

# constructor
mixtureCopula <- function (param = c(0.2, 0.2, 0.5), memberCops = c(normalCopula(0), claytonCopula(1))) {
  stopifnot(length(memberCops) == 2)
  stopifnot(memberCops[[1]]@dimension == memberCops[[2]]@dimension)
  
  cop1.nPar <- length(memberCops[[1]]@parameters)
  cop2.nPar <- length(memberCops[[2]]@parameters)
  
  if (missing(param))
    param <- 0.5
  if (length(param) == 1)
    param <- c(memberCops[[1]]@parameters, memberCops[[2]]@parameters, param)
  else {
    stopifnot(length(param) == cop1.nPar + cop2.nPar + 1)
  
    memberCops[[1]]@parameters <- param[1:cop1.nPar]
    memberCops[[2]]@parameters <- param[(1:cop2.nPar)+cop1.nPar]
  }
  
  new("mixtureCopula", dimension = memberCops[[1]]@dimension, parameters = param, memberCops = memberCops,
      param.names = c(memberCops[[1]]@param.names, memberCops[[2]]@param.names, "mixLambda"),
      param.lowbnd = c(memberCops[[1]]@param.lowbnd, memberCops[[2]]@param.lowbnd, 0),
      param.upbnd = c(memberCops[[1]]@param.upbnd, memberCops[[2]]@param.upbnd, 1), 
      fullname = paste("mixture of a", describeCop(memberCops[[1]], "very short"),
                       "and a", describeCop(memberCops[[2]], "very short")))
}

## density ##
setMethod("dCopula", signature(copula = "mixtureCopula"), 
          function(u, copula, log, ...) {
            mixLambda <- tail(copula@parameters, 1)
            res <- (1-mixLambda) * dCopula(u, copula@memberCops[[1]], ...) + mixLambda * dCopula(u, copula@memberCops[[2]], ...)
            if (log)
              return(log(res))
            else 
              return(res)
          })

## jcdf ##
setMethod("pCopula", signature( copula = "mixtureCopula"),
          function(u, copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * pCopula(u, copula@memberCops[[1]]) + mixLambda * pCopula(u, copula@memberCops[[2]])
          })

## partial derivatives ##
## ddu

setMethod("dduCopula", signature(copula = "mixtureCopula"),
          function(u, copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * dduCopula(u, copula@memberCops[[1]]) + mixLambda * dduCopula(u, copula@memberCops[[2]])
          })

# ddv
setMethod("ddvCopula", signature(copula = "mixtureCopula"),
          function(u, copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * ddvCopula(u, copula@memberCops[[1]]) + mixLambda * ddvCopula(u, copula@memberCops[[2]])
          })

## inverse partial derivative 
# invddu
invdduMixCop <- function (u, copula, y) {
  stopifnot(length(u) == length(y)) 
  
  opti <- function(ind) {
    optFun <- function(v) {
      (dduCopula(cbind(u[ind], v), copula) - y[ind])^2
    }
    optimise(optFun, c(0,1))$minimum
  }
  
  sapply(1:length(y), opti)  
}

setMethod("invdduCopula", 
          signature("numeric", "mixtureCopula", "numeric"), 
          invdduMixCop)

# invddv
invddvMixCop <- function (v, copula, y) {
  stopifnot(length(v) == length(y)) 
  
  opti <- function(ind) {
    optFun <- function(u) {
      (dduCopula(cbind(u, v[ind]), copula) - y[ind])^2
    }
    optimise(optFun, c(0,1))$minimum
  }
  
  sapply(1:length(y), opti)  
}

setMethod("invddvCopula", 
          signature("numeric", "mixtureCopula", "numeric"),
          invddvMixCop)

## random number generator

rMixCop <- function(n, copula, ...) {
  u <- runif(n)
  y <- runif(n)
  
  cbind(u, invdduCopula(u, copula, y))
}

setMethod("rCopula", signature(copula = "mixtureCopula"), rMixCop)

## fitment
fitMixCop <- function(copula, data, start, method="mpl",
                      lower = NULL, upper = NULL, 
                      optim.method = "L-BFGS-B", 
                      optim.control = list(maxit = 1000)){
  if (missing(start))
    start <- copula@parameters
  stopifnot(method %in% c("ml", "mpl"))
  
  if(any(is.na(start)))
    stop("Copula parameters or 'start' contains an 'NA' value.")
  
  if(is.null(lower))
    lower <- copula@param.lowbnd
  if(is.null(upper))
    upper <- copula@param.upbnd
  
  optFun <- function(parSet) {
    cop <- mixtureCopula(parSet, copula@memberCops)
    -sum(log(dCopula(data, cop)))
  }
  
  optOut <- optim(start, optFun, method = optim.method, 
                  lower = lower, upper = upper, 
                  control = optim.control)
  
  new("fitCopula",
      copula = mixtureCopula(optOut$par, copula@memberCops),
      estimate = optOut$par,
      var.est = matrix(NA),
      loglik = -optOut$value,
      nsample = nrow(data),
      method = method,
      fitting.stats = append(optOut[c("convergence","counts","message")],
                             optim.control))
}

setMethod(fitCopula, 
          signature = c(copula = "mixtureCopula"), 
          fitMixCop)

## 

setMethod("tau", signature = c(copula = "mixtureCopula"),
          function(copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * tau(copula@memberCops[[1]], ...) + mixLambda * tau(copula@memberCops[[2]], ...)
          })

setMethod("rho", signature = c(copula = "mixtureCopula"),
          function(copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * rho(copula@memberCops[[1]], ...) + mixLambda * rho(copula@memberCops[[2]], ...)
          })

setMethod("lambda", signature = c(copula = "mixtureCopula"),
          function(copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * lambda(copula@memberCops[[1]], ...) + mixLambda * lambda(copula@memberCops[[2]], ...)
          })