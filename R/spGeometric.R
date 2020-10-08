# geometric means of copulas

## dists is expected to contain 0 and components is expected to contain a copula 
## for 0-distance, either perfect dependence (nugget free), or some strongly 
## correlated copula; components needs to contain one copula that is used beyond
## the range, typically the product copula

spGeomCopula <- function(components, distances, unit=NULL) {

  param       <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@parameters}))
  param.names <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@param.names}))
  param.low   <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@param.lowbnd}))
  param.up    <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@param.upbnd}))
  
  new("spGeomCopula", dimension=as.integer(2), 
      parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      components=components, distances=distances, unit=unit)
}

## show method
showCopula <- function(object) {
  cat("Spatial Copula based on geometric means.", "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@components)) {
    cmpCop <- object@components[[i]]
    cat("  ", describeCop(cmpCop, "very short"), "at", object@distances[i], 
        paste("[",object@unit,"]",sep=""), "\n")
  }
  if(!is.null(object@calibMoa(normalCopula(0),0))) 
    cat("A spatial dependence function is used. \n")
}

setMethod("show", signature("spCopula"), showCopula)


# CDF #
#######

pCop.spGeomCop <- function(u, copula, h, ...) {
  dists <- copula@distances

  resLow <- numeric(nrow(u))
  resHigh <- numeric(nrow(u))
  p <- numeric(nrow(u))
  
  pairsInt <- findInterval(h, dists, , left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) 
      next;
    if (int >= length(copula@components)) {
      resLow[sel] <- resHigh[sel] <- pCopula(u[sel,,drop=FALSE],
                                             copula@components[[int]], ...)
      
      p[sel] <- 0.5
    } else {
      resLow[sel] <- pCopula(u[sel,,drop=FALSE], 
                             copula@components[[int]], ...)
      resHigh[sel] <- pCopula(u[sel,,drop=FALSE], 
                              copula@components[[int+1]], ...)
    
      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
  resLow^p * resHigh^(1-p)
}

setMethod(pCopula, signature("numeric","spGeomCopula"), 
          function(u, copula, ...) pCop.spGeomCop(matrix(u, ncol=2), copula, ...))
setMethod(pCopula, signature("matrix", "spGeomCopula"), 
          pCop.spGeomCop)

# PDF #
#######

dCop.spGeomCop <- function(u, copula, h, do.logs=F, ...) {
  dists <- copula@distances

  dCopLow <- matrix(NA, nrow(u), 4)
  dCopHigh <- matrix(NA, nrow(u), 4)
  p <- numeric(nrow(u))
  
pairsInt <- findInterval(h, dists, left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) 
      next;
    if (int >= length(copula@components)) {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      dCopHigh[sel,] <- dCopLow[sel,] <- cbind(pCopula(tmpPairs, copLow, ...),
                                               dduCopula(tmpPairs, copLow, ...),
                                               ddvCopula(tmpPairs, copLow, ...),
                                               dCopula(tmpPairs, copLow, ...))
      
      p[sel] <- 0.5
    } else {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      copHigh <- copula@components[[int+1]]
      dCopLow[sel,] <- cbind(pCopula(tmpPairs, copLow, ...),
                             dduCopula(tmpPairs, copLow, ...),
                             ddvCopula(tmpPairs, copLow, ...),
                             dCopula(tmpPairs, copLow, ...))
      dCopHigh[sel,] <- cbind(pCopula(tmpPairs, copHigh, ...),
                             dduCopula(tmpPairs, copHigh, ...),
                             ddvCopula(tmpPairs, copHigh, ...),
                             dCopula(tmpPairs, copHigh, ...))

      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
# f[u,v]^(-1-p)        g[u,v]^(-2+p)         ((-1+p)   (p   (g[u,v]       f^(0,1)[u,v] -f[u,v]      g^(0,1)[u,v])   (g[u,v]       f^(1,0)[u,v]- f[u,v]      g^(1,0)[u,v]) - f[u,v]        g[u,v]^2       f^(1,1)[u,v])+ p   f[u,v]^2      g[u,v]       g^(1,1)[u,v])
  res <- dCopLow[,1]^(-1-p) * dCopHigh[,1]^(-2+p) * ((-1+p) * (p * (dCopHigh[,1]*dCopLow[,3] - dCopLow[,1]*dCopHigh[,3]) * (dCopHigh[,1]*dCopLow[,2] - dCopLow[,1]*dCopHigh[,2]) - dCopLow[,1] * dCopHigh[,1]^2*dCopLow[,4]) + p * dCopLow[,1]^2*dCopHigh[,1]*dCopHigh[,4])
  
  if (do.logs)
    return(log(res))
  else 
    return(res)
}

setMethod(dCopula, signature("numeric","spGeomCopula"), 
          function(u, copula, log, ...) dCop.spGeomCop(matrix(u,ncol=2), 
                                                       copula, log=log, ...))
setMethod(dCopula, signature("matrix","spGeomCopula"), 
          dCop.spGeomCop)

# partial derivative d/du #
###########################

dduCop.spGeomCop <- function(u, copula, h, do.logs=F, ...) {
  dists <- copula@distances
  
  dCopLow <- matrix(NA, nrow(u), 4)
  dCopHigh <- matrix(NA, nrow(u), 4)
  p <- numeric(nrow(u))
  
  pairsInt <- findInterval(h, dists, left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) next;
    
    if (int >= length(copula@components)) {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      dCopHigh[sel, ] <- dCopLow[sel,] <- cbind(pCopula(tmpPairs, copLow, ...),
                                                dduCopula(tmpPairs, copLow, ...))
      
      p[sel] <- 0.5
    } else {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      copHigh <- copula@components[[int+1]] # check for max 
      dCopLow[sel,] <- cbind(pCopula(tmpPairs, copLow, ...),
                             dduCopula(tmpPairs, copLow, ...))
      dCopHigh[sel,] <- cbind(pCopula(tmpPairs, copHigh, ...),
                              dduCopula(tmpPairs, copHigh, ...))
      
      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
# (1-c)   f[u,v]^-c        g[u,v]^c         f^(1,0)[u,v]+ c   f[u,v]^(1-c)        g[u,v]^(-1+c)        g^(1,0)[u,v]
  (1-p) * dCopLow[,1]^-p * dCopHigh[,1]^p * dCopLow[,2] + p * dCopLow[,1]^(1-p) * dCopHigh[,1]^(p-1) * dCopHigh[,2]
}

setMethod("dduCopula", signature("matrix","spGeomCopula"), 
          dduCop.spGeomCop)
setMethod("dduCopula", signature("numeric","spGeomCopula"), 
          function(u, copula, ...) dduCop.spGeomCop(matrix(u, ncol=copula@dimension),
                                                    copula, ...) )

# invddu #
##########

invdduCop.spGeomCop <- function(u, copula, y, h, tol=.Machine$double.eps^0.5) {
  message("invdduCopula is numerically evalauted.")
  
  nElem <- length(u)
  stopifnot(nElem == length(y))
  stopifnot(length(h) == 1 | length(h)==nElem)
  
  optFun <- function(u, v, y, h) abs(dduCop.spGeomCop(cbind(rep(u, length(v)), v), copula, h)-y)
  
  optMe <- function(aU, aY, aH) optimise(function(v) optFun(u=aU, v, y=aY, h=aH), c(0,1))$minimum
  
  if(length(h) == 1 & nElem > 1)
    h <- rep(h, nElem)
  
  rV <- numeric(nElem)
  for (i in 1:nElem) {
    rV[i] <- optMe(u[i], y[i], h[i])
  }
  
  return(rV)
}

setMethod("invdduCopula", signature("numeric", "spGeomCopula"), invdduCop.spGeomCop)

# partial derivative d/dv #
###########################

ddvCop.spGeomCop <- function(u, copula, h) {
  dists <- copula@distances
  
  dCopLow <- matrix(NA, nrow(u), 2)
  dCopHigh <- matrix(NA, nrow(u), 2)
  p <- numeric(nrow(u))
  
  pairsInt <- findInterval(h, dists, left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) next;
    
    if (int >= length(copula@components)) {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      dCopHigh[sel, ] <- dCopLow[sel,] <- cbind(pCopula(tmpPairs, copLow),
                                                ddvCopula(tmpPairs, copLow))
      
      p[sel] <- 0.5
    } else {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      copHigh <- copula@components[[int+1]]
      dCopLow[sel,] <- cbind(pCopula(tmpPairs, copLow),
                             ddvCopula(tmpPairs, copLow))
      dCopHigh[sel,] <- cbind(pCopula(tmpPairs, copHigh),
                              ddvCopula(tmpPairs, copHigh))
      
      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
  # (1-c) f[u,v]^-c          g[u,v]^c         f^(0,1)[u,v]+ c f[u,v]^(1-c)      g[u,v]^(-1+c)        (g^(0,1))[u,v]
  (1-p) * dCopLow[,1]^(-p) * dCopHigh[,1]^p * dCopLow[,2] + p*dCopLow[,1]^(1-p)*dCopHigh[,1]^(p-1) * dCopHigh[,2]
  
}

setMethod("ddvCopula", signature("matrix","spGeomCopula"), 
          ddvCop.spGeomCop)
setMethod("ddvCopula", signature("numeric","spGeomCopula"), 
          function(u, copula, ...) ddvCop.spGeomCop(matrix(u, ncol=copula@dimension),
                                                    copula, ...) )

# invddv #
##########

invddvCop.spGeomCop <- function(v, copula, y, h, tol=.Machine$double.eps^0.5) {
  message("invddvCopula is numerically evalauted.")
  
  nElem <- length(v)
  stopifnot(nElem == length(y))
  stopifnot(length(h) == 1 | length(h) == nElem)
  
  optFun <- function(u, v, y, h) abs(ddvCop.spGeomCop(cbind(u, rep(v, length(u))), copula, h)-y)
  
  optMe <- function(aV, aY, aH) optimise(function(u) optFun(u, v=aV, y=aY, h=aH), c(0,1))$minimum
  
  if(length(h) == 1 & nElem > 1)
    h <- rep(h, nElem)
  
  rU <- numeric(nElem)
  for (i in 1:nElem) {
    rU[i] <- optMe(v[i], y[i], h[i])
  }
  
  return(rU)
}

setMethod("invddvCopula", signature("numeric", "spGeomCopula"), invddvCop.spGeomCop)


## simulation

rCop.spGeomCop <- function(n, copula, h) {
  u <- runif(n)
  v <- invdduCop.spGeomCop(u, copula, y=runif(n), h=h)
  
  return(cbind(u, v))
}

setMethod("rCopula", signature("numeric", "spGeomCopula"), rCop.spGeomCop)
