# constructor
# dimension = "integer"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# components="list"         list of copulas (will be automatically supplemented 
#			      by the independent copula)
# distances="numeric"       the linking distances + the range (will be assigned
#			      to the independent copula)
# unit="character"          measurement unit of distance
# depFun="function"         a optional spatial dependence function providing 
#                             Kendalls tau or Spearman's rho to calib* or exact 
#                             parameters


spCopula <- function(components, distances, spDepFun, unit="m") {
  if (missing(spDepFun)) { 
    calibMoa <- function(copula, h) return(NULL)
  } else {
    if (is.na(match(spDepFun(NULL), c("kendall","spearman","id")))) 
      stop("spDepFun(NULL) must return 'spearman', 'kendall' or 'id'.")
    cat("The parameters of the components will be recalculated according to the provided spDepFun where possible. \nIn case no 1-1 relation is known, the copula as in components is used. \n")
    calibMoa <- switch(spDepFun(NULL), 
                       kendall=function(copula, h) iTau(copula, spDepFun(h)),
                       spearman=function(copula, h) iRho(copula, spDepFun(h)),
                       id=function(copula, h) return(h))
    
    for (i in 1:length(components)) {
      if(class(components[[i]]) != "indepCopula")
        param <- try(calibMoa(components[[i]], distances[i]),T)
      if (class(param) == "numeric")
        components[[i]]@parameters[1:length(param)] <- param
    }
  }

#   components <- append(components,indepCopula())
  
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
  
  new("spCopula", dimension=as.integer(2), parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      fullname="Spatial Copula: distance dependent convex combination of bivariate copulas",
      components=components, distances=distances, calibMoa=calibMoa, unit=unit)
}

## show method
showCopula <- function(object) {
  cat(object@fullname, "\n")
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

## spatial copula jcdf ##


# for spatial copulas with a spatial dependece function
spDepFunCop <- function(fun, copula, pairs, h, do.logs=F, ...) {
  dists <- copula@distances
  n.dists <- length(dists)
  calibPar <- copula@calibMoa
  
  # data is sorted to be ascending in distance
  res <- numeric(0)
  sel <- which(h < dists[1])
  if(sum(sel)>0) {
    tmpH <- h[sel]
    tmpCop <- copula@components[[1]]
    tmpPairs <- pairs[sel,,drop=FALSE]
    if(class(tmpCop) == "indepCopula")
      res <- fun(tmpPairs, tmpCop, ...)
    else {
      for (j in 1:length(tmpH)) {
        tmpParam <- calibPar(tmpCop, tmpH[j])
        tmpCop@parameters[1:length(tmpParam)] <- tmpParam
        res <- c(res, fun(tmpPairs[j,], tmpCop, ...))
      }
    }
  }
  
  if (n.dists >= 2) {
    for ( i in 2:n.dists ) {
      low  <- dists[i-1]
      high <- dists[i]
      sel <- which(h >= low & h < high)
      if (sum(sel)>0) {
        tmpH <- h[sel]
        tmpPairs <- pairs[sel,,drop=FALSE]
        lowerCop <- copula@components[[i-1]]
        upperCop <- copula@components[[i]]
        if (class(lowerCop) != class(upperCop)) {
          lowerVals <- numeric(0)
          upperVals <- numeric(0)
          for (j in 1:length(tmpH)) {
            if (class(lowerCop) != "indepCopula") {
              lowerParam <- calibPar(lowerCop,  tmpH[j])
              lowerCop@parameters[1:length(lowerParam)] <- lowerParam
            }
              
            if (class(upperCop) != "indepCopula") {
              upperParam <- calibPar(upperCop, tmpH[j])
              upperCop@parameters[1:length(upperParam)] <- upperParam
            }
            lowerVals <- c(lowerVals, fun(tmpPairs[j,], lowerCop))
            upperVals <- c(upperVals, fun(tmpPairs[j,], upperCop))
          }
          if(do.logs)
            res <- c(res,
                     log((high-tmpH)/(high-low)*lowerVals+(tmpH-low)/(high-low)*upperVals))
          else 
            res <- c(res,(high-tmpH)/(high-low)*lowerVals+(tmpH-low)/(high-low)*upperVals)
        } else {
          if (class(lowerCop) == "indepCopula")
            newVals <- fun(tmpPairs, lowerCop, ...)
          else {
            newVals <- numeric(0)
            for (j in 1:length(tmpH)) {
              lowerParam <- calibPar(lowerCop, tmpH[j])
              lowerCop@parameters[1:length(lowerParam)] <- lowerParam
              newVals <- c(newVals, fun(tmpPairs[j,], lowerCop, ...))
            }
          }
          res <- c(res, newVals)
        }
      }
    }
  }
  
  sel <- which(h >= dists[n.dists])
  if(sum(sel)>0) {
    res <- c(res, fun(pairs[which(h >= dists[n.dists]),],
                      copula@components[[n.dists]], ...))
  }

  return(res)
}

# for spatial copulas with a spatial dependece function and a single distance but many pairs
spDepFunCopSnglDist <- function(fun, copula, pairs, h, do.logs=F, ...) {
  dists <- copula@distances
  n.dists <- length(dists)
  calibPar <- copula@calibMoa

  if(h < dists[1]) {
    tmpCop <- copula@components[[1]]
    if(class(tmpCop) != "indepCopula") {
      tmpParam <- calibPar(tmpCop, h)
      tmpCop@parameters[1:length(tmpParam)] <- tmpParam
    }
      
    return(fun(pairs, tmpCop, ...))
  }
  
  if(h >= dists[n.dists]) {
    return(fun(pairs, copula@components[[n.dists]], ...))
  }
  
  for (i in 2:n.dists) {
    low  <- dists[i-1]
    high <- dists[i]
    if(low <= h & h < high) {
      lowerCop <- copula@components[[i-1]]
      upperCop <- copula@components[[i]]
      if (class(lowerCop) != class(upperCop)) {
        if (class(lowerCop) != "indepCopula") {
          lowerParam <- calibPar(lowerCop,  h)
          lowerCop@parameters[1:length(lowerParam)] <- lowerParam
        }
        if (class(upperCop) != "indepCopula") {
          upperParam <- calibPar(upperCop, h)
          upperCop@parameters[1:length(upperParam)] <- upperParam
        }

        lowerVals <- fun(pairs, lowerCop)
        upperVals <- fun(pairs, upperCop)

        res <- (high-h)/(high-low)*lowerVals + (h-low)/(high-low)*upperVals
        if(do.logs)
          return(log(res))
        return(res)
      } else {
        if(class(lowerCop) != "indepCopula") {
          lowerParam <- calibPar(lowerCop, h)
          lowerCop@parameters[1:length(lowerParam)] <- lowerParam
        }
        return(fun(pairs, lowerCop, ...))
      }
    }
  }
}


# for static convex combinations of copulas
spConCop <- function(fun, copula, pairs, h, do.logs=F, ...) {
  dists <- copula@distances
  n.dists <- length(dists)
  
  res <- numeric(nrow(pairs))
  sel <- which(h < dists[1])
  if(sum(sel)>0) {
    res[sel] <- fun(pairs[sel,,drop=FALSE],copula@components[[1]], ...)
  }
  
  if (n.dists >= 2) {
    for ( i in 2:n.dists ) {
      low  <- dists[i-1]
      high <- dists[i]
      sel <- which(h >= low & h < high)
      if (sum(sel)>0) {
        tmpH <- h[sel]
        tmpPairs <- pairs[sel,,drop=FALSE]

        lowerVals <- fun(tmpPairs[,], copula@components[[i-1]])
        upperVals <- fun(tmpPairs[,], copula@components[[i]])

        if(do.logs)
          res[sel] <- log((high-tmpH)/(high-low)*lowerVals+(tmpH-low)/(high-low)*upperVals)
        else
          res[sel] <- (high-tmpH)/(high-low)*lowerVals+(tmpH-low)/(high-low)*upperVals
      }
    }
  }
  
  sel <- which(h >= dists[n.dists])
  if(sum(sel)>0) {
    res[sel] <- fun(pairs[which(h >= dists[n.dists]),],
                    copula@components[[n.dists]], ...)
  }

  return(res)
}


## spatial copula CDF
######################

pSpCopula <- function (u, copula, h, block=1) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance \"h\".")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  
  if(is.null(copula@calibMoa(normalCopula(0),0)))
    return(spConCop(pCopula, copula, u, rep(h,length.out=nrow(u))))
  
  if(length(h)>1) {
    ordering <- order(h)
        
    # ascending sorted pairs allow for easy evaluation
    u <- u[ordering,,drop=FALSE] 
    h <- h[ordering]
        
    res <- spDepFunCop(pCopula, copula, u, h)
        
    # reordering the values
    return(res[order(ordering)])
  } else
      return(spDepFunCopSnglDist(pCopula, copula, u, h))
}

setMethod(pCopula, signature("numeric","spCopula"), 
          function(u, copula, ...) pSpCopula(matrix(u,ncol=2),copula, ...))
setMethod(pCopula, signature("matrix","spCopula"), pSpCopula)

## spatial Copula density 
##########################

dSpCopula <- function (u, copula, log, h) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance \"h\".")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  
  if(is.null(copula@calibMoa(normalCopula(0),0)))
    return(spConCop(dCopula, copula, u, rep(h, length.out=nrow(u)), do.logs=log,
                    log=log))

  if(length(h)>1) {
    ordering <- order(h)
    
    # ascending sorted pairs allow for easy evaluation
    u <- u[ordering,,drop=FALSE] 
    h <- h[ordering]
        
    res <- spDepFunCop(dCopula, copula, u, h, do.logs=log, log=log)
        
    # reordering the values
    return(res[order(ordering)])
  } else
      return(spDepFunCopSnglDist(dCopula, copula, u, h, do.logs=log, log=log))
}

setMethod(dCopula, signature("numeric","spCopula"), 
          function(u, copula, log, ...) dSpCopula(matrix(u,ncol=2), copula, log=log, ...))
setMethod(dCopula, signature("matrix","spCopula"), dSpCopula)

## partial derivatives ##
## dduSpCopula
###############

dduSpCopula <- function (u, copula, h) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance h.")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")

  if(is.null(copula@calibMoa(normalCopula(0),0)))
    return(spConCop(dduCopula, copula, u, rep(h, length.out=nrow(u))))
  
  if(length(h)>1) {
    ordering <- order(h)
        
    # ascending sorted pairs allow for easy evaluation
    u <- u[ordering, , drop=FALSE] 
    h <- h[ordering]
        
    res <- spDepFunCop(dduCopula, copula, u, h)      
        
    # reordering the values
    return(res[order(ordering)])
  } else
    return(spDepFunCopSnglDist(dduCopula, copula, u, h))
}

setMethod("dduCopula", signature("matrix","spCopula"), dduSpCopula)
setMethod("dduCopula", signature("numeric","spCopula"), 
          function(u, copula, ...) dduSpCopula(matrix(u,ncol=copula@dimension),copula, ...) )

invdduSpCopula <- function(u, copula, y, h, tol=.Machine$double.eps^0.5) {
  message("invdduCopula is numerically evalauted.")
  
  nElem <- length(u)
  stopifnot(nElem == length(y))
  stopifnot(length(h) == 1 | length(h)==nElem)
  
  optFun <- function(u, v, y, h) abs(dduSpCopula(cbind(rep(u, length(v)), v), copula, h)-y)
  
  optMe <- function(aU, aY, aH) optimise(function(v) optFun(u=aU, v, y=aY, h=aH), c(0,1))$minimum
  
  if(length(h) == 1 & nElem > 1)
    h <- rep(h, nElem)
  
  rV <- numeric(nElem)
  for (i in 1:nElem) {
    rV[i] <- optMe(u[i], y[i], h[i])
  }
  
  return(rV)
}

setMethod("invdduCopula", signature("numeric", "spCopula"), invdduSpCopula)

## ddvSpCopula
###############

ddvSpCopula <- function (u, copula, h) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance h.")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
    
  if(is.null(copula@calibMoa(normalCopula(0),0)))
    return(spConCop(ddvCopula, copula, u, rep(h, length.out=nrow(u))))
  
  if(length(h)>1) {
    ordering <- order(h)
        
    # ascending sorted pairs allow for easy evaluation
    u <- u[ordering,,drop=FALSE] 
    h <- h[ordering]
        
    res <- spDepFunCop(ddvCopula, copula, u, h)      
        
    # reordering the values
    return(res[order(ordering)])
  } else
      return(spDepFunCopSnglDist(ddvCopula, copula, u, h))
}

setMethod("ddvCopula", signature("matrix","spCopula"), ddvSpCopula)
setMethod("ddvCopula", signature("numeric","spCopula"), 
          function(u, copula, ...) ddvSpCopula(matrix(u,ncol=copula@dimension),copula, ...) )

invddvSpCopula <- function(v, copula, y, h, tol=.Machine$double.eps^0.5) {
  message("invddvCopula is numerically evalauted.")
  
  nElem <- length(v)
  stopifnot(nElem == length(y))
  stopifnot(length(h) == 1 | length(h)==nElem)
  
  optFun <- function(u, v, y, h) abs(ddvSpCopula(cbind(u, rep(v, length(u))), copula, h)-y)
  
  optMe <- function(aV, aY, aH) optimise(function(u) optFun(u, v=aV, y=aY, h=aH), c(0,1))$minimum
  
  if(length(h) == 1 & nElem > 1)
    h <- rep(h, nElem)
  
  rU <- numeric(nElem)
  for (i in 1:nElem) {
    rU[i] <- optMe(v[i], y[i], h[i])
  }
  
  return(rU)
}

setMethod("invddvCopula", signature("numeric", "spCopula"), invddvSpCopula)

## simulation

spCop.rCop <- function(n, copula, h) {
  u <- runif(n)
  v <- invdduCopula(u, copula, y=runif(n), h=h)
  
  return(cbind(u, v))
}

setMethod("rCopula", signature("numeric", "spCopula"), spCop.rCop)

#############
##         ##
## FITTING ##
##         ##
#############

# two models: 
# 1) Kendall's tau driven:
#    fit curve through emp. Kendall's tau values, identify validity ranges for
#    copula families deriving parameters from the fit, fade from one family to 
#    another at borders
# 2) convex-linear combination of copulas: 
#    fit one per lag, fade from one to another

# towards the first model:

# INPUT: the stBinning
# steps
# a) fit a curve
# b) estimate bivariate copulas per lag (limited to those with some 1-1-relation 
#    to Kendall's tau')
# INTERMEDIATE RESULT
# c) select best fits based on ... e.g. log-likelihood, visual inspection
# d) compose bivariate copulas to one spatial copula
# OUTPUT: a spatial copula parametrised by distance through Kendall's tau

# towards a)
# bins   -> typically output from calcBins
# degree -> the degree of the polynominal
# cutoff -> maximal distance that should be considered for fitting
# bounds -> the bounds of the correlation function (typically c(0,1))
# method -> the measure of association, either "kendall" or "spearman"
fitCorFunSng <- function(bins, degree, cutoff, bounds, cor.method, weighted) {
  if (weighted) {
    bins <- as.data.frame(bins[c("np","meanDists","lagCor")])
    if(!is.na(cutoff)) 
      bins <- bins[bins$meanDists <= cutoff,]
    fitCor <- lm(lagCor ~ poly(meanDists, degree), data = bins, weights=bins$np)
  } else {
    bins <- as.data.frame(bins[c("meanDists","lagCor")])
    if(!is.na(cutoff)) 
      bins <- bins[bins$meanDists <= cutoff,]
    fitCor <- lm(lagCor ~ poly(meanDists, degree), data = bins)
  }
  
  print(fitCor)
  cat("Sum of squared residuals:",sum(fitCor$residuals^2),"\n")
  
  if(cor.method=="fasttau") 
    cor.method <- "kendall"
  
  function(x) {
    if (is.null(x)) return(cor.method)
    return(pmin(bounds[2], pmax(bounds[1], 
                                eval(predict(fitCor, data.frame(meanDists=x))))))
  }
}

fitCorFun <- function(bins, degree=3, cutoff=NA, tlags, bounds=c(0,1), 
                      cor.method=NULL, weighted=FALSE){
  if(is.null(cor.method)) {
    if(is.null(attr(bins,"cor.method")))
      stop("Neither the bins arguments has an attribute cor.method nor is the parameter cor.method provided.") 
    else 
      cor.method <- attr(bins,"cor.method")
  } else {
    if(!is.null(attr(bins,"cor.method")) && cor.method != attr(bins,"cor.method"))
      stop("The cor.method attribute of the bins argument and the argument cor.method do not match.")
  }
  
  if(is.null(nrow(bins$lagCor))) # the spatial case
    return(fitCorFunSng(bins, degree, cutoff, bounds, cor.method, weighted))
    
  # the spatio-temporal case
  degree <- rep(degree, length.out = nrow(bins$lagCor))
  calcKTau <- list()
  for (j in 1:nrow(bins$lagCor)) {
    calcKTau[[paste("fun",j,sep="")]] <- fitCorFunSng(data.frame(np=bins$lagNp[j,],
                                                                 meanDists=bins$meanDists, 
                                                                 lagCor=bins$lagCor[j,]),
                                                      degree[j], cutoff, bounds, 
                                                      cor.method, weighted)
  }
  
  tlsort <- sort(tlags,decreasing=TRUE)
  
  corFun <- function(h, time, tlags=tlsort) {
    t <- which(tlags==time)
    calcKTau[[time]](h)
  }
  
  attr(corFun, "tlags") <- sort(tlags, decreasing=TRUE)
  return(corFun)
}


# towards b)
  
## loglikelihoods for a dynamic spatial copula
loglikByCopulasLags.dyn <- function(bins, lagData, families, calcCor) {
  moa <- switch(calcCor(NULL),
                kendall=function(copula, h) iTau(copula, calcCor(h)),
                spearman=function(copula, h) iRho(copula, calcCor(h)),
                id=function(copula, h) calcCor(h))
  
  loglik <- NULL
  copulas <- list()
  for (cop in families) {
    cat(describeCop(cop, "very short"),"\n")
    tmploglik <- NULL
    tmpCop <- list()
    
    pb <- txtProgressBar(0, length(bins$meanDists), style=3)
    for(i in 1:length(bins$meanDists)) {
      if(class(cop)!="indepCopula") {
        if(class(cop) == "asCopula") {
          cop <- switch(calcCor(NULL),
                        kendall=fitASC2.itau(cop, lagData[[i]], 
                                              tau=calcCor(bins$meanDists[i]))@copula,
                        spearman=fitASC2.irho(cop, lagData[[i]],
                                              rho=calcCor(bins$meanDists[i]))@copula,
                        stop(paste(calcCor(NULL), "is not yet supported.")))
          param <- cop@parameters
        } else {
          if(class(cop) == "cqsCopula") {
            cop <- switch(calcCor(NULL),
                          kendall=fitCQSec.itau(cop, lagData[[i]], 
                                                tau=calcCor(bins$meanDists[i]))@copula,
                          spearman=fitCQSec.irho(cop, lagData[[i]],
                                                rho=calcCor(bins$meanDists[i]))@copula,
                          stop(paste(calcCor(NULL), "is not yet supported.")))
            param <- cop@parameters
          } else {
            param <- moa(cop, bins$meanDists[i])
            if(!is.na(param))
              cop@parameters[1:length(param)] <- param
          }
        }
      }
      
      if(any(is.na(param)))
        tmploglik <- c(tmploglik, NA)
      else 
        tmploglik <- c(tmploglik, sum(dCopula(lagData[[i]], cop, log=T)))
      tmpCop <- append(tmpCop, cop)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    loglik <- cbind(loglik, tmploglik)
    copulas[[class(cop)]] <- tmpCop
  }

  colnames(loglik) <- sapply(families, function(x) class(x)[1])

  return(list(loglik=loglik, copulas=copulas))
}

## loglikelihoods for a static spatial copula
loglikByCopulasLags.static <- function(lagData, families) {
  
  fits <-lapply(families, 
                function(cop) {
                  cat(describeCop(cop, "very short"), "\n")
                  lapply(lagData,
                         function(x) {
                           tryCatch(fitCopula(cop, x, estimate.variance = FALSE),
                                    error=function(e) return(NA))
                         })
                })
  
  loglik <- lapply(fits, function(x) sapply(x, function(fit) {
    if(class(fit)=="fitCopula")
      return(fit@loglik)
    else
      return(NA)
  }))
  
  loglik <- matrix(unlist(loglik),ncol=length(loglik),byrow=F)
  colnames(loglik) <- sapply(families, function(x) class(x)[1])
  
  copulas <- lapply(fits, function(x) sapply(x, function(fit) {
    if(class(fit)=="fitCopula")
      return(fit@copula)
    else
      return(NULL)
  }))

  names(copulas) <- colnames(loglik)
  
  return(list(loglik=loglik, copulas=copulas))
}

##

loglikByCopulasLags <- function(bins, data, families=c(normalCopula(), 
                                                       tCopula(),
                                                       claytonCopula(), frankCopula(), 
                                                       gumbelCopula()),
                                calcCor, lagSub=1:length(bins$meanDists)) {
  var <- attr(bins, "variable")
  
  if(missing(data)) {
    lagData <- bins$lagData
  }
  else {
    lagData <- lapply(bins$lags[lagSub], 
                      function(x) {
                        as.matrix((cbind(data[x[, 1], var, drop=FALSE]@data,
                                         data[x[, 2], var, drop=FALSE]@data)))
    })
  }
  
  lagData <- lapply(lagData, 
                    function(pairs) {
                      bool <- !is.na(pairs[,1]) & !is.na(pairs[,2])
                      pairs[bool,]
                    })
  
  if(missing(calcCor))
    return(loglikByCopulasLags.static(lagData, families))
  else
    return(loglikByCopulasLags.dyn(lapply(bins, function(x) x[lagSub]),
                                          lagData, families, calcCor))
}



# towards d)
composeSpCopula <- function(bestFit, families, bins, calcCor, range=max(bins$meanDists)) {
  nFits <- length(bestFit)
  if(nFits > length(bins$meanDists))
    stop("There may not be less bins than best fits.\n")
  rangeIndex <- min(nFits, max(which(bins$meanDists <= range)))
  
  if (missing(calcCor)) {
    return(spCopula(components = as.list(families[bestFit[1:rangeIndex]]),
                    distances = bins$meanDists[1:rangeIndex], 
                    unit = "m"))
  }
  
  else {
    rangeIndex <- min(rangeIndex, which(calcCor(bins$meanDists) <= 0))
    
    return(spCopula(components = as.list(families[bestFit[1:rangeIndex]]),
                    distances = bins$meanDists[1:rangeIndex], 
                    unit = "m", spDepFun = calcCor))
  }
}

# in once

# bins   -> typically output from calcBins
# cutoff -> maximal distance that should be considered for fitting
# families -> a vector of dummy copula objects of each family to be considered
#             DEFAULT: c(normal, t_df=4, clayton, frank, gumbel
# ...
# type   -> the type of curve (by now only polynominals are supported)
# degree -> the degree of the polynominal
# bounds -> the bounds of the correlation function (typically c(0,1))
# method -> the measure of association, either "kendall" or "spearman"
fitSpCopula <- function(bins, data, cutoff=NA, 
                        families=c(normalCopula(), tCopula(),
                                   claytonCopula(), frankCopula(),
                                   gumbelCopula()), ...) {
  calcCor <- fitCorFun(bins, cutoff=cutoff, ...)
  loglik <- loglikByCopulasLags(bins, data, families, calcCor)
  
  bestFit <- apply(apply(loglik$loglik, 1, rank),2, 
                   function(x) which(x==length(families)))
  
  return(composeSpCopula(bestFit, families, bins, calcCor, range=cutoff))
}

## dropping a spatial tree, returning a conditional neighbourhood
dropSpTree <- function(neigh, dataLocs, spCop) {
  u1 <- matrix(NA,nrow(neigh@data),ncol(neigh@data)-1)
  h1 <- matrix(NA,nrow(neigh@distances),ncol(neigh@distances)-1)

  for(i in 1:ncol(neigh@distances)) {
    u1[,i] <- dduCopula(as.matrix(neigh@data[,c(1,1+i)]), spCop, 
                     neigh@distances[,i])
    if (i < ncol(neigh@distances)) {
      h1[,i] <- apply(neigh@index[,c(2,2+i)],1, 
                   function(x) spDists(dataLocs[x,])[1,2])
    }
  }
  
  varSplit <- strsplit(neigh@var,"|",fixed=TRUE)[[1]]
  cond <- suppressWarnings(as.numeric(varSplit[length(varSplit)]))
  coVar <- neigh@coVar
  if(is.na(cond)) {
    var <- paste(neigh@var,"|0",sep="")
    if(length(coVar)>0)
      coVar <- paste(neigh@coVar,"|0",sep="")
    colnames(u1) <- paste(paste("N", rep(1:(ncol(u1)), each=length(var)), sep=""),
                          rep(var,ncol(u1)),sep=".")
  } else {
    var <- paste(neigh@var,cond+1,sep="")
    colnames(u1) <- paste(paste("N", rep(cond:(ncol(u1)+cond-1)+2,
                                         each=length(var)), sep=""),
                          rep(var,ncol(u1)),sep=".")
  }
  return(neighbourhood(data=u1, distances=h1, index=neigh@index[,-1],
                       var=var, coVar=coVar, prediction=neigh@prediction))
}