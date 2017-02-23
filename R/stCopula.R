######################################
## Spatio-Temporal Bivariate Copula ##
######################################

## constructor ##
#################

stCopula <- function(components, tlags, distances=NA, stDepFun, unit="m", tres="day") {
  if(all(sapply(components, function(x) class(x)=="spCopula"))) {
    if(length(unique(sapply(components, function(x) x@unit))) >1 )
      stop("All spatial copulas need to have the same distance unit.")
    stopifnot(length(tlags) == length(components))
    spCopList <- components
  } else {
    spCopList <- list()
    
    if(!missing(stDepFun)) {
      getSpCop <- function(comp,dist,time) spCopula(comp, dist,
                                                    spDepFun=function(h) stDepFun(h, time, 1:length(tlags)), unit)
      for(i in 1:length(tlags)){
        spCopList <- append(spCopList, getSpCop(components[[i]], distances[[i]], i))
      }
    } else {
      for(i in 1:length(tlags)){
        spCopList <- append(spCopList, spCopula(components[[i]], distances[[i]], unit=unit))
      }
    }
  }
  
  param       <- unlist(lapply(spCopList, function(x) x@parameters))
  param.names <- unlist(lapply(spCopList, function(x) x@param.names))
  param.low   <- unlist(lapply(spCopList, function(x) x@param.lowbnd))
  param.up    <- unlist(lapply(spCopList, function(x) x@param.upbnd))
  
  new("stCopula", dimension=as.integer(2), parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      fullname="Spatio-Temporal Copula: distance and time dependent convex combination of bivariate copulas",
      spCopList=spCopList, tlags=tlags, tres=tres)
}

## show method ##
#################

showStCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@spCopList)) {
    cmpCop <- object@spCopList[[i]]
    cat("  ", describeCop(cmpCop, "very short"), "at", object@tlags[i], 
        paste("[",object@tres,"]",sep=""), "\n")
    show(cmpCop)
  }
}

setMethod("show", signature("stCopula"), showStCopula)

## spatial copula cdf ##
########################

pStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- pSpCopula(u, copula@spCopList[[match(tDist, copula@tlags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- pSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod(pCopula, signature("numeric","stCopula"), 
          function(u, copula, log, ...) pStCopula(matrix(u,ncol=2), copula, ...))
setMethod(pCopula, signature("matrix","stCopula"), pStCopula)

## spatial Copula density ##
############################

dStCopula <- function (u, copula, log, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- dSpCopula(u, copula@spCopList[[match(tDist, copula@tlags)]], log, h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- dSpCopula(u[tmpInd,,drop=F], tmpCop, log, h[tmpInd,1])
    }
  }
  res
}

setMethod(dCopula, signature("numeric","stCopula"), 
          function(u, copula, log, ...) dStCopula(matrix(u,ncol=2), copula, log=log, ...))
setMethod(dCopula, signature("matrix","stCopula"), dStCopula)


## partial derivatives ##

## dduSpCopula ##
#################

dduStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not match the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- dduSpCopula(u, copula@spCopList[[match(tDist, copula@tlags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- dduSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod("dduCopula", signature("numeric","stCopula"), 
          function(u, copula, ...) dduStCopula(matrix(u,ncol=2), copula, ...))
setMethod("dduCopula", signature("matrix","stCopula"), dduStCopula)

invdduStCopula <- function(u, copula, y, h, tol=.Machine$double.eps^0.5) {
  message("invdduCopula is numerically evalauted.")
  
  if(!is.matrix(h)) 
    h <- matrix(h,ncol=2)
  
  nElem <- length(u)
  stopifnot(nElem == length(y))
  stopifnot(nrow(h) == 1 | nrow(h)==nElem)
  
  optFun <- function(u, v, y, h) abs(dduStCopula(cbind(rep(u, length(v)), v), copula, h)-y)
  
  optMe <- function(aU, aY, aH) optimise(function(v) optFun(u=aU, v, y=aY, h=aH), c(0,1))$minimum
  
  if(nrow(h) == 1 & nElem > 1)
    h <- matrix(rep(h,nElem),ncol=2, byrow=T)
  
  rV <- numeric(nElem)
  for (i in 1:nElem) {
    rV[i] <- optMe(u[i], y[i], h[i,,drop=FALSE])
  }
  
  return(rV)
}

setMethod("invdduCopula", signature("numeric", "stCopula"), invdduStCopula)

## ddvSpCopula ##
#################

ddvStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not match the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- ddvSpCopula(u, copula@spCopList[[match(tDist,copula@tlags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- ddvSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod("ddvCopula", signature("numeric","stCopula"), 
          function(u, copula, ...) ddvStCopula(matrix(u,ncol=2), copula, ...))
setMethod("ddvCopula", signature("matrix","stCopula"), ddvStCopula)

invddvStCopula <- function(v, copula, y, h, tol=.Machine$double.eps^0.5) {
  message("invdduCopula is numerically evalauted.")
  
  if(!is.matrix(h)) 
    h <- matrix(h,ncol=2)
  
  nElem <- length(v)
  stopifnot(nElem == length(y))
  stopifnot(nrow(h) == 1 | nrow(h)==nElem)
  
  optFun <- function(u, v, y, h) abs(ddvStCopula(cbind(u, rep(v, length(u))), copula, h)-y)
  
  optMe <- function(aV, aY, aH) optimise(function(u) optFun(u, v=aV, y=aY, h=aH), c(0,1))$minimum
  
  if(nrow(h) == 1 & nElem > 1)
    h <- matrix(rep(h,nElem),ncol=2, byrow=T)
  
  rU <- numeric(nElem)
  for (i in 1:nElem) {
    rU[i] <- optMe(v[i], y[i], h[i,,drop=FALSE])
  }
  
  return(rU)
}

setMethod("invddvCopula", signature("numeric", "stCopula"), invddvStCopula)

# log-likelihood by copula for all spatio-temporal lags


loglikByCopulasStLags <- function(stBins, data, families = c(normalCopula(),
                                                             tCopula(),
                                                             claytonCopula(),
                                                             frankCopula(),
                                                             gumbelCopula()),
                                  calcCor, lagSub=1:length(stBins$meanDists)) {
  nTimeLags <- dim(stBins$lagCor)[1]
  if(is.null(nTimeLags))
    nTimeLags <- 1
  var <- attr(stBins, "variable")
  
  retrieveData <- function(spIndex, tempIndices) {
    binnedData <- NULL
    for (i in 1:(ncol(tempIndices)/2)) {
      binnedData <- cbind(binnedData, 
                          as.matrix((cbind(data[spIndex[,1], tempIndices[,2*i-1], var, drop=FALSE]@data, 
                                           data[spIndex[,2], tempIndices[,2*i], var, drop=FALSE]@data))))
    }
    return(binnedData)
  }
  
  lagData <- lapply(stBins$lags[[1]][lagSub], retrieveData, tempIndices=stBins$lags[[2]])
  
  tmpBins <- list(meanDists=stBins$meanDists[lagSub])
  attr(tmpBins, "variable") <- var
  
  loglikTau <- list()
  for(j in 1:nTimeLags) {
    tmpLagData <- lapply(lagData, function(x) x[,c(2*j-1,2*j)])
    tmpLagData <- lapply(tmpLagData, function(pairs) {
      bool <- !is.na(pairs[,1]) & !is.na(pairs[,2])
      pairs[bool,]
    })
    
    if(missing(calcCor))
      res <- loglikByCopulasLags.static(tmpLagData, families)
    else
      res <- loglikByCopulasLags.dyn(tmpBins, tmpLagData, families, 
                                     function(h) calcCor(h, j, 1:nTimeLags))
    loglikTau[[paste("loglik",j,sep="")]] <- res
  }
  
  return(loglikTau)
}

# dropping a spatio-temporal tree
dropStTree <- function (stNeigh, dataLocs, stCop) {
  stopifnot(class(stNeigh) == "stNeighbourhood")
  
  u0 <- as.matrix(stNeigh@data)
  h0 <- stNeigh@distances
  u1 <- matrix(NA, nrow(u0), ncol(u0)-1-length(stNeigh@coVar))
  h1 <- array(dim = c(nrow(u0), ncol(h0)-1, 2))
  
  pb <- txtProgressBar(0,dim(h0)[2],style=3)
  for (i in 1:dim(h0)[2]) {
    u1[,i] <- dduCopula(u0[, c(1, i + 1)], stCop, h = h0[, i, ])
    if (i < ncol(h0)) {
      h1[,i,1] <- apply(stNeigh@index[, c(1, i + 1), 1], 1, 
                        function(x) spDists(dataLocs@sp[x, ])[1, 2])
      h1[,i,2] <- apply(stNeigh@index[, c(1, i + 1), 2], 1, 
                        function(x) diff(x))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #   # add covariate to the conditioned neighbourhood?
  #   if (length(stNeigh@coVar) > 0)
  #     u1[,ncol(u0)-(1:length(stNeigh@coVar))] <- u0[,ncol(u0) + 1 - (1:length(stNeigh@coVar))]
  
  varSplit <- strsplit(stNeigh@var, "|", fixed = TRUE)[[1]]
  cond <- suppressWarnings(as.numeric(varSplit[length(varSplit)]))
  
  if (is.na(cond)) {
    #     coVar <- paste(stNeigh@coVar, "|0", sep = "")
    cond <- paste(stNeigh@var, "|0", sep = "")
  }
  else {
    #     coVar <- stNeigh@coVar
    cond <- paste(stNeigh@var, cond + 1, sep = "")
  }
  
  return(stNeighbourhood(data = u1, distances = h1, index = stNeigh@index[, -1, ],
                         var = cond, prediction = stNeigh@prediction))
}