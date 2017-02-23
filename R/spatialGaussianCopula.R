## spatial Gaussian Copula

# "density" evaluation
spGaussLogLik <- function(corFun, neigh, dataLocs, log=T) {
  neighDim <- ncol(neigh@data)
  
  allDataDists <- spDists(dataLocs)
  
  pb <- txtProgressBar(0, nrow(neigh@data), 0, width = getOption("width") - 10, style = 3)
  
  loglik <- 0
  
  for(i in 1:nrow(neigh@data)) { # i <- 2
    setTxtProgressBar(pb, i)
    tmpDists <- allDataDists[neigh@index[i,], neigh@index[i,]]
    
    tmpCor <- corFun(tmpDists)
    
    tmpGaussCop <- normalCopula(tmpCor[lower.tri(tmpCor)], neighDim, dispstr="un")
    
    loglik <- loglik + dCopula(as.numeric(neigh@data[i,]), tmpGaussCop, log=T)
  }
  close(pb)
  
  if(log)
    return(loglik)
  else
    return(exp(loglik))
}

# interpolation based on a valid corelogram function
spGaussCopPredict <- function(corFun, predNeigh, dataLocs, predLocs, margin, p=0.5, ..., n=1000) {
  stopifnot(is.list(margin))
  stopifnot(is(margin$q, "function"))
  
  neighDim <- ncol(predNeigh@data)
  allDataDists <- spDists(dataLocs)
  
  pb <- txtProgressBar(0, nrow(predNeigh@data), 0, width = getOption("width") - 10, style = 3)
  
  predQuantile <- NULL
  
  for(i in 1:nrow(predNeigh@data)) { # i <- 2
    setTxtProgressBar(pb, i)
    tmpDataDists <- allDataDists[predNeigh@index[i,-1], predNeigh@index[i,-1]]
    
    tmpDists <- rbind(c(0,predNeigh@distances[i,]),
                      cbind(predNeigh@distances[i,], tmpDataDists))
    tmpCor <- corFun(tmpDists)
    
    tmpGaussCop <- normalCopula(tmpCor[lower.tri(tmpCor)], neighDim, dispstr="un")
    rat <- 50:1 %x% c(1e-06, 1e-05, 1e-04, 0.001)
    xVals <- unique(sort(c(rat, 1 - rat, 1:(n - 1)/n)))
    nx <- length(xVals)
    
    condGausCop <- dCopula(cbind(xVals, matrix(rep(as.numeric(predNeigh@data[i,-1]),
                                                   length(xVals)), ncol=neighDim-1, byrow=T)),
                           tmpGaussCop)
    
    condGausCop <- c(max(0, 2 * condGausCop[1] - condGausCop[2]), condGausCop, 
                     max(0, 2 * condGausCop[nx] - condGausCop[nx - 1]))
    int <- sum(diff(c(0, xVals, 1)) * (0.5 * diff(condGausCop) + condGausCop[-(nx + 2)]))
    
    condVineFun <- approxfun(c(0, xVals, 1), condGausCop/int)
    
    condGausCop <- condVineFun(xVals)
    int <- cumsum(c(0, diff(xVals) * (0.5 * diff(condGausCop) + condGausCop[-nx])))
    lower <- max(which(int <= p))
    m <- (condGausCop[lower + 1] - condGausCop[lower])/(xVals[lower + 1] - xVals[lower])
    b <- condGausCop[lower]
    xRes <- -b/m + sign(m) * sqrt(b^2/m^2 + 2 * (p - int[lower])/m)
    predQuantile <- c(predQuantile, margin$q(xVals[lower] + xRes))
  }
  close(pb)
  
  if ("data" %in% slotNames(predLocs)) {
    res <- predLocs
    res@data[[paste("quantile.", p, sep = "")]] <- predQuantile
    return(res)
  }
  else {
    predQuantile <- data.frame(predQuantile)
    colnames(predQuantile) <- paste("quantile.", p, sep = "")
    return(addAttrToGeom(predLocs, predQuantile, 
                         match.ID = FALSE))
  }
}
