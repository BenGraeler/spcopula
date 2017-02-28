## Truncated copulas exhibting a crisp boundary, often induced by lower bounds.
## Points below the boundary are shifted "upwards" onto the boundary. Hence, 
## considerable mass is concentrated on the boundary yielding a mixed density 
## analougously to mixed discrete continuous distributions in the univariate 
## case.

# class truncated copula
validTrunCop <- function(object) {
  if(any(object@trunFamily@parameters != object@parameters[-c(length(object@parameters)-(1:0))])) {
    warning("Missmatch of parameters between the parameter slot and the parameter slot of the \"trunFamily\".")
  }
  
  ifelse(object@dimension == 2, TRUE, FALSE)
}

# Slots:
#   
#   Name:    trunFamily      contPar       .tools    dimension   parameters  param.names param.lowbnd  param.upbnd
# Class:       copula      numeric         list      integer      numeric    character      numeric      numeric
# 
# Name:      fullname
# Class:    character

setClass("trunCopula", 
         list("copula", trunFamily = "copula", contPar = "numeric", .tools = "list"),
         validity = validTrunCop,
         contains = list("copula"))

trunCopula <- function(copula, contPar, approx.u=1:1000/1000) {
  
  # setting helper functions
  contFun <- function(x) x^contPar
  invContFun <- function(x) x^(1/contPar)
  
  trunFun <- approxfun(c(0, approx.u), 
                       c(0, qCopula_u(copula, contFun(approx.u), approx.u)[,2]))

  invTrunFun <- approxfun(trunFun(c(0, approx.u)), c(0, approx.u))
  
  CDF <-  approxfun(c(0, approx.u), 
                    c(0, pCopula(cbind(invTrunFun(approx.u), approx.u), copula)))
  invCDF <- approxfun(CDF(c(0, approx.u)), c(0, approx.u))

  # calculate density along the contour line
  dCont <- function(u) {
    v <- trunFun(u)
    (dduCopula(cbind(u,v), copula) - dduCopula(cbind(u,0), copula))
  }
  
  new("trunCopula", 
      dimension = dim(copula),
      parameters = c(copula@parameters, contPar),
      param.names = c(copula@param.names, "truncation"),
      param.lowbnd = c(copula@param.lowbnd, -Inf),
      param.upbnd = c(copula@param.upbnd, Inf),
      fullname = "truncated copula",
      trunFamily = copula,
      contPar = contPar, 
      .tools = list(trunFun = trunFun,
                    invTrunFun = invTrunFun,
                    CDF = CDF,
                    invCDF = invCDF,
                    contFun = contFun,
                    invContFun = invContFun,
                    dCont = dCont))
}

## console printing
setMethod("describeCop", c("trunCopula", "character"),
          function(x, kind = c("short", "very short", "long"), prefix = "", ...) {
            kind <- match.arg(kind)
            if(kind == "very short") # e.g. for show() which has more parts
              return(paste0(prefix, "truncated copula"))
            
            name <- paste("truncated", describeCop(x@trunFamily, "very short"))
            d <- dim(x)
            ch <- paste0(prefix, name, ", dim. d = ", d)
            switch(kind <- match.arg(kind),
                   short = ch,
                   long = paste0(ch, "\n", prefix, " param.: ",
                                 capture.output(str(x@parameters,
                                                    give.head=FALSE))),
                   stop("invalid 'kind': ", kind))
          })

## density

dTrunCop <- function(u, copula, log=FALSE, ..., tol=1e-3) {
  if (log) {
    res <- rep(NA, nrow(u))
  } else {
    res <- rep(0, nrow(u))
  }
  
  contVals <- copula@.tools$contFun(u[,1])
  diffContVals <- u[,2] - contVals
  
  # split in above and on contour
  boolAbove <- diffContVals >= tol
  boolContour <- abs(diffContVals) < tol
  
  # shift back
  u[,2] <- sapply(u[,2], function(v) copula@.tools$invCDF(v))
  
  res[boolAbove] <- dCopula(u[boolAbove,], copula@trunFamily, log, ...)
  
  if (any(boolContour)) {
    res[boolContour] <- copula@.tools$dCont(u[boolContour,1])
    if (log)
      res[boolContour] <- log(res[boolContour])
  }
  
  return(res)
}

# setMethod(dCopula, c("matrix", "trunCopula"), dTrunCop)
# 
# setMethod(dCopula, c("numeric", "trunCopula"), 
#           function(u, copula, log, ...) {
#             dTrunCop(matrix(u, ncol=2), copula, log, ...)
#           })

## sampling from the trunCopula

rTrunCop <- function(n, copula, ...) {
  smpl <- rCopula(n, copula@trunFamily, ...)
  smpl[,2] <- pmax(copula@.tools$CDF(smpl[,2]),
                   copula@.tools$contFun(smpl[,1]))
  
  return(smpl)
}

setMethod(rCopula, c("numeric", "trunCopula"), rTrunCop)

## CDF of the trunCopula

pTrunCop <- function(u, copula, ...) {
  res <- u[,1]
  boolu11 <- u[,1] == 1
  res[boolu11] <- u[boolu11,2]
  
  boolu21 <- u[,2] == 1
  res[boolu21] <- u[boolu21,1]
  
  contVals <- copula@.tools$contFun(u[,1])
  boolBelow <- u[,2] < contVals

  u[boolBelow, 1] <- copula@.tools$invContFun(u[boolBelow,2])
  
  u[,2] <- copula@.tools$invCDF(u[,2])
  
  res[!(boolu11 | boolu21)] <- pCopula(u[!(boolu11 | boolu21),], copula@trunFamily)# , ...)
  return(res)
}

setMethod(pCopula, c("numeric", "trunCopula"), 
          function(u, copula, ...) pTrunCop(matrix(u, ncol = dim(copula)), copula, ...))

setMethod(pCopula, c("matrix", "trunCopula"), pTrunCop)

### CDF version ###
fitTrunCop <- function(copula, data, ..., method, lower, upper, tol=1e-3) {
  if (missing(method))
    method <- ifelse(length(copula@trunFamily@parameters) > 1, "Nelder-Mead", "Brent")
  if (missing(lower))
    lower <- ifelse(is.infinite(copula@trunFamily@param.lowbnd), -1e3, copula@trunFamily@param.lowbnd)
  if (missing(upper))
    upper <- ifelse(is.infinite(copula@trunFamily@param.upbnd), 1e3, copula@trunFamily@param.upbnd)
  
  pEmpCop <- pCopula(data, empiricalCopula(data))
  
  optFun <- function(par) {
    cat(par, "\n")
    innerCop <- copula@trunFamily
    innerCop@parameters <- par
    cop <- trunCopula(innerCop, copula@contPar)
    
    mae <- mean(abs(pCopula(data, cop) - pEmpCop))
    cat(mae, "\n")
    mae
  }
  
  optOut <- optim(copula@trunFamily@parameters, optFun, 
                  method = method, lower = lower, upper = upper, ...)
  
  innerCop <- copula@trunFamily
  innerCop@parameters <- optOut$par
  cop <- trunCopula(innerCop, copula@contPar)
  
  new("fitCopula", 
      copula=cop, 
      estimate = c(optOut$par, copula@contPar),
      var.est = matrix(NA),
      loglik = sum(dCopula(data, cop, log=T, tol=tol)),
      nsample = as.integer(nrow(data)),
      method = "Copula CDF optimisation with fixed boundary.",
      call = match.call(),
      fitting.stats = optOut)
}

setMethod("fitCopula", c("trunCopula", "matrix"), fitTrunCop)

# ## sample along contour
rTrunCop_y <- function(y, copula, n=1, n.disc = 1000) {
  stopifnot(copula@dimension == 2)
  n.y <- length(y)
  stopifnot(n.y == 1 | n == 1)
  
  uIntSec <- copula@.tools$invContFun(y)

  smpl <- matrix(NA, n.y*n, 2)
  
  for (i in 1:n.y) { # i <- 1 i <- i+1
    condVals <- seq(y[i], 1-(1-y[i])/n.disc^2, length.out = n.disc)
    uv <- qCopula_v(copula, rep(y[i], n.disc-1), condVals[-1])
    uv <- rbind(uv, qCopula_u(copula, rep(y[i], sum(condVals < uIntSec[i])),
                              condVals[condVals < uIntSec[i]]))
    uv <- uv[order(uv[,1]),]
    
    dSeq <- cumsum(c(0, apply((uv[-nrow(uv),]-uv[-1,])^2, 1, function (x) sqrt(sum(x)))))
    probs <- dTrunCop(uv, copula)
    
    apFun <- approxfun(dSeq, probs, rule = 2)
    probLine <- copula@.tools$dCont(uIntSec[i])
    probCont <- integrate(apFun, 0, max(dSeq))$value
    
    rContour <- runif(n, 0, probCont + probLine)
    
    funAppConPoint <- function(rCont) {
      invCDFContour <- function(x) {
        abs(integrate(apFun, 0, x)$value - rCont)
      } 
      
      lContour <- optimise(invCDFContour, c(0, max(dSeq)))$minimum
      
      dSeqInt <- findInterval(lContour, dSeq)
      
      lSeq <- sqrt(sum((uv[dSeqInt,]-uv[dSeqInt+1,])^2))
      
      uv[dSeqInt,] + (lContour - dSeq[dSeqInt])/lSeq * (uv[dSeqInt+1,]-uv[dSeqInt,])
    }
  
    if (n == 1) {
      if (rContour <= probLine) {
        smpl[i,] <- c(uIntSec[i], y[i])
        next;
      }
      rContour <- rContour - probLine
      
      appConPoint <- funAppConPoint(rContour)
      
      if (appConPoint[1] > appConPoint[2]) {
        smpl[i,] <- qCopula_u(copula, y[i], appConPoint[1])
      } else {
        smpl[i,] <- qCopula_v(copula, y[i], appConPoint[2])
      }
    } else {
      boolLine <- rContour <= probLine
      smpl <- cbind(rep(uIntSec, n),
                           rep(y, n))
      rContour <- rContour - probLine
      smpl[!boolLine,] <- t(sapply(rContour[!boolLine], funAppConPoint))
      
      boolLower <- smpl[,1] > smpl[,2]
      if (any(boolLower & !boolLine))
        smpl[boolLower & !boolLine,] <- qCopula_u(copula, rep(y[i], sum(boolLower & !boolLine)),
                                                  smpl[boolLower & !boolLine, 1])
      if (any(!boolLower & !boolLine))
        smpl[!boolLower & !boolLine,] <- qCopula_v(copula, rep(y[i], sum(!boolLower & !boolLine)),
                                                   smpl[!boolLower & !boolLine, 2])
    }
  }
  
  return(smpl)
}

setMethod(rCopula_y, signature = c("numeric", "trunCopula"), rTrunCop_y)

## cond inverse

qTrunCop_v <- function(copula, p, v, tol=.Machine$double.eps^.5) { # sample=NULL
  stopifnot(length(p) == length(v)) 

  cbind(sapply(1:length(p), 
                function(ind) {
                  if (v[ind] < p[ind]) 
                    return(NA)
                  if (v[ind] == 1)
                    return(p[ind])
                  optimise(function(u) abs(pCopula(cbind(u, v[ind]), copula) - p[ind]),
                           c(p[ind], copula@.tools$invContFun(v[ind])), tol=tol)$minimum
                }), v)
}

setMethod("qCopula_v", signature = c("trunCopula"), qTrunCop_v)
