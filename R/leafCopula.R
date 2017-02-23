## constructor
leafCopula <- function (param=c(1.446923, -1.722742)) {
  val <- new("leafCopula", dimension = as.integer(2), parameters = param, 
             param.names = c("a3", "a2"), param.lowbnd = c(-Inf, -Inf),
             param.upbnd = c(Inf, Inf), fullname = "leaf copula")
  return(val)
}

# param <- c(1.446923, -1.722742)

# weak lower border, two-place parameter
weakBorderPoly <- function(x, par) {
  par[1]*x^3+par[2]*x^2+x
} 

# ##
# curve(weakBorderPoly,asp=1,ylim=c(0,1),xlim=c(0,1))
# ##

ddxweakBorderPoly <- function(x, par) {
  3*par[1]*x^2+2*par[2]*x+1
}

invWeakBor <- function(v, par) {
  optFun <- function(u) {
    abs(weakBorderPoly(u, par)-v)
  }
  
  optimise(optFun,c(0,1),tol=.Machine$double.eps^0.5)$minimum
}

# strong upper border, no parameter
strongBorderPoly <- function(x) {
  -x^3+x^2+x
}

# ##
# curve(strongBorderPoly,add=T)
# ##

ddxstrongBorderPoly <- function(x) {
  -3*x^2+2*x+1
}

invStrongBor <- function(v) {
  optFun <- function(u) {
    sqrt(mean((strongBorderPoly(u)-v)^2))
  }
  
  optimise(optFun,c(0,1),tol=.Machine$double.eps^0.5)$minimum
}

# precalculate ellipse parameters
solveQ <- function(u) {
  sqrt(0.5)*(strongBorderPoly(u)-u)
}

ddxsolveQ <- function(u) {
  sqrt(0.5)*(ddxstrongBorderPoly(u)-1)
}

# ## double check
# curve(solveQ)
# curve(ddxsolveQ, add=T) 
# abline(h=0)
# abline(solveQ(0.5)-ddxsolveQ(0.5)*0.5,ddxsolveQ(0.5))
# abline(solveQ(0.9)-ddxsolveQ(0.9)*0.9,ddxsolveQ(0.9))
# abline(v=c(0.5,0.9),col="grey")
# ##

solveXb <- function(u, par) {
  sqrt(2)*(u-weakBorderPoly(u, par))+solveQ(u)
}

ddxsolveXb <- function(u, par) {
  wBor <- weakBorderPoly(u, par)
  sqrt(2)*(1-ddxweakBorderPoly(u, par))+ddxsolveQ(u)
}

# ## double check
# curve(solveXb)
# curve(ddxsolveXb, add=T)
# abline(h=0)
# abline(solveXb(0.5)-ddxsolveXb(0.5)*0.5,ddxsolveXb(0.5))
# abline(solveXb(0.9)-ddxsolveXb(0.9)*0.9,ddxsolveXb(0.9))
# abline(v=c(0.5,0.9), col="grey")
# ## okay

###
# ellipse(xb)=b-q:
# b = -q/(sqrt((a^2-x^2)/a^2)-1)
# 
# ellipse'(xb)=1:
# b = (a^2 sqrt((a^2-x^2)/a^2))/x
#
# solve for a:
# -q/(sqrt((a^2-x^2)/a^2)-1) = (a^2 sqrt((a^2-x^2)/a^2))/x
###

solveA <- function(u, par) {
  xb <- solveXb(u, par)
  q <- solveQ(u)
  sqrt(xb^2 - q^2*xb/(2*q-xb))
}

ddusolveA <- function(u, par) {
  xb <- solveXb(u, par)
  q <- solveQ(u)
  dduXb <- ddxsolveXb(u, par)
  dduQ <- ddxsolveQ(u)
  1/(2*solveA(u, par)) * (-2*(q-xb)*(q*xb*(dduQ-3*dduXb)+q^2*dduXb+xb^2*dduXb))/(-2*q + xb)^2
}

# ## double check
# curve(solveA)
# curve(ddusolveA,add=T)
# abline(h=0)
# abline(solveA(0.5)-ddusolveA(0.5)*0.5,ddusolveA(0.5))
# abline(solveA(0.9)-ddusolveA(0.9)*0.9,ddusolveA(0.9))
# abline(v=c(0.5, 0.9), col="grey")
# ## okay

solveB <- function(u, par) {
  a <- solveA(u, par)
  xb <- solveXb(u, par)
  a*sqrt(a^2-xb^2)/xb
} 

ddusolveB <- function(u, par) {
  a <- solveA(u, par)
  xb <- solveXb(u, par)
  dduA <- ddusolveA(u, par)
  dduXb <- ddxsolveXb(u, par)
  (2*a^2*xb*dduA - xb^3*dduA - a^3*dduXb)/(a*xb^2*sqrt(1 - xb^2/a^2)) 
}

# ## double check
# curve(solveB)
# curve(ddusolveB,add=T)
# abline(h=0)
# abline(solveB(0.5)-ddusolveB(0.5)*0.5,ddusolveB(0.5))
# abline(solveB(0.9)-ddusolveB(0.9)*0.9,ddusolveB(0.9))
# abline(v=c(0.5,0.9), col="grey")
# ## okay

ellipse <- function(v,a,b) {
  b*sqrt(1-(v/a)^2)-b
}

# ## double check
# ellipseInV <- function(v) ellipse(v,a,b)
# 
# curve(ellipseInV,0,0.00001,asp=1)
# abline(v=solveXb(0.9)+c(0,-solveQ(0.9)))
# abline(solveXb(0.9)-solveQ(0.9),-1)
# abline(h=0)
# ## okay

dduellipse <- function(v, a, b, dduA, dduB) {
  dduB*sqrt(1-(v/a)^2)+b/(2*sqrt(1-(v/a)^2))*2*v^2/a^3*dduA-dduB
}

# ## double check
# dduellipseInU <- function(u) dduellipse(0.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))
# ellipseInU <- function(u) ellipse(0.5,solveA(u),solveB(u))
# 
# curve(ellipseInU,0.5,1,ylim=c(-0.2,.5))
# curve(dduellipseInU,0.5,1,add=T)
# abline(h=0)
# abline(ellipseInU(0.6)-dduellipseInU(0.6)*0.6,dduellipseInU(0.6))
# abline(ellipseInU(0.9)-dduellipseInU(0.9)*0.9,dduellipseInU(0.9))
# abline(v=c(0.6,0.9), col="grey")
# ## okay

ddvellipse <- function(v, a, b) {
  -b*v/sqrt(1-(v/a)^2)/a^2
}

# ## double check
# curve(ellipseInV,0,0.6,ylim=c(-.2,0))
# curve(ddvellipse,0,0.6,asp=2,add=T)
# abline(h=0)
# abline(ellipseInV(0.2)-ddvellipse(0.2)*0.2,ddvellipse(0.2))
# abline(ellipseInV(0.4)-ddvellipse(0.4)*0.4,ddvellipse(0.4))
# abline(v=c(0.2,0.4), col="grey")
# ## okay

ddvuEllipse <- function(v, a, b, dduA, dduB) {
  (2*a^2*dduA*b*v - a^3*dduB*v - dduA*b*v^3 + a*dduB*v^3)/(a^3*(a^2 - v^2)*sqrt(1 - v^2/a^2))
}

# ## double checking
# curve(dduellipse,ylim=c(-1,.5))
# curve(ddvuEllipse,add=T) #,0,0.6,ylim=c(-.2,0))
# abline(h=0)
# abline(dduellipse(0.2)-ddvuEllipse(0.2)*0.2,ddvuEllipse(0.2))
# abline(dduellipse(0.4)-ddvuEllipse(0.4)*0.4,ddvuEllipse(0.4))
# abline(v=c(0.2,0.4), col="grey")
# ## okay

ddvvEllipse <- function(v, a, b) {
  -b/((a^2 - v^2)*sqrt(1 - v^2/a^2))
}

# ## double checking
# curve(ddvellipse,ylim=c(-1,.5))
# curve(ddvvEllipse,add=T) #,0,0.6,ylim=c(-.2,0))
# abline(h=0)
# abline(ddvellipse(0.2)-ddvvEllipse(0.2)*0.2,ddvvEllipse(0.2))
# abline(ddvellipse(0.4)-ddvvEllipse(0.4)*0.4,ddvvEllipse(0.4))
# abline(v=c(0.2,0.4), col="grey")
# ## okay

vNorm <- function(vo, a, b) {
  (-a^2*b + 2*sqrt(0.5)*a^2*vo + 2*a*b*sqrt(0.25*a^2 + sqrt(0.5)*b*vo - 0.5*vo^2))/(a^2 + b^2)
}

dduvNorm <- function(vo, a, b, dduA, dduB) {
  cRoot <- sqrt(-0.5*vo^2 + 0.25*a^2 + sqrt(0.5)*vo*b) 
  
  low <- a^2+b^2
  high <- -a^2*b+sqrt(2)*a^2*vo+2*a*b*cRoot
  
  dduLow <- 2*a*dduA + 2*b*dduB
  dduHigh <- (-2*a*dduA*b - a^2*dduB + 2*sqrt(2)*vo*a*dduA
              + 2*(dduA*b + a*dduB)*cRoot
              + a*b/cRoot*(0.5*a*dduA+sqrt(0.5)*vo*dduB))

  (low*dduHigh - high*dduLow)/low^2
}

# ## double check
# vNormInU <- function(u) vNorm(.5,solveA(u),solveB(u))
# dduvNormInU <- function(u) dduvNorm(.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))
# 
# curve(vNormInU,0.6,0.85,ylim=c(0.5,.7))
# curve(dduvNormInU,add=T)
# abline(vNormInU(0.65)-dduvNormInU(0.65)*0.65,dduvNormInU(0.65))
# abline(vNormInU(0.8)-dduvNormInU(0.8)*0.8,dduvNormInU(0.8))
# abline(v=c(0.65,0.8), col="grey")
# ## okay

ddvvNorm <- function(vo, a, b) {
  c <- sqrt(0.5)
  (2*a^2*(c + (b^2*(0.5*b*c - 0.5*vo))/sqrt(a^2*b^2*(0.25*a^2 + (b*c - 0.5*vo)*vo))))/(a^2 + b^2)
}

# ## double check
# curve(vNorm,0,0.5,ylim=c(0,1))
# curve(ddvvNorm,add=T)
# abline(vNorm(0.4)-ddvvNorm(0.4)*0.4,ddvvNorm(0.4))
# abline(vNorm(0.1)-ddvvNorm(0.1)*0.1,ddvvNorm(0.1))
# abline(v=c(0.1,0.4), col="grey")
# ## okay

ddvuvNorm <-  function(vo, a, b, dduA, dduB) {
  cRoot <- sqrt(-0.5*vo^2 + 0.25*a^2 + sqrt(0.5)*vo*b) 
  
  low <- a^2+b^2
  
  dduLow <- 2*a*dduA + 2*b*dduB
  ddvHigh <- sqrt(2)*a^2+(a*b*(-vo + sqrt(0.5)*b))/cRoot
  
  ddvuHigh <- (2*sqrt(2)*a*dduA 
               + ((dduA*b + a*dduB)*(-vo+sqrt(0.5)*b)+sqrt(1/2)*a*b*dduB)/cRoot
               - (a*b*(0.5*a*dduA+sqrt(0.5)*dduB*vo)*(-vo+sqrt(0.5)*b))/(2*cRoot^3))
    
  # ddv (low*dduHigh - high*dduLow):
  (low * ddvuHigh - ddvHigh*dduLow)/low^2
}

# ## double check
# curve(dduvNorm,0,0.5,ylim=c(-2,1))
# curve(ddvuvNorm,add=T)
# abline(h=0)
# abline(dduvNorm(0.2)-ddvuvNorm(0.2)*0.2,ddvuvNorm(0.2))
# abline(dduvNorm(0.4)-ddvuvNorm(0.4)*0.4,ddvuvNorm(0.4))
# abline(v=c(0.2,0.4), col="grey")
# ## okay

ddvvvNorm <- function(vo, a, b) {
  (-2*a^4*b^4)/((a^2*b^2*(a^2 + 2*(2*sqrt(0.5)*b - vo)*vo))^(3/2))
}

# ## double check
# curve(ddvvNorm,0,0.5)#,ylim=c(0,1))
# curve(ddvvvNorm,add=T)
# abline(ddvvNorm(0.4)-ddvvvNorm(0.4)*0.4,ddvvvNorm(0.4))
# abline(ddvvNorm(0.1)-ddvvvNorm(0.1)*0.1,ddvvvNorm(0.1))
# abline(v=c(0.1,0.4), col="grey")
# ## okay

yOrig <- function(vn, a, b) {
  sqrt(0.5)*(vn+ellipse(vn, a, b))
}

dduyOrig <- function(vn, a, b, dduA, dduB) {
  sqrt(0.5)*dduellipse(vn, a, b, dduA, dduB)
}

# ## double check
# yOrigInU <- function(u) yOrig(.5,solveA(u),solveB(u))
# dduyOrigInU <- function(u) dduyOrig(.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))
# 
# curve(yOrigInU,0.5,1,ylim=c(0.2,.4))
# curve(dduyOrigInU,add=T)
# abline(yOrigInU(0.6)-dduyOrigInU(0.6)*0.6,dduyOrigInU(0.6))
# abline(yOrigInU(0.9)-dduyOrigInU(0.9)*0.9,dduyOrigInU(0.9))
# abline(v=c(0.6,0.9), col="grey")
# ## okay

ddvyOrig <- function(vn, a, b) {
  sqrt(0.5)*(1+ddvellipse(vn, a, b))
}

# ## double check
# curve(yOrig,0,0.5,ylim=c(0,1))
# curve(ddvyOrig,add=T)
# abline(yOrig(0.475)-ddvyOrig(0.475)*0.475,ddvyOrig(0.475))
# abline(yOrig(0.1)-ddvyOrig(0.1)*0.1,ddvyOrig(0.1))
# abline(v=c(0.1,0.475), col="grey")
# ## okay

ddvuyOrig <- function(vn, a, b, dduA, dduB) {
  sqrt(0.5)*ddvuEllipse(vn, a, b, dduA, dduB)
}

# ## double check
# curve(dduyOrig,0,0.5,ylim=c(-0.2,0.4))
# curve(ddvuyOrig,add=T)
# abline(dduyOrig(0.475)-ddvuyOrig(0.475)*0.475,ddvuyOrig(0.475))
# abline(dduyOrig(0.1)-ddvuyOrig(0.1)*0.1,ddvuyOrig(0.1))
# abline(v=c(0.1,0.475), col="grey")
# ## okay

ddvvyOrig <- function(vn, a, b) {
  sqrt(0.5)*(ddvvEllipse(vn, a, b))
}

# ## double check
# curve(ddvyOrig,0,0.5)#,ylim=c(-0.2,0.4))
# # curve(ddvvyOrig,add=F)
# abline(ddvyOrig(0.475)-ddvvyOrig(0.475)*0.475,ddvvyOrig(0.475))
# abline(ddvyOrig(0.1)-ddvvyOrig(0.1)*0.1,ddvvyOrig(0.1))
# abline(v=c(0.1,0.475), col="grey")
# ## okay

ellipseOrig <- function(vo, a, b) {
  yOrig(vNorm(vo, a, b), a, b)
}

dduellipseOrig <- function(vo, a, b, dduA, dduB) {
  ddvyOrig(vNorm(vo, a, b), a, b)*dduvNorm(vo, a, b, dduA, dduB)+dduyOrig(vNorm(vo, a, b), a, b, dduA, dduB)
}

# ## double check
# ellipseOrigInU <- function(u) ellipseOrig(.5,solveA(u),solveB(u))
# dduellipseOrigInU <- function(u) dduellipseOrig(.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))
# 
# curve(ellipseOrigInU,0.6,0.85,ylim=c(0.2,0.4))
# curve(dduellipseOrigInU,add=T)
# abline(ellipseOrigInU(0.8)-dduellipseOrigInU(0.8)*0.8,dduellipseOrigInU(0.8))
# abline(ellipseOrigInU(0.65)-dduellipseOrigInU(0.65)*0.65,dduellipseOrigInU(0.65))
# abline(v=c(0.8,0.65), col="grey")
# ## okay

ddvellipseOrig <- function(vo, a, b) {
  ddvyOrig(vNorm(vo, a, b), a, b)*ddvvNorm(vo, a, b)
}

# ## double check
# curve(ellipseOrig,0,0.6,ylim=c(0,1))
# curve(ddvellipseOrig,add=T,col="green")
# abline(ellipseOrig(0.4)-ddvellipseOrig(0.4)*0.4,ddvellipseOrig(0.4))
# abline(ellipseOrig(0.1)-ddvellipseOrig(0.1)*0.1,ddvellipseOrig(0.1))
# abline(v=c(0.1,0.4), col="grey")
# abline(h=0)
# ## okay

ddvuEllipseOrig <- function(vo, a, b, dduA, dduB) {
  cvNorm <- vNorm(vo, a, b)
  (ddvuvNorm(vo, a, b, dduA, dduB)*ddvyOrig(cvNorm, a, b)
   + ddvvNorm(vo, a, b)*(dduvNorm(vo, a, b, dduA, dduB)*ddvvyOrig(cvNorm, a, b)
                         + ddvuyOrig(cvNorm, a, b, dduA, dduB)))
}

# ## double check
# curve(dduellipseOrig,0,0.5,ylim=c(-1,1.5))
# curve(ddvuEllipseOrig,add=T)
# abline(dduellipseOrig(0.4)-ddvuEllipseOrig(0.4)*0.4,ddvuEllipseOrig(0.4))
# abline(dduellipseOrig(0.1)-ddvuEllipseOrig(0.1)*0.1,ddvuEllipseOrig(0.1))
# abline(v=c(0.1,0.4), col="grey")
# ## okay

ddvvEllipseOrig <- function(vo, a, b) {
  cvNorm <- vNorm(vo, a, b)
  (ddvvNorm(vo, a, b)^2*ddvvyOrig(cvNorm, a, b)
   +ddvvvNorm(vo, a, b)*ddvyOrig(cvNorm, a, b))
}

# ## double check
# curve(ddvellipseOrig,0,0.5,ylim=c(-1,1.5))
# curve(ddvvEllipseOrig,add=T)
# abline(ddvellipseOrig(0.4)-ddvvEllipseOrig(0.4)*0.4,ddvvEllipseOrig(0.4))
# abline(ddvellipseOrig(0.1)-ddvvEllipseOrig(0.1)*0.1,ddvvEllipseOrig(0.1))
# abline(v=c(0.1,0.4), col="grey")
# ## okay

# cdf 
cdf.leaf <- function(uvabw) {
    uvabw[5]+ellipseOrig(uvabw[2]-uvabw[5], uvabw[3], uvabw[4])
}

pLeafCopula <- function(u, copula) {
  aVec <- solveA(u[,1], copula@parameters)
  bVec <- solveB(u[,1], copula@parameters)
  wBor <- weakBorderPoly(u[,1], copula@parameters)
  sBor <- strongBorderPoly(u[,1])
  
  ret <- apply(u,1,min)
  bool <- (u[,2] > wBor & u[,2] < sBor)
  if(any(bool)) 
    ret[bool] <- apply(cbind(u,aVec,bVec, wBor)[bool,,drop=F], 1, cdf.leaf)
  return(ret)
}

setMethod("pCopula",signature("matrix","leafCopula"), pLeafCopula)
setMethod("pCopula",signature("numeric","leafCopula"), 
          function(u, copula) pLeafCopula(matrix(u,ncol=copula@dimension), copula))

# partial derivative d/dv
ddv.cdf.leaf <- function(uvabw) {
  ddvellipseOrig(uvabw[2]-uvabw[5], uvabw[3], uvabw[4])
}

ddvLeafCopula <- function(u, copula) {
  aVec <- solveA(u[,1], copula@parameters)
  bVec <- solveB(u[,1], copula@parameters)
  wBor <- weakBorderPoly(u[,1], copula@parameters)
  sBor <- strongBorderPoly(u[,1])
  
  ret <- numeric(nrow(u))
  ret[u[,2] <= wBor] <- 1
  bool <- (u[,2] > wBor & u[,2] < sBor)
  if(any(bool)) 
    ret[bool] <- apply(cbind(u,aVec,bVec, wBor)[bool,,drop=F], 1, ddv.cdf.leaf)
  return(ret)
}

setMethod("ddvCopula",signature("matrix","leafCopula"), ddvLeafCopula)
setMethod("ddvCopula",signature("numeric","leafCopula"), 
          function(u, copula) ddvLeafCopula(matrix(u,ncol=copula@dimension), copula))

# partial derivative d/du
ddu.cdf.leaf <- function(uvabw, par) {
  (ddxweakBorderPoly(uvabw[1], par)*(1-ddvellipseOrig(uvabw[2]-uvabw[5], uvabw[3], uvabw[4]))
   + dduellipseOrig(uvabw[2]-uvabw[5], uvabw[3], uvabw[4], ddusolveA(uvabw[1], par), ddusolveB(uvabw[1], par)))
}

dduLeafCopula <- function(u, copula) {
  aVec <- solveA(u[,1], copula@parameters)
  bVec <- solveB(u[,1], copula@parameters)
  wBor <- weakBorderPoly(u[,1], copula@parameters)
  sBor <- strongBorderPoly(u[,1])
  
  ret <- numeric(nrow(u))
  ret[u[,2] >= sBor] <- 1
  bool <- (u[,2] > wBor & u[,2] < sBor)
  if(any(bool)) 
    ret[bool] <- apply(cbind(u,aVec,bVec, wBor)[bool,,drop=F], 1, ddu.cdf.leaf, 
                       par=copula@parameters)
  return(ret)
}

setMethod("dduCopula",signature("matrix","leafCopula"), dduLeafCopula)
setMethod("dduCopula",signature("numeric","leafCopula"), 
          function(u, copula) dduLeafCopula(matrix(u, ncol=copula@dimension), copula))

# density
pdf.leaf <- function(uvabw, par) {
  (ddvuEllipseOrig(uvabw[2]-uvabw[5], uvabw[3], uvabw[4], 
                   ddusolveA(uvabw[1], par), ddusolveB(uvabw[1], par))
   - ddxweakBorderPoly(uvabw[1], par)*ddvvEllipseOrig(uvabw[2]-uvabw[5],
                                                       uvabw[3], uvabw[4]))
}

dLeafCopula <- function(u, copula) {
  aVec <- solveA(u[,1], copula@parameters)
  bVec <- solveB(u[,1], copula@parameters)
  wBor <- weakBorderPoly(u[,1], copula@parameters)
  sBor <- strongBorderPoly(u[,1])
  
  ret <- numeric(nrow(u))
  bool <- (u[,2] > wBor & u[,2] < sBor)
  if(any(bool))
    ret[bool] <- apply(cbind(u,aVec,bVec,wBor)[bool,,drop=F],1,pdf.leaf,
                       par=copula@parameters)
  return(ret)
}

setMethod("dCopula",signature("matrix","leafCopula"), dLeafCopula)
setMethod("dCopula",signature("numeric","leafCopula"), 
          function(u, copula) dLeafCopula(matrix(u,ncol=copula@dimension), copula))

## random pair generator

# inverse of partial derivative d/dv
invddvLeafCop <- function(vy, par) {
  
  optFun <- function(u) {
    y <- ddv.cdf.leaf(c(u,vy[1], solveA(u, par), solveB(u, par), weakBorderPoly(u, par)))
    ret <- abs(y-vy[2])
    if(is.nan(ret) | is.infinite(ret)) {
      cat("The warning occured for the points: \n")
      cat("u:", u,"\n")
      cat("y:", y,"\n")
      cat("vy", vy,"\n")
    }
    return(ret)
  }
  
  invWBor <- invWeakBor(vy[1], par)
  invSBor <- invStrongBor(vy[1])
  
  if(invSBor == invWBor)
    return(invSBor)
  optimise(optFun,c(min(invSBor,invWBor),max(invSBor,invWBor)))$minimum
}

r.leaf <- function(n, copula) {
  v <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
    
  res <- cbind(apply(cbind(v,y), 1, invddvLeafCop, par=copula@parameters), v)
  colnames(res) <- c("u","v")
    
  return(res)
}

setMethod("rCopula",signature("numeric","leafCopula"),r.leaf)

# ## example section
# leafCop <- leafCopula()
# plot(rCopula(500,leafCop))
# 
# persp(leafCop,dCopula,zlim=c(0,10))
# persp(leafCop,pCopula)