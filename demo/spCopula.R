## librarys ##
library("sp")

## meuse - spatial poionts data.frame ##
data("meuse")
coordinates(meuse) = ~x+y

spplot(meuse,"zinc", col.regions=bpy.colors(5))

## margins ##
hist(meuse[["zinc"]],freq=F,n=30,ylim=c(0,0.0035), 
     main="Histogram of zinc", xlab="zinc concentration")

# gevEsti <- fgev(meuse[["zinc"]])$estimate
meanLog <- mean(log(meuse[["zinc"]]))
sdLog <- sd(log(meuse[["zinc"]]))
# curve(dgev(x,gevEsti[1], gevEsti[2], gevEsti[3]),add=T,col="red")
curve(dlnorm(x,meanLog,sdLog),add=T,col="green")

pMar <- function(q) plnorm(q, meanLog, sdLog)
qMar <- function(p) qlnorm(p, meanLog, sdLog)
dMar <- function(x) dlnorm(x, meanLog, sdLog)

meuse$marZinc <- pMar(meuse$zinc)

## lag classes ##
bins <- calcBins(meuse, var="marZinc", nbins=10, cutoff=800)

## calculate parameters for Kendall's tau function ##
# either linear
calcKTauLin <- fitCorFun(bins, degree=1, cutoff=600)
curve(calcKTauLin,0, 1000, col="red",add=TRUE)

# or polynomial (used here)
calcKTauPol <- fitCorFun(bins, degree=3)
curve(calcKTauPol,0, 1000, col="purple",add=TRUE)

## find best fitting copula per lag class
loglikTau <- loglikByCopulasLags(bins, meuse, calcKTauPol,
                                 families=c(normalCopula(0), tCopula(0),
                                            claytonCopula(0), frankCopula(1), 
                                            gumbelCopula(1), joeBiCopula(1.5),
                                            indepCopula()))
bestFitTau <- apply(apply(loglikTau$loglik, 1, rank, na.last=T), 2, 
                    function(x) which(x==7))
colnames(loglikTau$loglik)[bestFitTau]

## set-up a spatial Copula ##
spCop <- spCopula(components=list(normalCopula(0), tCopula(0),
                                  normalCopula(1), tCopula(0), 
                                  claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

## compare spatial copula loglik by lag:
lagData <- lapply(bins$lags, function(x) {
                               as.matrix((cbind(meuse[x[,1], "marZinc"]@data,
                                                meuse[x[,2], "marZinc"]@data)))
                             })

spLoglik <- NULL
for(i in 1:length(bins$lags)) { # i <- 7
  cat("Lag",i,"\n")
  spLoglik <- c(spLoglik,
                sum((dCopula(u=lagData[[i]], spCop,log=T,
                            h=bins$lags[[i]][,3]))))
}

plot(spLoglik, ylab="log-likelihood", xlim=c(1,11)) 
points(loglikTau$loglik[cbind(1:10,bestFitTau)], col="green", pch=16)
points(loglikTau$loglik[,1], col="red", pch=5)
legend(6, 50,c("Spatial Copula", "best copula per lag", "Gaussian Copula",
               "number of pairs"), 
       pch=c(1,16,5,50), col=c("black", "green", "red"))
text(x=(1:10+0.5), y=spLoglik, lapply(lagData,length))

##
# spatial vine
vineDim <- 5L
meuseNeigh <- getNeighbours(meuse,var="marZinc",size=vineDim)

vineCop <- vineCopula(4L)

meuseSpVine <- fitCopula(spVineCopula(spCop, vineCop),
                         list(meuseNeigh, meuse), method="none")

# log-likelihood:
meuseSpVine@loglik

meuseSpVine <- meuseSpVine@copula

##
# leave-one-out x-validation
time <- proc.time()  # ~100 s
predMedian <- NULL
predMean <- NULL
for(loc in 1:nrow(meuseNeigh@data)) { # loc <- 145
  cat("Location:",loc,"\n")
  condSecVine <- condSpVine(condVar=as.numeric(meuseNeigh@data[loc,-1]), 
                          dists=list(meuseNeigh@distances[loc,,drop=F]),meuseSpVine)

  predMedian <- c(predMedian, qMar(optimise(function(x) abs(integrate(condSecVine,0,x)$value-0.5),c(0,1))$minimum))
  
  condExp <-  function(x) {
    condSecVine(pMar(x))*dMar(x)*x
  }
  
  predMean <- c(predMean, integrate(condExp,0,3000,subdivisions=1e6)$value)
}
proc.time()-time

mean(abs(predMean-meuse$zinc))
mean(predMean-meuse$zinc)
sqrt(mean((predMean-meuse$zinc)^2))

mean(abs(predMedian-meuse$zinc))
mean(predMedian-meuse$zinc)
sqrt(mean((predMedian-meuse$zinc)^2))

plot(predMean,meuse$zinc)
abline(0,1)

plot(predMedian,meuse$zinc)
abline(0,1)

## kriging results:
# same neighbourhood size 5L:
# MAE:  158.61
# BIAS:  -4.24
# RMSE: 239.85
#
# global kriging:
# MAE:  148.85
# BIAS:  -3.05
# RMSE: 226.15