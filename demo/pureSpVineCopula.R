## librarys ##
library("VineCopula")
library("sp")
par(mfrow=c(1,1))
## meuse - spatial poionts data.frame ##
data("meuse")
coordinates(meuse) = ~x+y

spplot(meuse,"zinc", col.regions=bpy.colors(5))

## margins ##
hist(meuse[["zinc"]],freq=F,n=30,ylim=c(0,0.0035), 
     main="Histogram of zinc", xlab="zinc concentration")

 #gevEsti <- fgev(meuse[["zinc"]])$estimate
meanLog <- mean(log(meuse[["zinc"]]))
sdLog <- sd(log(meuse[["zinc"]]))
# curve(dgev(x,gevEsti[1], gevEsti[2], gevEsti[3]),add=T,col="red")
curve(dlnorm(x,meanLog,sdLog),add=T,col="green")

pMar <- function(q) plnorm(q, meanLog, sdLog)
qMar <- function(p) qlnorm(p, meanLog, sdLog)
dMar <- function(x) dlnorm(x, meanLog, sdLog)

# pMar <- function(q) pgev(q, gevEsti[1], gevEsti[2], gevEsti[3])
# qMar <- function(p) qgev(p, gevEsti[1], gevEsti[2], gevEsti[3])
# dMar <- function(x) dgev(x, gevEsti[1], gevEsti[2], gevEsti[3])

meuse$rtZinc <- rank(meuse$zinc)/(length(meuse)+1)
hist(meuse$rtZinc)
## lag classes ##
bins <- calcBins(meuse, var="rtZinc", nbins=10, cutoff=800)


## calculate parameters for Kendall's tau function ##
calcKTau <- fitCorFun(bins, degree=2)
curve(calcKTau,0, 1000, col="purple",add=T)

families <- list(normalCopula(0.3), tCopula(0.3,df=2.15), claytonCopula(0.3),
                 gumbelCopula(2), frankCopula(1), joeBiCopula(1.5),
                 surClaytonCopula(1), surGumbelCopula(1), surJoeBiCopula(1.5))

## find best fitting copula per lag class
loglikTau <- loglikByCopulasLags(bins, meuse, families, calcKTau)
bestFitTau <- apply(apply(loglikTau$loglik, 1, rank, na.last=T), 2, 
                    function(x) which(x==9))
colnames(loglikTau$loglik)[bestFitTau]

## set up a first bivariate spatial Copula
###########################################
spCop <- spCopula(c(families[bestFitTau[1:8]],indepCopula()),
                  distances=bins$meanDists[1:9],
                  spDepFun=calcKTau, unit="m")

## estimation neighbourhood for the pure spatial vine copula
#############################################################
vineDim <- 20L
meuseNeigh1 <- getNeighbours(dataLocs=meuse,var="rtZinc",size=vineDim)

## second spatial tree
#######################
meuseNeigh2 <- dropSpTree(meuseNeigh1, meuse, spCop)
bins2 <- calcBins(meuseNeigh2, "rtZinc", boundaries=c(0,2:15)*50, plot=F)
points(bins2$meanDists, bins2$lagCor, pch=2)
calcKTau2 <- fitCorFun(bins2, degree=3,cutoff=500)
curve(calcKTau2,0, 800, col="green",add=TRUE)

loglikTau2 <- loglikByCopulasLags(bins2, families = families, calcCor =  calcKTau2)
bestFitTau2 <- apply(apply(loglikTau2$loglik, 1, rank, na.last=T), 2, 
                    function(x) which(x==9))
colnames(loglikTau2$loglik)[bestFitTau2]

## set up the second bivariate spatial Copula
##############################################
spCop2 <- spCopula(c(families[bestFitTau2[1:6]], indepCopula()),
                   distances=bins2$meanDists[1:7],
                   spDepFun=calcKTau2, unit="m")

## third spatial tree
######################
meuseNeigh3 <- dropSpTree(meuseNeigh2, meuse, spCop2)
bins3 <- calcBins(meuseNeigh3, "rtZinc", plot=F)
points(bins3$meanDists, bins3$lagCor, pch=3)
calcKTau3 <- fitCorFun(bins3, degree=1, cutoff=500)
curve(calcKTau3, 0, 500, col="red", add=TRUE)

loglikTau3 <- loglikByCopulasLags(bins3, families = families, calcCor = calcKTau3)
bestFitTau3 <- apply(apply(loglikTau3$loglik, 1, rank, na.last=T), 2, 
                     function(x) which(x==9))
colnames(loglikTau3$loglik)[bestFitTau3]

## set up the third bivariate spatial Copula
#############################################
spCop3 <- spCopula(c(families[bestFitTau3[1:5]],indepCopula()),
                   distances=bins3$meanDists[1:6],
                   spDepFun=calcKTau3, unit="m")

## fourth spatial tree
#######################
meuseNeigh4 <- dropSpTree(meuseNeigh3, meuse, spCop3)
bins4 <- calcBins(meuseNeigh4, "rtZinc", plot=F)
points(bins4$meanDists, bins4$lagCor, pch=4)
calcKTau4 <- fitCorFun(bins4, degree=1,cutoff=400)
curve(calcKTau4,0, 500, col="blue",add=TRUE)

legend("topright",c("1. spatial cop.", "2. spatial cop.",
                    "3. spatial cop.", "4. spatial cop."),
       pch=1:4,col=c("purple","green","red","blue"),lty=1)

loglikTau4 <- loglikByCopulasLags(bins4, families = families,calcCor =  calcKTau4)
bestFitTau4 <- apply(apply(loglikTau4$loglik, 1, rank, na.last=T), 2, 
                     function(x) which(x==9))
colnames(loglikTau4$loglik)[bestFitTau4]

## set up the fourth bivariate spatial Copula
#############################################
spCop4 <- spCopula(c(families[bestFitTau4[1:3]], normalCopula(0),indepCopula()),
                   distances=bins4$meanDists[1:5],
                   spDepFun=calcKTau4, unit="m")

## pure spatial vine
#####################
meuseSpVine <- spVineCopula(list(spCop, spCop2, spCop3, spCop4))

## neighbourhood for cross-validation using a 5 dim pure spatial vine copula
#############################################################################
vineDim <- 5L

meuse$lnZinc <- pMar(meuse$zinc)
meuseNeigh <- getNeighbours(dataLocs=meuse, predLocs=meuse, prediction=T, 
                            min.dist=10, var="lnZinc", size=vineDim)

# meuse$evZinc <- pMar(meuse$zinc)
# meuseNeigh <- getNeighbours(dataLocs=meuse, predLocs=meuse, prediction=T, 
#                             min.dist=10, var="evZinc", size=vineDim)

## leave-one-out x-validation
##############################

time <- proc.time()  # ~160 s
predMedian <- spCopPredict(meuseNeigh, meuse, meuse, meuseSpVine, margin=list(q=qMar), "quantile")
predMean <- spCopPredict(meuseNeigh, meuse, meuse, meuseSpVine, margin=list(q=qMar), "expectation")
proc.time() - time

c(mean(abs(predMean$expect-meuse$zinc)),
  mean(predMean$expect-meuse$zinc),
  sqrt(mean((predMean$expect-meuse$zinc)^2)))

c(mean(abs(predMedian$quantile.0.5-meuse$zinc),na.rm = T),
  mean(predMedian$quantile.0.5-meuse$zinc, na.rm = T),
  sqrt(mean((predMedian$quantile.0.5-meuse$zinc)^2, na.rm = T)))

par(mfrow=c(1,2))
plot(predMean$expect, meuse$zinc,asp=1)
abline(0,1)

plot(predMedian$quantile.0.5, meuse$zinc,asp=1)
abline(0,1)

boxplot(predMean@data[c("zinc","expect")])
boxplot(predMedian@data[c("zinc","quantile.0.5")])

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