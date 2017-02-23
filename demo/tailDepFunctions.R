library("VineCopula")
library("spcopula")
data("simulatedTriples")

rtPair <- 1-as.matrix(rankTransform(triples[,c(1,3)]))

plot(rtPair,asp=1)

tdfEmp <- empBivJointDepFun(rtPair)
plot(tdfEmp,ylim=c(0,1), ylab="tail index", xlab="u")
abline(v=0.5, col="grey")

gaussCop <- fitCopula(normalCopula(0), rtPair)@copula
tdfGauss <- bivJointDepFun(gaussCop)
curve(tdfGauss, add=T,col="green",n=500)

gumbelCop <- fitCopula(gumbelCopula(2),rtPair)@copula
tdfGumbel <- bivJointDepFun(gumbelCop)
curve(tdfGumbel,add=T, col="blue",n=500)

BB6Cop <- fitCopula(BB6Copula(), rtPair)@copula
tdfBB6 <- bivJointDepFun(BB6Cop)
curve(tdfBB6, add=T,col="red",n=500)

legend("bottomright",
       c("empirical", "Gaussian", "Gumbel", "BB6"),
       col=c("black", "green", "blue", "red"), lty=1)