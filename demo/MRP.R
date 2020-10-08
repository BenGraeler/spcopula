## get the data
data("simulatedTriples")

## rank order transformation
peakVol <- rankTransform(triples[,1],triples[,3])
colnames(peakVol) <- c("Qp","Vp")
plot(peakVol, asp=1)

# Kendall's tau correlation
cor(triples,method="kendall")

# estiamte the BB7 copula by means of maximum likelihood
copQV <- fitCopula(BB7Copula(param=c(2,14)), peakVol, method="mpl", 
                   start=c(2,14), estimate.variance=F)@copula
copQV

# we use a design return period of 10 years
# the MAR-case: given the 0.9 quantile of one marginal, what is the 
# corresponding quantile of the second variable with the given joint return 
# period of 1/(1-0.9)
v_MAR <- c(0.9, invdduCopula(0.9, copQV, 0.9))
v_MAR
pCopula(v_MAR, copQV)

## the anlytical kendall distribution
kendallFunQV <- getKendallDistr(copQV)

# the Kendall distribution value (fraction of pairs having a smaller copula value than "t")
kendallFunQV(t=0.9)

curve(kendallFunQV, from=0, to=1, asp=1)

# the critical level of the KEN2-RP for 10 years
t_KEN2 <- criticalLevel(kendallFunQV, 10, mu=1)
t_KEN2

# the corresponmding KEN2-RP for the OR-RP
kendallRP(kendallFunQV, cl=0.9, mu=1, copula=copQV)

# illustrating the critical lines (empirically)
contour(copQV,pCopula,levels=c(0.9,t_KEN2),
        xlim=c(0.8,1), ylim=c(0.8,1), n=100, asp=1, col="blue")