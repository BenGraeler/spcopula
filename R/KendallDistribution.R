# derivation of kendall distributions in higher dimensions for Archimedean copulas

###########
## Frank ##
###########

# generator
genFrank <- function(t, theta) copFrank@iPsi(t, theta)
# -log( (exp(-theta*u)-1) / (exp(-theta)-1) )

# use series expansion for small u?
d1genFrank <- function(t, theta) {
  theta / (1 - exp(theta * t))
}

d2genFrank <- function(t, theta) {
  (theta^2 * exp(theta * t)) / (1 - exp(theta * t))^2
}

d3genFrank <- function(t, theta) {
  -(theta^3 * exp(theta*t) * (exp(theta * t) + 1))/(exp(theta * t)-1)^3
}

## inverse generator
igenFrank <- function(s, theta) copFrank@psi(s, theta)
# -log(1-(1-exp(-theta))*exp(-t))/theta

d1igenFrank <- function(s, theta) {
  eth <- exp(theta)
  (1 - eth) / (theta * (-eth + 1 + exp(theta + s)))
}

d2igenFrank <- function(s, theta) {
  eth <- exp(theta)
  eths <- exp(theta + s)
  ((eth-1) * eths) / (theta * (-eth + eths + 1)^2)
}

d3igenFrank <- function(s, theta) {
  eth <- exp(theta)
  eths <- exp(theta + s)
  ((eth - 1) * eths) / (theta * (-eth + eths + 1)^2)-(2*(eth-1) * exp(2 * theta + 2 * s))/(theta * (-eth + eths + 1)^3)  
}

kdFrank <- function(t, cop) {
  stopifnot(cop@dimension <=4)
  .theta <- cop@parameters
  .val <- 0 < t & t < 1
  
  gt <- genFrank(t[.val], .theta)
  sum1 <- 0
  for (i in 1:(cop@dimension-1)) {
    digen <- switch(i,
                    d1igenFrank(gt, .theta),
                    d2igenFrank(gt, .theta),
                    d3igenFrank(gt, .theta))
    sum1 <- sum1 + (-1)^i/factorial(i) * gt^i * digen
  }
  
  res <- t
  res[.val] <- res[.val] + sum1
  res
}

setMethod("kendall", signature("numeric", "frankCopula"), function(t, copula) kdFrank(t, copula))

###################
## Gumbel Copula ##
###################
# generator
genGumbel <- function(t, theta) copGumbel@iPsi(t, theta)
# -log(t)^theta

# d1genGumbel <- function(t, theta) {}
# d2genGumbel <- function(t, theta) {}
# d3genGumbel <- function(t, theta) {}

## inverse generator
igenGumbel <- function(s, theta) copGumbel@psi(s, theta)
# exp(-s^(1/theta))

d1igenGumbel <- function(s, theta) {
  -(exp(-s^(1/theta)) * s^(1/theta-1))/theta
}

d2igenGumbel <- function(s, theta) {
  s1th <- s^(1/theta)
  (exp(-s1th) * s^(1/theta-2) * (theta+s1th-1))/theta^2
}

d3igenGumbel <- function(s, theta) {
  s1th <- s^(1/theta)
  ems1th <- exp(-s1th)
  s2th3 <- s^(2 / theta - 3)
  -(ems1th * (theta + s1th - 1) * s2th3) / (theta^3 + ems1th * s2th3) / theta^3 + ( (1 / theta - 2) * ems1th * (theta + s1th-1) * s^(1 / theta - 3)) / theta^2
}

kdGumbel <- function(t, cop) {
  stopifnot(cop@dimension <=4)
  .theta <- cop@parameters
  .val <- 0 < t & t < 1
  
  gt <- genGumbel(t[.val], .theta)
  sum1 <- 0
  for (i in 1:(cop@dimension-1)) {
    digen <- switch(i,
                    d1igenGumbel(gt, .theta),
                    d2igenGumbel(gt, .theta),
                    d3igenGumbel(gt, .theta))
    sum1 <- sum1 + (-1)^i/factorial(i) * gt^i * digen
  }
  
  res <- t
  res[.val] <- res[.val] + sum1
  res
}

setMethod("kendall", signature("numeric", "gumbelCopula"), function(t, copula) kdGumbel(t, copula))

#############
## Clayton ##
#############

# generator
genClayton <- function(t, theta) copClayton@iPsi(t, theta)
# u^(-theta) - 1

# d1genClayton <- function(t, theta) {}
# d2genClayton <- function(t, theta) {}
# d3genClayton <- function(t, theta) {}

## inverse generator
igenClayton <- function(s, theta) copClayton@psi(s, theta)
# (1 + t)^(-1/theta)

d1igenClayton <- function(s, theta) {
  -(s+1)^(-(theta+1)/theta)/theta
}

d2igenClayton <- function(s, theta) {
  ((theta + 1) * (s + 1)^(-1 / theta - 2)) / theta^2
}

d3igenClayton <- function(s, theta) {
  -((theta + 1) * (2 * theta + 1) * (s + 1)^(-1 / theta - 3)) / theta^3
}

kdClayton <- function(t, cop) {
  stopifnot(cop@dimension <=4)
  .theta <- cop@parameters
  .val <- 0 < t & t < 1
  
  gt <- genClayton(t[.val], .theta)
  sum1 <- 0
  for (i in 1:(cop@dimension-1)) {
    digen <- switch(i,
                    d1igenClayton(gt, .theta),
                    d2igenClayton(gt, .theta),
                    d3igenClayton(gt, .theta))
    sum1 <- sum1 + (-1)^i/factorial(i) * gt^i * digen
  }
  
  res <- t
  res[.val] <- res[.val] + sum1
  res
}

setMethod("kendall", signature("numeric", "claytonCopula"), function(t, copula) kdClayton(t, copula))

#########
## Joe ##
#########


# generator
genJoe <- function(t, theta) copJoe@iPsi(t, theta)
# -log(1-(1 - t)^theta)

# d1genClayton <- function(t, theta) {}
# d2genClayton <- function(t, theta) {}
# d3genClayton <- function(t, theta) {}

## inverse generator
igenJoe <- function(s, theta) copJoe@psi(s, theta)
# 1 - (1 - exp(-s))^(1 / theta)

d1igenJoe <- function(s, theta) {
  ( -expm1(-s))^(1/theta)/(theta-theta * exp(s))
}

d2igenJoe <- function(s, theta) {
  exps <- exp(s)
  expm1ms <- expm1(-s)
  expm1s <- expm1(s)
  
  ((-expm1ms)^(1/theta) * (theta * exps - 1))/(theta * expm1s)^2
}

# (e^(-s) (1-e^(-s))^(1/theta-1))/(theta (theta-theta e^s))
# (theta e^s (1-e^(-s))^(1/theta))/(theta-theta e^s)^2

d3igenJoe <- function(s, theta) {
  exps <- exp(s)
  expm1ms <- expm1(-s)
  expm1s <- expm1(s)
  
  ds1 <- (-expm1ms)^(1 / theta) * (2 * theta * exps - 1)
  ds2 <- -theta * (exps * (-expm1ms)^(1/theta) * (theta + theta * exps - 1))
  
  (ds1 + ds2) / (theta * expm1s)^3 
}

kdJoe <- function(t, cop) {
  stopifnot(cop@dimension <= 4)
  .theta <- cop@parameters
  .val <- 0 < t & t < 1
  
  gt <- genJoe(t[.val], .theta)
  sum1 <- 0
  for (i in 1:(cop@dimension-1)) {
    digen <- switch(i,
                    d1igenJoe(gt, .theta),
                    d2igenJoe(gt, .theta),
                    d3igenJoe(gt, .theta))
    sum1 <- sum1 + (-1)^i/factorial(i) * gt^i * digen
  }
  
  res <- t
  res[.val] <- res[.val] + sum1
  res
}

setMethod("kendall", signature("numeric", "joeCopula"), function(t, copula) kdJoe(t, copula))
setMethod("kendall", signature("numeric", "joeBiCopula"), function(t, copula) kdJoe(t, copula))