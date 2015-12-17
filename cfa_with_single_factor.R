### CFA with a single factor

## Fit values for a CFA with a single factor can only be computed
## with models with  more than 3 observed variables.
## With 3 obeserved variables and one factor you have to compute 6 free
## parameters (with one fixed). But your covariance matrix has also 6
## nonredundant elements. Hence you get df=0 degress of freedom.
## A chi-squared distribution with df=0 is constant 0, hence we
## can expect the test statistic to be 0.

library(lavaan)

## data/model/fit with 3 and with 4 variables
data1 <- HolzingerSwineford1939[, 7:9]
data2 <- HolzingerSwineford1939[, 7:10]

model1 <- 'trait =~ x1 + x2 + x3'
model2 <- 'trait =~ x1 + x2 + x3 + x4'

fit1 <- cfa(model1, data=data1)
fit2 <- cfa(model2, data=data2)

## observed covariance matrix (with denominator n)
s1 <- cov(data1) * (nrow(data1)-1)/nrow(data1)
s2 <- cov(data2) * (nrow(data2)-1)/nrow(data2)

## estimated model implied covariance matrix
si1 <- fit1@implied$cov[[1]]
si2 <- fit2@implied$cov[[1]]

## number of observed variables
p1 <- ncol(data1)
p2 <- ncol(data2)

## discrepancy function for the ML estimation procedure (Jöreskog, 1967)
F.hat <- function(s, si, p) {
    log(det(si)) + sum(diag(s %*% solve(si))) - log(det(s)) - p
}

F.hat(s1, si1, p1)
F.hat(s2, si2, p2)

## test statistics are chi-square distributed
(t1 <- F.hat(s1, si1, p1) * nrow(data1))
(t2 <- F.hat(s2, si2, p2) * nrow(data2))

## compare with lavaan output
t1.l <- fit1@test[[1]]$stat
t2.l <- fit2@test[[1]]$stat

abs(t1 - t1.l)
abs(t2 - t2.l)

## cfi fit measure
cfi <- function(chi2.base, df.base, chi2.user, df.user) {
    d.base <- chi2.base - df.base
    d.user <- chi2.user - df.user
    cfi <- (d.base - d.user) / d.base
    names(cfi) <- "cfi"
    return(cfi)
}

base1 <- fitmeasures(fit1, c("baseline.chisq", "baseline.df"))
user1 <- fitmeasures(fit1, c("chisq", "df"))

base2 <- fitmeasures(fit2, c("baseline.chisq", "baseline.df"))
user2 <- fitmeasures(fit2, c("chisq", "df"))

cfi(base1[1], base1[2], user1[1], user1[2])  # equals 1 b/c d.user==0
cfi(base2[1], base2[2], user2[1], user2[2])

