######### estimate r from our data
library(dplyr)


############## from Ganyani's code to estimate R0, use bootstrap parameters of GT to estimate R0 and CI
# load("forwardGT1000_bootpars_lnorm_MAY06.rda")

GTboootpars_firstwindow <- parsboot[[1]]

save(GTboootpars_firstwindow, file = "forwardGT_bootpars_firstwindow.rda")


# use the formula as R = 1/int_{0}{inf}exp(-ra)g(a)da, where r is growth rate g(a) is GT dist

# here we use our estimated GT based on lnorm dist

calculateR <- function(r, mlog, slog){
  integrand <- function(x) {dlnorm(x, meanlog = mlog, sdlog = slog)*exp(-r*x)}
  intg <- integrate(integrand, lower = 0, upper = Inf)
  1/intg$value
}

df.estGT.ln <- read.csv("forwardGT_lnorm_May06.csv")

mu <- df.estGT.ln$mu[1]
sd <- df.estGT.ln$sd[1]


# Tim's estimation of r outside Hubei before Jan 23: 0.10 (0.08, 0.12)

R.est.main <- calculateR(r = 0.10, mlog = log(mu^2/sqrt(mu^2+sd^2)),
                         slog = sqrt(log(sd^2/mu^2+1)))

R.est.main # 1.947185

parsboot <- GTboootpars_firstwindow

meanlog <- log(parsboot[1,]^2/sqrt(parsboot[1,]^2+parsboot[2,]^2))
sdlog <- sqrt(log(parsboot[2,]^2/parsboot[1,]^2+1))

set.seed(123)

rsim = rnorm(1000, 0.10, 0.02/1.96)
tmpR = matrix(ncol = 1000, nrow = 1000, byrow = T)

progressbar <- txtProgressBar(min = 0, max = 1000, style = 3)
for(i in 1:1000){
  for(j in 1:1000){
    tmpR[i, j] <- calculateR(r = rsim[i], mlog = meanlog[j], slog = sdlog[j])
  }
  setTxtProgressBar(progressbar, i)
}


quantile(as.vector(tmpR), 0.025) # 1.698191
quantile(as.vector(tmpR), 0.975) # 2.256993


# save(tmpR, file = "simulateR0_May04.RData")


# R0 by backward


# load("backwardGT1000_bootpars_gamma_MAY06.rda")

GTboootpars_firstwindow.back <- parsboot[[1]]

save(GTboootpars_firstwindow.back, file = "backwardGT_bootpars_firstwindow.rda")


df.GTback.gamma <- read.csv("backwardGT_est_gamma.csv")

calculateR <- function(r, shape.b, rate.b){
  integrand <- function(x) {dgamma(x,  shape = shape.b, rate = rate.b)*exp(-r*x)}
  intg <- integrate(integrand, lower = 0, upper = Inf)
  1/intg$value
}

mu <- df.GTback.gamma$mu[1]
sd <- df.GTback.gamma$sd[1]

R.est.main <- calculateR(r = 0.10, shape.b = mu^2/sd^2, rate.b = mu/sd^2)

R.est.main # 1.573701

parsboot <- GTboootpars_firstwindow.back 

bootshape1 <-  parsboot[1,]^2/parsboot[2,]^2
bootrate1 <- parsboot[1,]/parsboot[2,]^2


set.seed(123)
rsim = rnorm(1000, 0.10, 0.02/1.96)
tmpR = matrix(ncol = 1000, nrow = 1000, byrow = T)

progressbar <- txtProgressBar(min = 0, max = 1000, style = 3)
for(i in 1:1000){
  for(j in 1:1000){
    tmpR[i, j] <- calculateR(r = rsim[i], shape.b = bootshape1[j], rate.b = bootrate1[j])
  }
  setTxtProgressBar(progressbar, i)
}


quantile(as.vector(tmpR), 0.025) # 1.434052
quantile(as.vector(tmpR), 0.975) # 1.741322

# save(tmpR, file = "simulateR0_bybackwardGT_May06.rda")







