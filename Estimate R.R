######### estimate r from our data
library(plyr)


allpairs <- read.csv("C:/Users/Dongxuan Chen/SPH Dropbox/Dongxuan Chen/Transmissibility_Methods/Data/China/allpairs0908.csv")

colexpo <- grep(c("expo"), colnames(allpairs)) 
for(i in colexpo){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

colonset <- grep(c("Onset"), colnames(allpairs)) 
for(i in colonset){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

length(unique(allpairs$Infector.ID))
length(unique(allpairs$Infectee.ID))

myfreq <- count(allpairs$Infector.ID)
summary(myfreq$freq)
hist(myfreq$freq, breaks = 16)



ind.ISboth <- which(allpairs$Infectee.ID %in% allpairs$Infector.ID)

intersect(allpairs$Infector.ID, allpairs$Infectee.ID)

# View(allpairs)

# note in this file some infectors are repeated if they infected more than one case

ind.infector <- numeric(length(unique(allpairs$Infector.ID)))

for(i in 1:length(unique(allpairs$Infector.ID))){
  ind.infector[i] <- which(allpairs$Infector.ID == unique(allpairs$Infector.ID)[i])[1]
}

df.onset <- data.frame(
  dates = c(allpairs$Onset_Infector[ind.infector], allpairs$Onset_Infectee),
  group = c(rep("infector", length(unique(allpairs$Infector.ID))), rep("infectee", nrow(allpairs)))
)


df.onset.infector <- df.onset[df.onset$group == "infector",]


library(incidence)

View(df.onset.infector)

infector.incidence <- incidence(df.onset.infector$dates)
infector.incidence$counts

fit.r.infector <- fit(infector.incidence[1:20])

# r = 0.1499 (0.0825, 0.2172)

plot(infector.incidence[1:20], fit = fit.r.infector, border = "white", color = "grey40")

df.onset.alt <- data.frame(
  dates = c(allpairs$Onset_Infector[ind.infector], allpairs$Onset_Infectee[-ind.ISboth]),
  group = c(rep("infector", length(unique(allpairs$Infector.ID))), rep("infectee", nrow(allpairs[-ind.ISboth,])))
)

all.incidence <- incidence(df.onset.alt$dates)
all.incidence$date

fit.r.all <- fit(all.incidence[1:23])

plot(all.incidence, fit = fit.r.all, border = "white", color = "grey40")

# r = 0.1981 (0.1468, 0.2494)

#### CI too wide, refer to Ren et al's 0.14 (0.11, 0.17)

### assume r follows normal distribution with mu = 0.14, sigma = 0.03/1.96

############## from Ganyani's code to estimate R0, use bootstrap parameters of GT to estimate R0 and CI
load("forwardGT1000_bootpars_DEC30.RData")

# use the formula as R = 1/int_{0}{inf}exp(-ra)g(a)da, where r is growth rate g(a) is GT dist

# here we use our estimated GT based on lnorm dist

calculateR <- function(r, mlog, slog){
  integrand <- function(x) {dlnorm(x, meanlog = mlog, sdlog = slog)*exp(-r*x)}
  intg <- integrate(integrand, lower = 0, upper = Inf)
  1/intg$value
}

df.estGT.ln <- read.csv("GTfor_main_estlnorm_DEC30.csv")

mu <- df.estGT.ln$mu[1]
sd <- df.estGT.ln$sd[1]
# R.est.main <- calculateR(r = 0.14, mlog = log(mu^2/sqrt(mu^2+sd^2)),
#           slog = sqrt(log(sd^2/mu^2+1)))

# R.est.main # 2.474452 

# Tim's estimation of r outside Hubei before Jan 23: 0.10 (0.08, 0.12)

R.est.main <- calculateR(r = 0.10, mlog = log(mu^2/sqrt(mu^2+sd^2)),
                         slog = sqrt(log(sd^2/mu^2+1)))

R.est.main # 1.947185


bootmeanlog <- lapply(parsboot, function(x) log(x[1,]^2/sqrt(x[1,]^2+x[2,]^2)))
bootsdlog <- lapply(parsboot, function(x) sqrt(log(x[2,]^2/x[1,]^2+1)))

meanlog = bootmeanlog[[1]]
sdlog = bootsdlog[[1]]

set.seed(123)
# rsim = rnorm(1000, 0.14, 0.03/1.96)
# tmpR = matrix(ncol = 1000, nrow = 1000, byrow = T)

rsim = rnorm(1000, 0.10, 0.02/1.96)
tmpR = matrix(ncol = 1000, nrow = 1000, byrow = T)

progressbar <- txtProgressBar(min = 0, max = 1000, style = 3)
for(i in 1:1000){
  for(j in 1:1000){
    tmpR[i, j] <- calculateR(r = rsim[i], mlog = meanlog[j], slog = sdlog[j])
  }
  setTxtProgressBar(progressbar, i)
}


# quantile(as.vector(tmpR), 0.025) # 2.040997
# quantile(as.vector(tmpR), 0.975) # 3.046277

quantile(as.vector(tmpR), 0.025) # 1.698191
quantile(as.vector(tmpR), 0.975) # 2.256993


# save(tmpR, file = "simulateR0_DEC30.RData")
save(tmpR, file = "simulateR0_May04.RData")


# R0 by backward


load("backwardGT1000_bootpars_gamma_MAY06.rda")

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

bootshape <- lapply(parsboot, function(x) x[1,]^2/x[2,]^2)
bootrate <- lapply(parsboot, function(x) x[1,]/x[2,]^2)

bootshape1 = bootshape[[1]]
bootrate1 = bootrate[[1]]

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

save(tmpR, file = "simulateR0_bybackwardGT_May06.rda")







