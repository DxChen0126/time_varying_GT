
load("empirical_corr_ref.rda")

# expo = 4
load("GT_est_corrsims_medGT_expo4.rda")

# rho = 0

estsd <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sd)), ncol = 6, 
                     byrow = T)[-c(22, 24),]

estsdLB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdLB)), ncol = 6, 
                       byrow = T)[-c(22, 24),]

estsdUB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdUB)), ncol = 6, 
                       byrow = T)[-c(22, 24),]

# -22, - 24   48 sims left

estmu <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-c(22, 24),]

estmuLB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-c(22, 24),]

estmuUB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-c(22, 24),]

estrho <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$rho)), ncol = 6, 
                byrow = T)[-c(22, 24),]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$rhoLB)), ncol = 6, 
                  byrow = T)[-c(22, 24),]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$rhoUB)), ncol = 6, 
                  byrow = T)[-c(22, 24),]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[1]]$meanGTforonset
ref.sd <- df.sum[[1]]$sdGTforonset
ref.rho <- df.sum[[1]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo4GTsim_rho0 <- data.frame(
  empmean = df.sum[[1]]$meanGTforonset,
  empsd = df.sum[[1]]$sdGTforonset,
  emprho = df.sum[[1]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[1]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


# rho = 0.25

estsd <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sd)), ncol = 6, 
                byrow = T)[-47,]

estsdLB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sdLB)), ncol = 6, 
                  byrow = T)[-47,]

estsdUB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sdUB)), ncol = 6, 
                  byrow = T)[-47,]

# -47 48 sims left

estmu <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-47,]

estmuLB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-47,]

estmuUB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-47,]

estrho <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$rho)), ncol = 6, 
                 byrow = T)[-47,]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$rhoLB)), ncol = 6, 
                   byrow = T)[-47,]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$rhoUB)), ncol = 6, 
                   byrow = T)[-47,]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[2]]$meanGTforonset
ref.sd <- df.sum[[2]]$sdGTforonset
ref.rho <- df.sum[[2]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo4GTsim_rho0.25 <- data.frame(
  empmean = df.sum[[2]]$meanGTforonset,
  empsd = df.sum[[2]]$sdGTforonset,
  emprho = df.sum[[2]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[2]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


# rho = 0.5

estsd <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sd)), ncol = 6, 
                byrow = T)[-45,]

estsdLB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdLB)), ncol = 6, 
                  byrow = T)[-45,]

estsdUB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdUB)), ncol = 6, 
                  byrow = T)[-45,]

# -45 49 sims left

estmu <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-45,]

estmuLB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-45,]

estmuUB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-45,]

estrho <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$rho)), ncol = 6, 
                 byrow = T)[-45,]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$rhoLB)), ncol = 6, 
                   byrow = T)[-45,]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$rhoUB)), ncol = 6, 
                   byrow = T)[-45,]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[3]]$meanGTforonset
ref.sd <- df.sum[[3]]$sdGTforonset
ref.rho <- df.sum[[3]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo4GTsim_rho0.5 <- data.frame(
  empmean = df.sum[[3]]$meanGTforonset,
  empsd = df.sum[[3]]$sdGTforonset,
  emprho = df.sum[[3]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[3]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


# rho = 0.75

estsd <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sd)), ncol = 6, 
                byrow = T)[-c(15, 34),]

estsdLB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdLB)), ncol = 6, 
                  byrow = T)[-c(15, 34),]

estsdUB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdUB)), ncol = 6, 
                  byrow = T)[-c(15, 34),]

# -15 -34 48 sims left

estmu <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-c(15, 34),]

estmuLB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-c(15, 34),]

estmuUB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-c(15, 34),]

estrho <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$rho)), ncol = 6, 
                 byrow = T)[-c(15, 34),]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$rhoLB)), ncol = 6, 
                   byrow = T)[-c(15, 34),]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$rhoUB)), ncol = 6, 
                   byrow = T)[-c(15, 34),]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[4]]$meanGTforonset
ref.sd <- df.sum[[4]]$sdGTforonset
ref.rho <- df.sum[[4]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo4GTsim_rho0.75 <- data.frame(
  empmean = df.sum[[4]]$meanGTforonset,
  empsd = df.sum[[4]]$sdGTforonset,
  emprho = df.sum[[4]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[4]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


df.expo4GTsim <- rbind(df.expo4GTsim_rho0, df.expo4GTsim_rho0.25,
                       df.expo4GTsim_rho0.5, df.expo4GTsim_rho0.75)

write.csv(df.expo4GTsim, "df.expo4GTsim_extremedeleted.csv")


######## mean expo = 7
load("GT_est_corrsims_medGT_expo7.rda")

# rho = 0

estsd <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sd)), ncol = 6, 
                byrow = T)[-c(18, 20, 27, 38, 41),]

estsdLB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdLB)), ncol = 6, 
                  byrow = T)[-c(18, 20, 27, 38, 41),]

estsdUB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdUB)), ncol = 6, 
                  byrow = T)[-c(18, 20, 27, 38, 41),]

# -18, - 20, - 27, - 38, - 41   36 sims left

estmu <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-c(18, 20, 27, 38, 41),]

estmuLB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-c(18, 20, 27, 38, 41),]

estmuUB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-c(18, 20, 27, 38, 41),]

estrho <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$rho)), ncol = 6, 
                 byrow = T)[-c(18, 20, 27, 38, 41),]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$rhoLB)), ncol = 6, 
                   byrow = T)[-c(18, 20, 27, 38, 41),]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$rhoUB)), ncol = 6, 
                   byrow = T)[-c(18, 20, 27, 38, 41),]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[1]]$meanGTforonset
ref.sd <- df.sum[[1]]$sdGTforonset
ref.rho <- df.sum[[1]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo7GTsim_rho0 <- data.frame(
  empmean = df.sum[[1]]$meanGTforonset,
  empsd = df.sum[[1]]$sdGTforonset,
  emprho = df.sum[[1]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[1]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


# rho = 0.25

estsd <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sd)), ncol = 6, 
                byrow = T)[-c(16, 27),]

estsdLB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sdLB)), ncol = 6, 
                  byrow = T)[-c(16, 27),]

estsdUB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sdUB)), ncol = 6, 
                  byrow = T)[-c(16, 27),]

# -16 -27  33sims left

estmu <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-c(16, 27),]

estmuLB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-c(16, 27),]

estmuUB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-c(16, 27),]

estrho <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$rho)), ncol = 6, 
                 byrow = T)[-c(16, 27),]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$rhoLB)), ncol = 6, 
                   byrow = T)[-c(16, 27),]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$rhoUB)), ncol = 6, 
                   byrow = T)[-c(16, 27),]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[2]]$meanGTforonset
ref.sd <- df.sum[[2]]$sdGTforonset
ref.rho <- df.sum[[2]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo7GTsim_rho0.25 <- data.frame(
  empmean = df.sum[[2]]$meanGTforonset,
  empsd = df.sum[[2]]$sdGTforonset,
  emprho = df.sum[[2]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[2]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


# rho = 0.5

estsd <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sd)), ncol = 6, 
                byrow = T)[-c(30, 35, 36),]

estsdLB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdLB)), ncol = 6, 
                  byrow = T)[-c(30, 35, 36),]

estsdUB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdUB)), ncol = 6, 
                  byrow = T)[-c(30, 35, 36),]

# -30, 35, 36 42 sims left

estmu <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-c(30, 35, 36),]

estmuLB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-c(30, 35, 36),]

estmuUB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-c(30, 35, 36),]

estrho <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$rho)), ncol = 6, 
                 byrow = T)[-c(30, 35, 36),]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$rhoLB)), ncol = 6, 
                   byrow = T)[-c(30, 35, 36),]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$rhoUB)), ncol = 6, 
                   byrow = T)[-c(30, 35, 36),]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[3]]$meanGTforonset
ref.sd <- df.sum[[3]]$sdGTforonset
ref.rho <- df.sum[[3]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo7GTsim_rho0.5 <- data.frame(
  empmean = df.sum[[3]]$meanGTforonset,
  empsd = df.sum[[3]]$sdGTforonset,
  emprho = df.sum[[3]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[3]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


# rho = 0.75

estsd <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sd)), ncol = 6, 
                byrow = T)[-c(9, 15, 31, 43),]

estsdLB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdLB)), ncol = 6, 
                  byrow = T)[-c(9, 15, 31, 43),]

estsdUB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdUB)), ncol = 6, 
                  byrow = T)[-c(9, 15, 31, 43),]

# -9 15 31 43  42 sims left

estmu <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$mu)), ncol = 6, 
                byrow = T)[-c(9, 15, 31, 43),]

estmuLB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muLB)), ncol = 6, 
                  byrow = T)[-c(9, 15, 31, 43),]

estmuUB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muUB)), ncol = 6, 
                  byrow = T)[-c(9, 15, 31, 43),]

estrho <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$rho)), ncol = 6, 
                 byrow = T)[-c(9, 15, 31, 43),]

estrhoLB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$rhoLB)), ncol = 6, 
                   byrow = T)[-c(9, 15, 31, 43),]

estrhoUB <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$rhoUB)), ncol = 6, 
                   byrow = T)[-c(9, 15, 31, 43),]


mu.mean <- colMeans(estmu)  

sd.mean <- colMeans(estsd)  

rho.mean <- colMeans(estrho)  

ref.mu <- df.sum[[4]]$meanGTforonset
ref.sd <- df.sum[[4]]$sdGTforonset
ref.rho <- df.sum[[4]]$corres.GTIPback

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
ind.rho <- Bias.rho <- Width.rho <- numeric(6)

for(j in 1:6){
  ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
  ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
  ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
  
  Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
  Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
  Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
  
  Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
  Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
  Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
}

df.expo7GTsim_rho0.75 <- data.frame(
  empmean = df.sum[[4]]$meanGTforonset,
  empsd = df.sum[[4]]$sdGTforonset,
  emprho = df.sum[[4]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[4]]$samplesize.GTIPback,
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd,
  rhoCIrecover = ind.rho,
  estrho = rho.mean,
  biasrho = Bias.rho,
  widthrho = Width.rho
)


df.expo7GTsim <- rbind(df.expo7GTsim_rho0, df.expo7GTsim_rho0.25,
                       df.expo7GTsim_rho0.5, df.expo7GTsim_rho0.75)

write.csv(df.expo7GTsim, "df.expo7GTsim_extremedeleted.csv")






