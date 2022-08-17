library(mvtnorm)

# this function is used to evaluate log of mean of Likelihood using LogLikelihood
# https://mc-stan.org/docs/2_23/stan-users-guide/log-sum-of-exponentials.html
log.mean.exp.fx = function(log.x) {
  N = length(log.x)
  M = max(log.x)
  -log(N) + M + log(sum(exp(log.x - M)))
}

# log.mean.exp.fx(log(1:1000))
# log(mean(1:1000))



# #-- generate sample --
# set.seed(555)
# 
# mu = sigma = 5
# alpha = (mu/sigma)^2
# beta  = mu/sigma^2
# 
# y = rgamma(100, alpha, beta) # 100: sample size
# mean(y); sd(y) # 4.245797 4.114383
# 
# y.mc = lapply(
#   y, function(x) { 
#     x.mc = runif(1000, max(0, x - 0.1), x + 0.1) # 1000 monte carlo samples
#     return(x.mc)
#   }
# )
# 
# head(y.mc)

### likelihood of a bivariate normal dist


#-- estimation --

# log-normal
llh.fx.ln = function(pars.est, pars.fix, data.mc) {
  ### input
  ## pars to estimate: meanlog and sdlog of GT; corr coef rho
  ## pars fix: meanlog and sdlog of IPinfector
  ## each element in data.mc is a matrix that x[,1] is logGT x[,2] is logIP
 # mu1 = pars.est[1]
#  sd1 = pars.est[2]
#  rho = pars.est[3]
  
  rho.raw = pars.est[3]
  # to make sure rho is within(-1, 1)
  rho = tanh(rho.raw)
  
  # pars.est[1] is the mean of  GT
  # pars.est[2] is the sd of GT
  
  mean1 <- pars.est[1]
  sigma1 <- pars.est[2]
  
  mu1 = log(mean1^2/sqrt(mean1^2 + sigma1^2))
  sd1 = sqrt(log(sigma1^2/mean1^2 + 1))
  
  mu2 = pars.fix[1]
  sd2 = pars.fix[2]
  covmat = matrix(c(sd1^2, rho * sd1 * sd2, rho * sd1 * sd2, sd2^2), nrow = 2, byrow = T)
  
  log.mean.likelihood = sapply(
    data.mc, function(x) { log.mean.exp.fx(
      dmvnorm(x, c(mu1, mu2), covmat, log=TRUE))
    }
  )
#  print(c(mean1, sigma1, rho))
#  print(-sum(log.mean.likelihood))
  return( -sum(log.mean.likelihood))
  
}


llh.fx.ln.adj = function(pars.est, pars.fix, data.mc) {
  ### input
  ## pars to estimate: meanlog and sdlog of GT; corr coef rho
  ## pars fix: meanlog and sdlog of IPinfector
  ## each element in data.mc is a matrix that x[,1] is logGT x[,2] is logIP
  # mu1 = pars.est[1]
  #  sd1 = pars.est[2]
  #  rho = pars.est[3]
  
  rho = pars.fix[1]
  # to make sure rho is within(-1, 1)
  
  # pars.est[1] is the mean of  GT
  # pars.est[2] is the sd of GT
  
  mean1 <- pars.est[1]
  sigma1 <- pars.est[2]
  
  mu1 = log(mean1^2/sqrt(mean1^2 + sigma1^2))
  sd1 = sqrt(log(sigma1^2/mean1^2 + 1))
  
  mu2 = pars.fix[2]
  sd2 = pars.fix[3]
  covmat = matrix(c(sd1^2, rho * sd1 * sd2, rho * sd1 * sd2, sd2^2), nrow = 2, byrow = T)
  
  log.mean.likelihood = sapply(
    data.mc, function(x) { log.mean.exp.fx(
      dmvnorm(x, c(mu1, mu2), covmat, log=TRUE))
    }
  )
  #  print(c(mean1, sigma1, rho))
  #  print(-sum(log.mean.likelihood))
  return( -sum(log.mean.likelihood))
  
}
















# start.t = Sys.time()
# 
# test = optim(
#   par = c(0,0),
#   fn  = llh.fx.ln,
#   method = "Nelder-Mead",
#   hessian = T,
#   y.mc = y.mc
# )
# end.t = Sys.time()
# 
# end.t - start.t # 2.790537 secs in windows



#-- result --
# exp(test$par) # 4.244340 4.194196
# qlnorm(0.025, test$par, sqrt(diag(solve(test$hessian))))
# qlnorm(0.975, test$par, sqrt(diag(solve(test$hessian))))
# 
# 
# 
# #-- bootstrap --
# library(doSNOW)
# library(parallel)
# n.cores = detectCores() - 1
# 
# start.t = Sys.time()
# 
# cl = makeCluster(n.cores, type = "SOCK")
# registerDoSNOW(cl)
# 
# TMP = foreach(
#   b = 1:1000
# ) %dopar% {
#   
#   set.seed(b)
#   
#   ind.b = sample(1:length(y.mc), length(y.mc), replace = T)
#   y.mc.b = y.mc[ind.b]
#   
#   test = optim(
#     par = c(0,0),
#     fn  = llh.fx,
#     method = "Nelder-Mead",
#     hessian = T,
#     y.mc = y.mc.b
#   )
#   return(test)
# }
# 
# stopCluster(cl)
# 
# end.t = Sys.time()
# end.t - start.t # 3.536265 mins in windows with 19 cores
# 
# 
# 
# b.est = sapply(TMP, function(b) { exp(b$par) })
# quantile(b.est[1,], c(0.025, 0.975))
# quantile(b.est[2,], c(0.025, 0.975))
# 

