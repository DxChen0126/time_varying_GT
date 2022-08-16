

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


#-- estimation --
# gamma
llh.fx = function(pars, y.mc) {
  
  mean = exp(pars[1])
  sd   = exp(pars[2])
  
  alpha = (mean/sd)^2
  beta  = mean/(sd^2)
  
  log.mean.likelihood = sapply(
    y.mc, function(x) { 
      log.mean.exp.fx( dgamma(x, alpha, beta, log = T ) )
    }
  )
  return( -sum(log.mean.likelihood) )
}

# llh.fx(log(c(5,5)), y.mc)

# weibull
llh.fx.wb = function(pars, y.mc) {
  
  wbshape = exp(pars[1])
  wbscale = exp(pars[2])
  
  log.mean.likelihood = sapply(
    y.mc, function(x) { 
      log.mean.exp.fx( dweibull(x, shape = wbshape, scale = wbscale, log = T ) )
    }
  )
  return( -sum(log.mean.likelihood) )
}

# log-normal
llh.fx.ln = function(pars, y.mc) {
  
  mean = exp(pars[1])
  sd   = exp(pars[2])
  
  mu = log(mean^2/sqrt(mean^2 + sd^2))
  sigma = sqrt(log(sd^2/mean^2 + 1))
  
  log.mean.likelihood = sapply(
    y.mc, function(x) { 
      log.mean.exp.fx( dlnorm(x, meanlog = mu, sdlog = sigma, log = T ) )
    }
  )
  return( -sum(log.mean.likelihood) )
}

# start.t = Sys.time()
# 
# test = optim(
#   par = c(0,0),
#   fn  = llh.fx,
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

