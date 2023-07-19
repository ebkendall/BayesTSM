library(MASS)

# Simulates interval censored 3 state screening data assuming MV-normal covariates
# See demo.r for details
gendat_williams = function(n, p = 3, p.discrete = 0, r=.1, s = 1, 
                  sigma.X=1/2, mu.X = 1, beta.X = NULL, 
                  sigma.S=1/2, mu.S = 1, beta.S = NULL, 
                  cor.X.S = 0,
                  Tmax = 20, v.min = 1, v.max = 6, mean.rc = 40,
                  dist.X='weibull', dist.S='weibull',
                  do.b.XS = F, b.XS = NULL){
  
  if(p ==0){
    beta.X = as.matrix(mu.X)
    beta.S = as.matrix(mu.S)
  }   else{
  beta.X = as.matrix(c(mu.X, beta.X))
  beta.S = as.matrix(c(mu.S, beta.S))
  }
  
  # Sim Z
  R = matrix(r, p, p)
  diag(R) =  1
  S = rep(s, p)
  Sigma = cor2cov(R, S)
  if(p>0) {Z = mvrnorm(n, mu = rep(0,p), Sigma)} else Z = NULL
  # Sim discrete Z
  if(p.discrete == 1){
    Z.discrete = rbinom(n, 1, 0.5)
    Z = cbind( Z, Z.discrete)
    colnames(Z) = paste(1:ncol(Z))
  }  
  Z1 = cbind( as.matrix(rep(1,n)),Z)
  
  # Sim.X
  if( dist.X=='weibull' | dist.X=='loglog' | dist.X=='lognormal' ){
    if(dist.X=='weibull')    e.X = r.ev(n)
    if(dist.X=='loglog')     e.X = rlog(n)
    if(dist.X=='lognormal')  e.X = rnorm(n)
    X = Z1 %*% beta.X + sigma.X * e.X # log surv times
    X = exp(X) # surv times
    X  = as.numeric(X)
  }
  # Sim.S
  if( dist.S=='weibull' | dist.S=='loglog' | dist.S=='lognormal' ){
    if(dist.S=='weibull')   e.S = r.ev(n)
    if(dist.S=='loglog')    e.S = rlog(n)
    if(dist.S=='lognormal') e.S = rnorm(n)
    S = Z1 %*% beta.S + sigma.S * e.S # log surv times
    if(do.b.XS){
      #logstdX = scale(X)
      S = Z1 %*% beta.S + b.XS * log(X) + sigma.S * e.S # log surv times
    } 
    S = exp(S) # surv times
    S = as.numeric(S)
  }
  
  if(dist.X=='bv-lognormal'){
    Sigma.e  = matrix( c(1, cor.X.S, cor.X.S, 1), 2, 2)
    e = mvrnorm(n, c(0,0), Sigma.e)
    X = Z1 %*% rbind(mu.X, beta.X) + sigma.X * e[,1] # log surv times
    S = Z1 %*% rbind(mu.S, beta.S) + sigma.S * e[,2] # log surv times
    X = as.numeric(exp(X))
    S = as.numeric(exp(S))
  }
  
  # Create total time
  XS = as.numeric(X+S)  # Total time from baseline to event 2
  
  # Generate screening sequences
  t.rc = rexp(n, 1/mean.rc ) #runif(n, start.rc, vmax*(visitdist+leniancy))
  V = as.matrix(runif(n, v.min, v.max ))
  for(i in 2:Tmax){
    V = cbind(V, runif(n, V[,(i-1)]+v.min, V[,(i-1)]+v.max ))
  }
  t.rc = V[,1] + t.rc
  i = apply(V < t.rc, 1, sum)
  t.rc = apply( cbind(i, V), 1, function(x) x[x[1]+1] )
  t.rc = matrix( rep( t.rc, Tmax), nrow = n, ncol = Tmax)
  V    = (V < t.rc) * V + (V >= t.rc) * t.rc
  V    = cbind(0, V)
  
  # Find indices of events
  ind.h  = rowSums(V < X) # Last healthy time index
  ind.x  = ind.h+1        # First stage 1 
  ind.xs = rowSums(V < XS) + 1 # First stage 2
  ind.x[ind.x>ncol(V)] = ncol(V)
  
  # Define delta (event)
  d = rep(-99,n)
  d[ind.x < ind.xs ]  = 2
  d[ind.x == ind.xs  & ind.xs != (ncol(V)+1)] = 3
  d[ind.h == ind.x]  = 1
  
  # Define left and right interval bound
  L = V[ cbind(1:nrow(V), ind.h) ]
  R = V[ cbind(1:nrow(V), ind.x) ]
  R[ ind.h == ind.x ] = Inf
  
  if(p>0) dat = data.frame(L=L, R=R, d=d, X=X, S=S, Z=Z)
  if(p==0) dat = data.frame(L=L, R=R, d=d, X=X, S=S)
  
  
  time = state  = Z.1 = Z.2 =  numeric()
  copy = NULL
  i =1
  if(p>0){
    x = matrix(NA, ncol= 5+length(copy), nrow = 1)
  }else{
    x = matrix(NA, ncol= 3+length(copy), nrow = 1)
  }
  
  for(i in 1:nrow(dat)){
    time = numeric()
    state = numeric()
    Z.1 = numeric()
    Z.2 =  numeric()
    dd = dat$d
    if(dd[i] == 1){
      ns = length(unique(V[ i, 1: ind.x[i] ]))
      time[1:ns] = unique(V[ i, 1: ind.x[i] ])
      state[1:ns] = rep(1,ns)
      if(p>0){
        Z.1[1:ns]  = rep(Z[i,1],ns)
        Z.2[1:ns]  = rep(Z[i,2],ns)
      }
      
      
    }
    if(dd[i] == 2){
      ns = length(V[ i, 1: ind.x[i] ])
      time[1:ns] = unique(V[ i, 1: ind.x[i] ])
      state[1:ns] = rep(1,ns)
      state[ns] = 2
      if(p>0){
        Z.1[1:ns]  = rep(Z[i,1],ns)
        Z.2[1:ns]  = rep(Z[i,2],ns)
      }
    }
    if(dd[i] == 3){
      ns = length(V[ i, 1: ind.x[i] ])
      time[1:ns] = unique(V[ i, 1: ind.x[i] ])
      state[1:ns] = rep(1,ns)
      state[ns] = 3
      if(p>0){
        Z.1[1:ns]  = rep(Z[i,1],ns)
        Z.2[1:ns]  = rep(Z[i,2],ns)
      }
    }
    if(p>0){
      x.t = cbind(i, time, state, Z.1, Z.2)
    }else{
      x.t = cbind(i, time, state)
    }
    x.t = cbind(x.t, matrix( rep( as.matrix(dat[i,copy]), dim(x.t)[1]), nrow = dim(x.t)[1], ncol =length(copy), T))
    x = rbind(x, x.t)
  }
  x = x[-1,]
  
  
  if(p>0){
    data= data.frame( id = x[,1] ,time = x[,2], state = x[,3], Z.1 = x[,4], Z.2 =  x[,5])
  }else{
    data = data.frame( id = x[,1] ,time = x[,2], state = x[,3])
  }
  
  data
}

# Transforms correlation matrix to variance-covariance matrix using SDs
cor2cov <- function(R,S){
  diag(S) %*% R %*% diag(S)
}

# Quantile and random number generation for the extreme value distribution
q.ev = function(p) log(-log(1-p))
r.ev = function(n){
  u = runif(n)
  q.ev(u)
}

# Calculates marginal predictive CIF from model mod fit by bts_survreg

unlist.par = function(par.list, Z, dist, pst.samples, s) {
  par = as.matrix(par.list[1])
  p = ncol(par)
  beta = par[,1:p-1, drop=F]
  sigma = par[,p]
  for(i in 1:length(par.list)){
    par = as.matrix(par.list[i])
    beta = rbind(beta,par[,1:(p-1), drop=F])
    sigma = c(sigma, par[,p])
  }
  linterm = cbind(1,Z) %*% t(beta[s,])
  if(dist != 'lognormal'){
    a = sigma[s]^-1
    b = exp(linterm)
  }
  if(dist == 'lognormal'){
    a = sigma[s] 
    b = linterm 
  }
  x = apply(rbind(a,b),2, function(x) rdist(n = length(x)-1, par = cbind(x[2:length(x)],x[1]), dist = dist ) )
  as.matrix(x)
}

get.ppd = function(mod, pst.samples=10^3, perc = seq(0, 1, 0.01)) {
  par.list.X = (mod$par.X.bi)
  par.list.S = (mod$par.S.bi)
  if( !is.null(mod$Z.X) ) Z.X = as.matrix(mod$Z.X) else Z.X = mod$Z.X # mod$Z.X[sample(1:nrow(mod$Z.X),nrow(mod$Z.X), replace=T),] 
  if( !is.null(mod$Z.X) ) Z.S = as.matrix(mod$Z.S) else Z.S = mod$Z.S # Z.X
  dist.X = mod$dist.X
  dist.S = mod$dist.S
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  ret = list()
  x = unlist.par(par.list=par.list.X, Z.X, dist.X, pst.samples, s)
  s = unlist.par(par.list.S, Z.S, dist.S, pst.samples, s)
  q.x  = apply( x, 2, quantile, perc )
  q.s  = apply( s, 2, quantile, perc )
  q.xs = apply( x+s, 2, quantile, perc )
  ret$med.cdf.x     = apply(q.x, 1, median)
  ret$med.cdf.s     = apply(q.s, 1, median)
  ret$med.cdf.xs    = apply(q.xs, 1, median)
  ret$med.cdf.x.ci  = apply(q.x, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.s.ci  = apply(q.s, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.xs.ci = apply(q.xs, 1, quantile, c(0.025, 0.975))
  ret$perc = perc
  # ret$x.sample = as.numeric(ret$x)[sample(1:length(as.numeric(ret$x)), pst.samples, replace=F)]
  # ret$s.sample = as.numeric(ret$s)[sample(1:length(as.numeric(ret$s)), pst.samples, replace=F)]
  ret
}

get.ppd.grid = function(mod, pst.samples=10^3, q = seq(0, 20, 0.01)) {
  par.list.X = (mod$par.X.bi)
  par.list.S = (mod$par.S.bi)
  if( !is.null(mod$Z.X) ) Z.X = as.matrix(mod$Z.X) else Z.X = mod$Z.X # mod$Z.X[sample(1:nrow(mod$Z.X),nrow(mod$Z.X), replace=T),] 
  if( !is.null(mod$Z.X) ) Z.S = as.matrix(mod$Z.S) else Z.S = mod$Z.S # Z.X
  dist.X = mod$dist.X
  dist.S = mod$dist.S
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  ret = list()
  x = unlist.par(par.list=par.list.X, Z.X, dist.X, pst.samples, s)
  s = unlist.par(par.list.S, Z.S, dist.S, pst.samples, s)
  p.x  = apply(x, 2, function(y){ ecdf.x = ecdf(y); ecdf.x(q) })
  p.s  = apply(s, 2, function(y){ ecdf.s = ecdf(y); ecdf.s(q) })
  p.xs = apply(x + s, 2, function(y){ ecdf.s = ecdf(y); ecdf.s(q) })
  ret$med.cdf.x     = apply(p.x, 1, median)
  ret$med.cdf.s     = apply(p.s, 1, median)
  ret$med.cdf.xs    = apply(p.xs, 1, median)
  ret$med.cdf.x.ci  = apply(p.x, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.s.ci  = apply(p.s, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.xs.ci = apply(p.xs, 1, quantile, c(0.025, 0.975))
  ret$q = q
  ret
}

get.pCIF.q = get.ppd
get.pCIF.p = get.ppd.grid


# Function for use with bayes.2S only
get.ppd.2S = function(mod, pst.samples=10^3) {
  par.list.X = (mod$par.X.bi)
  Z.X = as.matrix(mod$Z.X)
  dist.X = mod$dist.X
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  ppd = unlist.par(par.list.X, Z.X, dist.X, pst.samples, s)
  ret = list()
  perc = seq(0, 1, 0.01)
  q = apply( ppd, 2, quantile, perc )
  ret$med.cdf = apply(q, 1, median)
  ret$med.cdf.ci = apply(q, 1, quantile, c(0.025, 0.975))
  ret$perc = perc
  ret$x = ppd
  ret$x.sample = as.numeric(ppd)[sample(1:length(as.numeric(ppd)), pst.samples, replace=F)]
  ret
}

# Running the data generating process
for(ttt in 1:10) {
  set.seed(ttt)
  print(ttt)
  bayestsm_dat = gendat_williams(n = 1000,     # Sample size
              p = 1,                 # Number of normally distributed covariates N(0,1)
              p.discrete = 1,        # Number of discrete covariates (Bernoulli(0.5))
              r = .1,                 # Correlation of the continuous covariates
              sigma.X = 0.2,         # True scale parameter
              mu.X    = 3,           # True intercept parameter
              beta.X  = c(0.5,0.5),  # True slope parameters of all covariates of X
              sigma.S = 0.3,         # True scale parameter
              mu.S    = 1.2,           # True intercept parameters
              beta.S  = c(0.5,0.5),  # True slope parameter of all covariates of S
              dist.X  = 'weibull',   # Distribution of X, alternatives are loglog and lognormal
              dist.S  = 'weibull',   # Distribution of S, alternatives are loglog and lognormal
              v.min   = 1,           # Minimum time between screening moments (c_min in the paper)
              v.max   = 7,           # Maximum time between screening moments (c_max in the paper)
              Tmax    = 2e2,         # Maximum number of screening times (this is set to high value; but too high values increase computation time)
              mean.rc = 20           # Mean time to right censoring (parameter of exponential distribution; theta in the paper)
  )
  save(bayestsm_dat, file = paste0("Data/bayestsm_dat", ttt, ".rda"))
}

# t1 = Sys.time()
# mod         = bts_survreg(d              = dat.bayestsm$d, # censoring indicator, assumes values 1, 2, 3
#                           L              = dat.bayestsm$L, # Time of left censoring; assumes time of last visit of d=1
#                           R              = dat.bayestsm$R, # Time of right censoring; assumes inf if d=1
#                           Z.X            = dat.bayestsm[,c('Z.1','Z.2')], # Covariates of X
#                           Z.S            = dat.bayestsm[,c('Z.1','Z.2')], # Covariates of S
#                           mc             = 5e4,       # MCMC draws (can be updated later, see below), half will be dropped for burn-in (other values can also be specified using brunin argument)
#                           chains         = 3,         # Number of parallel MCMC chains
#                           do.seperate.MH = F,         # Whether the metropolis step should be done jointly (F) or seperately for parameters of X and S
#                           prop.sd.X      = 0.006,      # The proposal standard deviation of the normal distribution used for the metropolis step (if do.seperate.MH = T also prop.sd.S should be tuned)
#                           beta.prior.X   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (X)
#                           beta.prior.S   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (S)
#                           sig.prior.X    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (X)
#                           sig.prior.S    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (S)
#                           dist.X         = 'weibull', # Distribution of X
#                           dist.S         = 'weibull', # Distribution of S
# )
# t2 = Sys.time()
# t2-t1

# # We obtain estimated CDF 
# cif.q = get.pCIF.q(mod = mod,             # The fitted model of bts_survreg
#                    pst.samples = 5e3,     # The number of posterior draws used in calculating the predictive CIF (computation time and accuracy scale with this number)
#                    perc = seq(0, 1, 0.01) # Vector of percentiles / probabilities at which we want to evaluate the inverse of the CDF (length scales computation time)
# )
# names(cif.q) # Returned elements are the medians and 95% credible intervals

# # Plots of obtained quantiles against specified percentiles
# par(mfrow=c(1,3))
# plot(cif.q$med.cdf.x, cif.q$perc, ty = 'l', main = 'Time X', xlim=c(0,150), lwd=2)
# lines(ecdf(dat.bayestsm$X),col=2, lwd=2) # True empirical CDF (data)
# lines( cif.q$med.cdf.x.ci[1,], cif.q$perc  ,col=3, lty=2, lwd=2)
# lines( cif.q$med.cdf.x.ci[2,], cif.q$perc  ,col=3, lty=2, lwd=2)

# plot(cif.q$med.cdf.s, cif.q$perc, ty = 'l', main = 'Time S', xlim=c(0,20), lwd=2)
# lines(ecdf(dat.bayestsm$S),col=2, lwd=2) # True empirical CDF (data)
# lines( cif.q$med.cdf.s.ci[1,], cif.q$perc  ,col=3, lty=2, lwd=2)
# lines( cif.q$med.cdf.s.ci[2,], cif.q$perc  ,col=3, lty=2, lwd=2)

# plot(cif.q$med.cdf.xs, cif.q$perc, ty = 'l', main = 'Time X+S', xlim=c(0,200), lwd=2)
# lines(ecdf(dat.bayestsm$X+dat.bayestsm$S),col=2, lwd=2) # True empirical CDF (data)
# lines( cif.q$med.cdf.xs.ci[1,], cif.q$perc  ,col=3, lty=2, lwd=2)
# lines( cif.q$med.cdf.xs.ci[2,], cif.q$perc  ,col=3, lty=2, lwd=2)