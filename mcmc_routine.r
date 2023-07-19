library(mvtnorm, quietly=T);
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(deSolve)
library(foreach)
library(doParallel, quietly=T)
sourceCpp("likelihood.cpp")

# Setting up the differential equations for deSolve
hmm_solver_P <-function(pars, par_index, t_i, x_i, k) {
    
    model_t <- function(t,p,parms) {
        
        betaMat <- matrix(parms$b, ncol = 4, byrow = T)
        
        q1  = exp( c(1,parms$x_ik,t) %*% betaMat[1,] ) 
        q2  = exp( c(1,parms$x_ik,t) %*% betaMat[2,] ) 
        
        # Vectorizing the matrix multiplication row-wise
        dP = rep(1,5)
        dP[1] = -p[1]*q1
        dP[2] = p[1]*q1 - p[2]*q2
        dP[3] = p[2]*q2
        
        dP[4] = -p[4]*q2
        dP[5] = p[4]*q2
        
        return(list(dP))
        
    }
    
    # NO par_index should be used
    beta <- pars[par_index$beta]
    
    p_ic <- c(p1=1,p2=0,p3=0,p4=1,p5=0)
    out <- deSolve::ode(p_ic, times = t_i[(k-1):k], func = model_t,
                        parms = list(b=beta, x_ik = x_i[k,]))

    P <- matrix(c(out[2,"p1"], out[2,"p2"], out[2,"p3"],
                            0, out[2,"p4"], out[2,"p5"], 
                            0,           0, 1), nrow = 3, byrow = T)
    return(P)
}

# Likelihood using HMM solver approach
fn_log_post_R <- function(pars, prior_par, par_index, x, y, t, id, eids) {
    init = matrix(c(1, 0, 0), ncol = 3)
    resp_fnc = diag(3);
    
    beta = pars[par_index$beta]
    
    log_total_val = foreach(i=unique(id), .combine='+', 
                            .export = c("hmm_solver_P"), 
                            .packages = c("deSolve", "expm")) %dopar% { 
                                
        y_i = y[id == i]                
        x_i = x[id == i, ,drop = F] 
        t_i = t[id == i]
                                
        D = diag(resp_fnc[,y_i[1]])
        f_i = init %*% D
                                
        for(k in 2:length(t_i)) { 
            
            # Time In-homogeneous Transition Matrix --------------------
            P = hmm_solver_P(pars, par_index, t_i, x_i, k)
            D = diag(resp_fnc[,y_i[k]])
            f_i = f_i %*% P %*% D
        }
                                
        return(log(sum(f_i)))
    }
    
    
    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=pars, mean=mean, sigma=sd, log=T)
    
    return(log_total_val + log_prior_dens)
}

# -----------------------------------------------------------------------------
# The mcmc routine for sampling the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y, x, t, id, init_par, prior_par, par_index, steps, burnin, disc, disc_t) {

  if(!disc) {
      cl <- makeCluster(15, outfile="")
      registerDoParallel(cl)      
  }
  
  pars = init_par
  n = length(y)
  n_par = length(pars)
  chain = matrix( 0, steps, n_par)

  # group = list(c(par_index$beta, par_index$misclass, par_index$pi_logit))
  group = list(c(par_index$beta))
  n_group = length(group)

  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
  pscale = rep( .0001, n_group)

  accept = rep( 0, n_group)
  eids = unique(id)

  # Evaluate the log posterior of the initial parameters
  if(disc) {
      log_post_prev = fn_log_post( pars, prior_par, par_index, x, y, t, id, eids, disc_t)
  } else {
      log_post_prev = fn_log_post_R( pars, prior_par, par_index, x, y, t, id, eids)
  }

  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }

  # Begin the MCMC algorithm --------------------------------------------------
  chain[1,] = pars
  for(ttt in 2:steps){
    for(j in 1:n_group){

      # Propose an update
      ind_j = group[[j]]
      proposal = pars
      proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])

      # Compute the log density for the proposal
      if(disc) {
          log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id, eids, disc_t)   
      } else {
          log_post = fn_log_post_R(proposal, prior_par, par_index, x, y, t, id, eids)   
      }

      # Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = pars
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                     sigma=pcov[[j]]*pscale[j])
          if(disc) {
              log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id, eids, disc_t)   
          } else {
              log_post = fn_log_post_R(proposal, prior_par, par_index, x, y, t, id, eids)   
          }
          
        }
      }

      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        pars[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
      }
      chain[ttt,ind_j] = pars[ind_j]

      # Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        # During the burnin period, update the proposal covariance in each step
        # to capture the relationships within the parameters vectors for each
        # transition.  This helps with mixing.
        if(ttt == 100)  pscale[j] = 1

        if(100 <= ttt & ttt <= 2000){
          temp_chain = chain[1:ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])

        } else if(2000 < ttt){
          temp_chain = chain[(ttt-2000):ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
        }
        if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )

        # Tune the proposal covariance for each transition to achieve
        # reasonable acceptance ratios.
        if(ttt %% 30 == 0){
          if(ttt %% 480 == 0){
            accept[j] = 0

          } else if( accept[j] / (ttt %% 480) < .4 ){ 
            pscale[j] = (.75^2)*pscale[j]

          } else if( accept[j] / (ttt %% 480) > .5 ){ 
            pscale[j] = (1.25^2)*pscale[j]
          }
        }
      }
      # -----------------------------------------------------------------------
    }
    # Restart the acceptance ratio at burnin.
    if(ttt == burnin)  accept = rep( 0, n_group)

    if(ttt%%1==0)  cat('--->',ttt,'\n')
  }
  # ---------------------------------------------------------------------------
  if(!disc) {
      stopCluster(cl)
  }
  print(accept/(steps-burnin))
  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale, pcov = pcov))
}
# -----------------------------------------------------------------------------