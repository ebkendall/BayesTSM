library(mvtnorm, quietly=T);library(foreach, quietly=T);library(msm, quietly=T)
library(doParallel, quietly=T);library(deSolve, quietly=T);library(expm, quietly=T)

library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("likelihood.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# Construct the transition rate matrix
Q_r <- function(x_ik,beta){
    
    betaMat = matrix(beta, ncol = 3, byrow = T) 
    q_x  = exp(c(1,x_ik) %*% betaMat[1,])  # Transition from state 1 to state 2
    q_t  = exp(c(1,x_ik) %*% betaMat[2,])  # Transition from state 2 to state 3
    
    qmat = matrix(c( 0,q_x,  0,
                     0,  0,q_t,
                     0,  0,  0), nrow=3, byrow = T)
    diag(qmat) = -rowSums(qmat)
    
    return(qmat)
}

fn_log_post_r <- function(pars, prior_par, par_index, x, y, t, id) {
    
    # Initial state probabilities
    # init_logit = c( 1, exp(pars[par_index$pi_logit][1]), exp(pars[par_index$pi_logit][2]))
    # init = init_logit / sum(init_logit)
    init = c(1,0,0)
    
    # Misclassification response matrix
    # resp_fnc = matrix(c(1, exp(pars[par_index$misclass][1]), exp(pars[par_index$misclass][2]),
    #                     exp(pars[par_index$misclass][3]), 1, exp(pars[par_index$misclass][4]),
    #                     exp(pars[par_index$misclass][4]), exp(pars[par_index$misclass][5]), 1),
    #                     ncol=3, byrow=TRUE)
    # resp_fnc = resp_fnc / rowSums(resp_fnc)
    resp_fnc = diag(3)
    
    beta <- pars[par_index$beta]
    # p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=0,p9=1) # initial condition for deSolve
    
    # Parallelized computation of the log-likelihood
    log_total_val = foreach(i=unique(id), .combine='+', .export = "Q_r", 
                            .packages = "expm") %dopar% {
                                
                                val = 1
                                y_i = y[id == i]                # the observed state
                                x_i = x[id == i, ,drop = F]     # Z1 and Z2
                                t_i = t[id == i]                # continuous time
                                
                                f_i = init %*% diag(resp_fnc[, y_i[1]])
                                # log_norm = 0
                                
                                for(k in 2:length(t_i)) {
                                    
                                    P   <- expm((t_i[k] - t_i[k-1]) * Q_r(x_i[k-1,], beta))
                                    f_i <- f_i %*% P %*% diag(resp_fnc[, y_i[k]])
                                    
                                    # out <- deSolve::ode(p_ic, times = t_i[(k-1):k], func = model_t,
                                    #                     parms = list(b=beta, x_ik = x_i[k,]))
                                    # P <- matrix(c(out[2,"p1"], out[2,"p2"], out[2,"p3"], 
                                    #               out[2,"p4"], out[2,"p5"], out[2,"p6"],
                                    #               out[2,"p7"], out[2,"p8"], out[2,"p9"]), nrow = 3, byrow = T)
                                    
                                    # P <-  expm((t_i[k] - t_i[k-1]) * Q(x_i[k-1,], beta))
                                    # val = f_i %*% P %*% diag(resp_fnc[, y_i[k]])
                                    # norm_val = sqrt(sum(val^2))
                                    # f_i = val / norm_val
                                    # log_norm = log_norm + log(norm_val)
                                }
                                
                                # return(log(sum(f_i)) + log_norm)
                                return(log(sum(f_i)))
                            }
    
    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=pars, mean=mean, sigma=sd, log=T)
    
    return(log_total_val + log_prior_dens)
    
}


init_par = c(  3,  0.5, 0.5,
               1.2,  0.5, 0.5) 

par_index = list( beta=1:6) 

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load(paste0("Data/hmm_sim1.rda"))
temp_data = as.matrix(rawData); rownames(temp_data) = NULL

id = temp_data[, "id"]
y  = temp_data[, "state"]
x  = temp_data[, c("Z.1", "Z.2"), drop=F]
t  = temp_data[, "time"]


# Testing
x_ik = c(2,2)
beta = c(3,0.5, 0.5, 1.2,0.5,0.5)

test_r = Q_r(x_ik, beta)

test_c = Q(x_ik, beta)

print(test_r)
print(test_c)

n_cores = 4
cl <- makeCluster(n_cores, outfile="")
registerDoParallel(cl)

likelihood_r = fn_log_post_r(init_par, prior_par, par_index, x, y, t, id)

stopCluster(cl)

eids = unique(id)
likelihood_c = fn_log_post(init_par, prior_par, par_index, x, y, t, id, eids)

print(likelihood_r)
print(likelihood_c)
