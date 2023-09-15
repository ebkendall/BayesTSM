# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(expm)

trial_num = 1

simulation = T
no_s2_s3 = F
disc = F

# Size of posterior sample from mcmc chains
n_post = 15000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps =  20000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = c(1:5) # the number of times the MCMC is run

par_index = list( beta=1:8)

true_vals = c(  3,  0.5, 0.5, 0.5, # beta(0, x), beta(1, x), beta(2, x), beta_time(x)
              1.2,  0.5, 0.5, 0.5) # beta(0, t), beta(1, t), beta(2, t), beta_time(t)

labels <- c('baseline(x)', 'betaZ1(x)', 'betaZ2(x)', 'beta_time(x)',
            'baseline(t)', 'betaZ1(t)', 'betaZ2(t)', 'beta_time(t)')

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))
chain_list <- NULL
ind = 0
for(seed in index_seeds){
    file_name = NULL
    if(simulation) {
        if(no_s2_s3) {
            if(disc) {
                print("A")
                file_name = paste0("Model_out/mcmc_out_extended_", 
                                   toString(seed), "_", trial_num, "_no_inhomog.rda")
            } else {
                print("B")
                file_name = paste0("Model_out/mcmc_out_", 
                                   toString(seed), "_", trial_num, "_no_inhomog.rda")
            }
        } else {
            if(disc) {
                print("C")
                file_name = paste0("Model_out/mcmc_out_extended_", 
                                   toString(seed), "_", trial_num, "_yes_inhomog.rda")
            } else {
                print("D")
                file_name = paste0("Model_out/mcmc_out_", 
                                   toString(seed), "_", trial_num, "_yes_inhomog.rda")
            }
        }
    } else {
        if(disc) {
            print("E")
            file_name = paste0("Model_out/mcmc_out_extended_", 
                               toString(seed), "_", trial_num, "_btsm_inhomog.rda")
        } else {
            print("F")
            file_name = paste0("Model_out/mcmc_out_", 
                               toString(seed), "_", trial_num, "_btsm_inhomog.rda")
        }
    }
    
    if(file.exists(file_name)){
        load(file_name) 
        ind = ind + 1
        print(mcmc_out$accept)
        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)
        # ind_keep = seq(1, nrow(main_chain), by=1)
        
        chain_list[[ind]] = main_chain[ind_keep, ]
        post_means[ind,] <- colMeans(main_chain[ind_keep, ])
    }
}

# Plot and save the mcmc trace plots and histograms.
stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))

plot_title = NULL
if(simulation) {
    if(no_s2_s3) {
        if(disc) {
            plot_title = paste0('Plots/mcmc_out_extended_', trial_num, '_no_inhomog.pdf')
        } else {
            plot_title = paste0('Plots/mcmc_out_', trial_num, '_no_inhomog.pdf')
        }
    } else {
        if(disc) {
            plot_title = paste0('Plots/mcmc_out_extended_', trial_num, '_yes_inhomog.pdf')
        } else {
            plot_title = paste0('Plots/mcmc_out_', trial_num, '_yes_inhomog.pdf')
        }
    }
} else {
    if(disc) {
        plot_title = paste0('Plots/mcmc_out_extended_', trial_num, '_btsm_inhomog.pdf')
    } else {
        plot_title = paste0('Plots/mcmc_out_', trial_num, '_btsm_inhomog.pdf')   
    }
}

pdf(plot_title)
par(mfrow=c(4, 2))
for(r in 1:length(labels)){
    
    xlabel = ""
    if(simulation) {xlabel = paste0('True Value: ', round(true_vals[r], 4))}
    
    plot( NULL, xlab=xlabel, ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
          ylim=range(stacked_chains[,r]) )
    
    for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)
    
    par_mean[r] = round( mean(stacked_chains[,r]), 4)
    par_median[r] = round( median(stacked_chains[,r]), 4)
    upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
    lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)
    
    print(paste(labels[r], ": [", lower[r], ", ", upper[r], "]"))
    
    hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
          freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
                              ' Median = ',toString(par_median[r])))
    abline( v=upper[r], col='red', lwd=2, lty=2)
    abline( v=lower[r], col='purple', lwd=2, lty=2)
    
    if(simulation) abline( v=true_vals[r], col='green', lwd=2, lty=2)
}

dev.off()

# ------------------------------------------------------------------------------
# Passage Probabilities --------------------------------------------------------
# ------------------------------------------------------------------------------
library(deSolve)

Q <- function(x_ik,beta){

    betaMat = matrix(beta, ncol = 4, byrow = T)
    q_x  = exp(c(1,x_ik) %*% betaMat[1,])  # Transition from state 1 to state 2
    q_t  = exp(c(1,x_ik) %*% betaMat[2,])  # Transition from state 2 to state 3

    qmat = matrix(c( 0,q_x,  0,
                     0,  0,q_t,
                     0,  0,  0), nrow=3, byrow = T)
    diag(qmat) = -rowSums(qmat)

    return(qmat)
}

pdf_title = NULL

if(simulation){
    if(no_s2_s3) {
        if(disc) {
            pdf_title = paste0('Plots/probEvol_extended_',trial_num,'_inhomog_no.pdf')
        } else {
            pdf_title = paste0('Plots/probEvol_',trial_num,'_inhomog_no.pdf')
        }
    } else {
        if(disc) {
            pdf_title = paste0('Plots/probEvol_extended_',trial_num,'_inhomog_yes.pdf')
        } else {
            pdf_title = paste0('Plots/probEvol_',trial_num,'_inhomog_yes.pdf')
        }
    }
} else {
    if(disc) {
        pdf_title = paste0('Plots/probEvol_extended_',trial_num,'_inhomog_btsm.pdf')
    } else {
        pdf_title = paste0('Plots/probEvol_',trial_num,'_inhomog_btsm.pdf')
    }
}

pdf(pdf_title)
par(mfrow=c(2, 2)) 

hmm_solver_P <-function(pars, t_i, x_i) {
    
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
    
    beta <- pars
    
    p_ic <- c(p1=1,p2=0,p3=0,p4=1,p5=0)
    out <- deSolve::ode(p_ic, times = t_i, func = model_t,
                        parms = list(b=beta, x_ik = x_i))
    
    P <- matrix(c(out[2,"p1"], out[2,"p2"], out[2,"p3"],
                  0, out[2,"p4"], out[2,"p5"], 
                  0,           0, 1), nrow = 3, byrow = T)
    return(P)
}

prob_evolution_hmmsolver <- function(state_from, state_to, z1, z2, par){
    
    t_seq = NULL
    if(state_to == 2) {
        t_seq = seq(0, 100, length.out=100)
    } else {
        t_seq = seq(0, 20, length.out=100)
    }

    probEvo = rep(NA,length(t_seq))
    
    P = rep(0,3)
    P[state_from] = 1
    index = 1
    
    for(t in 2:length(t_seq)) {
        beta = par
        new_P = hmm_solver_P(par, c(0,t_seq[t]), cov_vals)
        P_t = P %*% new_P
        
        # probEvo[index] = sum(P[state_to]) / sum(P[c(state_from, state_to)])
        if(state_to == 2) {
            probEvo[index] = P_t[state_to] + P_t[state_to+1]   
        } else {
            probEvo[index] = P_t[state_to]
        }
        
        index = index + 1
    }
    
    # cum_prob_evo = matrix(rep(t_seq, 2), ncol = 2)
    # cum_prob_evo[,2] = probEvo
    
    return(probEvo)
}

t.x = seq(0, 100, length.out=100)
t.y =  seq(0, 20, length.out=100)

# ------------------------------------------------------------------------------
# Load the dataset we are investigating ----------------------------------------
# ------------------------------------------------------------------------------
if(simulation) {
    if(no_s2_s3) {
        if(disc) {
            print("A")
            load(paste0("Data/hmm_sim_extended", ind, "_no_inhomog.rda"))
            myData = rawData
        } else {
            print("B")
            load(paste0("Data/hmm_sim", ind, "_no_inhomog.rda"))
            myData = rawData
        }
    } else {
        if(disc) {
            print("C")
            load(paste0("Data/hmm_sim_extended", ind, "_yes_inhomog.rda"))
            myData = rawData
        } else {
            print("D")
            load(paste0("Data/hmm_sim", ind, "_yes_inhomog.rda"))
            myData = rawData
        }
    }
} else {
    if(disc) {
        print("E")
        load(paste0("Data/bayestsm_dat_extended", ind, ".rda"))
        myData = bayestsm_dat
    } else {
        print("F")
        load(paste0("Data/bayestsm_dat", ind, ".rda"))
        myData = bayestsm_dat
    }
}

# create an empty matrix to store CDFs for each unique covariate combination
unique_ids = unique(myData[,'id'])
cdf.x = cdf.y = matrix(nrow = length(unique_ids), ncol = 100) 

par = colMeans(stacked_chains)
for(i in 1:nrow(cdf.x)) {
    subDat = myData[myData[,'id'] == unique_ids[i], ]
    z_1 = subDat[1,'Z.1']
    z_2 = subDat[1,'Z.2']
    cdf.x[i,] = prob_evolution_hmmsolver(1, 2, z_1, z_2, par)
    cdf.y[i,] = prob_evolution_hmmsolver(2, 3, z_1, z_2, par)
}


###  State 1 to 2
plot(t.x,colMeans(cdf.x),type="l",col=1,lty=1,lwd=2,ylim=c(0,1),
     xlab="time since entry in state 1",ylab="CDF",
     main = "Transition from state 1 to 2")

if(!simulation){
    lines(ecdf(dat.bayestsm$X),col=2,lwd=2)   # True empirical CDF (data)
    legend("bottomright", c( "Model estimate",  "True"), lwd=2, col=1:2, bty="n")
}   

###  State 2 to 3
plot(t.y,colMeans(cdf.y),type="l",col=1,lty=1,lwd=2,ylim=c(0,1),
     xlab="time since entry in state 2",ylab="CDF",
     main = "Transition from state 2 to 3")

if(!simulation) {
    lines(ecdf(dat.bayestsm$S),col=2,lwd=2)  # True empirical CDF (data)
    legend("bottomright", c( "Model estimate",  "True"), lwd=2, col=1:2, bty="n")
}