# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(expm)
library(deSolve)

trial_num = 1

no_s2_s3 = F
disc = F

# *****************************************************************************
# THESE NUMBERS WILL CHANGE IF THE NUMBER OF MCMC STEPS CHANGES ***************
# Size of posterior sample from mcmc chains
n_post = 15000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps =  20000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)
# *****************************************************************************

# *****************************************************************************
# NUMBER OF SIMULATIONS RUN (CORRESPONDS TO THE RANDOM SEED IN MY CASE) *******
index_seeds = c(1:5) 
# *****************************************************************************

par_index = list( beta=1:8)

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
    if(disc) {
        file_name = paste0("Model_out/mcmc_out_extended_", 
                            toString(seed), "_", trial_num, "_btsm_inhomog.rda")
    } else {
        file_name = paste0("Model_out/mcmc_out_", 
                            toString(seed), "_", trial_num, "_btsm_inhomog.rda")
    }
    
    if(file.exists(file_name)){
        load(file_name) 
        ind = ind + 1
        print(mcmc_out$accept)

        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)
        
        chain_list[[ind]] = main_chain[ind_keep, ]
        post_means[ind,] <- colMeans(main_chain[ind_keep, ])
    }
}

# Plot and save the mcmc trace plots and histograms.
stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))

plot_title = NULL
if(disc) {
    plot_title = paste0('Plots/mcmc_out_', trial_num, '_extended_btsm.pdf')
} else {
    plot_title = paste0('Plots/mcmc_out_', trial_num, '_btsm.pdf')   
}

pdf(plot_title)
par(mfrow=c(4, 2))
for(r in 1:length(labels)){
    
    xlabel = ""
    
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
}

dev.off()