# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(expm)

args = commandArgs(TRUE)

simulation = args[1]
no_s2_s3 = args[2]
trial_num = 1
time_inhomog = T

# Size of posterior sample from mcmc chains
n_post = 10000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps =  20000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = c(1:10) # the number of times the MCMC is run

par_index = list( beta=1:6)
# par_index = list( beta=1:6, misclass=7:12, pi_logit=13:14)

true_vals = c(  3,  0.5, 0.5, 0.5, # beta(0, x), beta(1, x), beta(2, x)
                1.2,  0.5, 0.5, 0.5) # beta(0, t), beta(1, t), beta(2, t)
#  -6, -6,          # P(obs S1 | S2), P(obs S1 | S3)
#  -6, -6,          # P(obs S2 | S1), P(obs S2 | S3)
#  -6, -6,          # P(obs S3 | S1), P(obs S3 | S2)
#  -6, -6)          # P(init S2), P(init S3)

labels <- c('baseline(x)', 'betaZ1(x)', 'betaZ2(x)', 'beta_time(x)',
            'baseline(t)', 'betaZ1(t)', 'betaZ2(t)', 'beta_time(t)')
# 'logit P( obs. state 1 | true state 2 )',
# 'logit P( obs. state 1 | true state 3 )',
# 'logit P( obs. state 2 | true state 1 )',
# 'logit P( obs. state 2 | true state 3 )',
# 'logit P( obs. state 3 | true state 1 )',
# 'logit P( obs. state 3 | true state 2 )',
# 'logit P( init. state 2 )','logit P( init. state 3 )')

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
cred_set = rep(list(matrix(ncol=2,nrow=length(index_seeds))), length(true_vals))

post_means = matrix(nrow = length(index_seeds), ncol = length(labels))
chain_list <- NULL
ind = 0
for(seed in index_seeds){
    file_name = NULL
    if(simulation) {
        if(no_s2_s3) {
            file_name = paste0('Model_out/mcmc_out_',toString(seed), '_', trial_num, '_no_inhomog.rda')
        } else {
            file_name = paste0('Model_out/mcmc_out_',toString(seed), '_', trial_num, '_yes_inhomog.rda')
        }
    } else {
        file_name = paste0('Model_out/mcmc_out_',toString(seed), '_', trial_num, '_btsm_inhomog.rda')
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
        
        # Calculating Credible Sets --------------------------------------------
        for(j in 1:length(true_vals)) {
            
            cred_set[[j]][ind, 1] = round(quantile( mcmc_out$chain[ind_keep,j],
                                                    prob=.025), 4)
            cred_set[[j]][ind, 2] = round(quantile( mcmc_out$chain[ind_keep,j],
                                                    prob=.975), 4)   
        }
    }
}

# Calculating Coverage --------------------------------------------------------
cov_df = rep(NA, length(true_vals))
for(i in 1:length(true_vals)) {
    val = true_vals[i]
    cov_df[i] = mean(cred_set[[i]][,1] <= val & val <= cred_set[[i]][,2], na.rm=T)
}

# Plot and save the mcmc trace plots and histograms.
stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))

plot_title = NULL
if(simulation) {
    if(no_s2_s3) {
        plot_title = paste0('Plots/mcmc_out_', trial_num, '_no_inhomog.pdf')
    } else {
        plot_title = paste0('Plots/mcmc_out_', trial_num, '_yes_inhomog.pdf')
    }
} else {
    plot_title = paste0('Plots/mcmc_out_', trial_num, '_btsm_inhomog.pdf')
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
plot.new()
plot.new()
# cred_set_cumulative = cbind(lower, upper)
# save(cred_set_cumulative, file = paste0('real_cav_analysis/Plots/cred_set_cumulative_', 
#                                             model_name[folder], '.rda'))

VP <- vector(mode="list", length = length(labels))
for(r in 1:length(labels)) {
    # Adding the boxplots
    yVar = post_means[,r]
    x_label = y_label = ""
    
    if(simulation) {
        x_label = paste0("Coverage is: ", round(cov_df[r], 3))
        y_label = paste0("Parameter Value: ", round(true_vals[r], 3))
    }
    
    plot_df = data.frame(yVar = yVar, xVar = rep(labels[r], length(yVar)))
    VP[[r]] = ggplot(plot_df, aes(y = yVar, x = xVar)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1) +
        ggtitle(labels[r]) +
        ylab(y_label) +
        xlab(x_label) +
        theme(text = element_text(size = 7))
}

grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]], VP[[6]], ncol=3, nrow =2)

dev.off()

# ------------------------------------------------------------------------------
# Probability evolution curves -------------------------------------------------
# ------------------------------------------------------------------------------
Q <- function(x_ik,beta){
    
    betaMat = matrix(beta, ncol = 3, byrow = T) 
    q_x  = exp(c(1,x_ik) %*% betaMat[1,])  # Transition from state 1 to state 2
    q_t  = exp(c(1,x_ik) %*% betaMat[2,])  # Transition from state 2 to state 3
    
    qmat = matrix(c( 0,q_x,  0,
                     0,  0,q_t,
                     0,  0,  0), nrow=3, byrow = T)
    diag(qmat) = -rowSums(qmat)
    
    return(qmat)
}

# Setup for the probability evolution curves ----------------------------------
# evaluate at the mean of Z1 and separate data based on Z2 values with posterior
# median values

# seq_min = seq_max = seq_by = NULL
# if(simulation) {
#     seq_min = 0
#     seq_max = 2
#     seq_by  = 0.002
# } else {
#     seq_min = 0
#     seq_max = 20
#     seq_by  = 0.02
# }
# 
# 
# prob_evolution <- function(state_from, state_to, z2, par){
# 	
# 	probEvo = rep(NA,length(seq(seq_min, seq_max, seq_by)))
# 	
# 	P = rep(0,3)
# 	P[state_from] = 1
# 	index = 1
# 
# 	for(t in seq(seq_min, seq_max, seq_by)) {
# 		beta = par[par_index$beta]
# 		P = P %*% expm(seq_by * Q(c(0, z2), beta))
# 	
# 		probEvo[index] = sum(P[state_to]) / sum(P[c(state_from, state_to)])
# 		
# 		index = index + 1
# 	}
# 	
# 	cum_prob_evo = matrix(rep(seq(seq_min, seq_max, seq_by), 2), ncol = 2)
# 	cum_prob_evo[,2] = probEvo
# 
# 	return(cum_prob_evo)
# }
# 
# pdf_title = NULL
# 
# if(simulation){
#     if(no_s2_s3) {
#         pdf_title = paste0('Plots/probEvol_',trial_num,'_no.pdf')
#     } else {
#         pdf_title = paste0('Plots/probEvol_',trial_num,'_yes.pdf')
#     }
# } else {
#     pdf_title = paste0('Plots/probEvol_',trial_num,'_btsm.pdf')
# }
# 
# pdf(pdf_title)
# par(mfrow=c(2, 2))  # mar = c(bottom, left, top, right)
# 
# # vals = colMeans(stacked_chains)
# vals = apply(stacked_chains, 2, median)
# 
# # Finding the quantiles ----------------------------------------------------
# sample_index = sample(1:nrow(stacked_chains), size = 1000)
# quantile_mat = cbind(sample_index, 0, 0, 0, 0)
# for(ss in 1:length(sample_index)) {
#     if (ss %% 100 == 0) print(ss)
#     temp1 = prob_evolution(1, 2, 0, stacked_chains[sample_index[ss], ])
#     temp2 = prob_evolution(1, 2, 1, stacked_chains[sample_index[ss], ])
#     temp3 = prob_evolution(2, 3, 0, stacked_chains[sample_index[ss], ])
#     temp4 = prob_evolution(2, 3, 1, stacked_chains[sample_index[ss], ])
#     quantile_mat[ss,2] = sum(temp1[,2])
#     quantile_mat[ss,3] = sum(temp2[,2])
#     quantile_mat[ss,4] = sum(temp3[,2])
#     quantile_mat[ss,5] = sum(temp4[,2])
# }
# 
# if(simulation){
#     if(no_s2_s3) {
#         save(quantile_mat, file = paste0('Plots/quantile_mat',trial_num,'_no.rda'))
#     } else {
#         save(quantile_mat, file = paste0('Plots/quantile_mat',trial_num,'_yes.rda'))
#     }
# } else {
#     save(quantile_mat, file = paste0('Plots/quantile_mat',trial_num,'_btsm.rda'))
# }
# 
# # Begin the plotting -------------------------------------------------------
# 
# # Plot 1 -------------------------------------------------------------------
# plot( NULL, ylab='Probability', xlab='time', ylim=c(0,1), xlim=c(seq_min, seq_max), 
# 	  main=paste0("S1 -> S2, Z2 = 0"))
# 
# z2_1 = prob_evolution(1, 2, 0, vals)
# z2_1_true = prob_evolution(1, 2, 0, true_vals)
# 
# lines(seq(seq_min, seq_max, seq_by), z2_1[,2], lwd=2)
# if(simulation) lines(seq(seq_min, seq_max, seq_by), z2_1_true[,2], lty = 4)
# 
# small_df1 = quantile_mat[,c(1,2)]
# small_df1 = small_df1[order(small_df1[,2]),]
# upper_q = small_df1[975, 1]
# lower_q = small_df1[25, 1]
# 
# z2_1_upp = prob_evolution(1, 2, 0, stacked_chains[upper_q, ])
# z2_1_low = prob_evolution(1, 2, 0, stacked_chains[lower_q, ])
# lines(seq(seq_min, seq_max, seq_by), z2_1_upp[,2], lwd=1, col = 'blue')
# lines(seq(seq_min, seq_max, seq_by), z2_1_low[,2], lwd=1, col = 'red')
# 
# if(simulation) {
#     legend("bottomright", legend=c("Estimated", "True"), lty = c(1,4), 
#        bty="n", horiz = T)
# }
# 
# # Plot 2 -------------------------------------------------------------------
# plot( NULL, ylab='Probability', xlab='time', ylim=c(0,1), xlim=c(seq_min, seq_max), 
#       main=paste0("S1 -> S2,  Z2 = 1"))
# z2_2 = prob_evolution(1, 2, 1, vals)
# z2_2_true = prob_evolution(1, 2, 1, true_vals)
# 
# lines(seq(seq_min, seq_max, seq_by), z2_2[,2], lwd=2)
# if(simulation) lines(seq(seq_min, seq_max, seq_by), z2_2_true[,2], lty = 4)
# 
# small_df2 = quantile_mat[,c(1,3)]
# small_df2 = small_df2[order(small_df2[,2]),]
# upper_q = small_df2[975, 1]
# lower_q = small_df2[25, 1]
# 
# z2_2_upp = prob_evolution(1, 2, 1, stacked_chains[upper_q, ])
# z2_2_low = prob_evolution(1, 2, 1, stacked_chains[lower_q, ])
# lines(seq(seq_min, seq_max, seq_by), z2_2_upp[,2], lwd=1, col = 'blue')
# lines(seq(seq_min, seq_max, seq_by), z2_2_low[,2], lwd=1, col = 'red')
# 
# if(simulation) {
#     legend("bottomright", legend=c("Estimated", "True"), lty = c(1,4), 
#        bty="n", horiz = T)
# }
# 
# # Plot 3 -------------------------------------------------------------------
# plot( NULL, ylab='Probability', xlab='time', ylim=c(0,1), xlim=c(seq_min, seq_max), 
# 	  main=paste0("S2 -> S3, Z2 = 0"))
# 
# z2_1b = prob_evolution(2, 3, 0, vals)
# z2_1b_true = prob_evolution(2, 3, 0, true_vals)
# 
# lines(seq(seq_min, seq_max, seq_by), z2_1b[,2], lwd=2)
# if(simulation) lines(seq(seq_min, seq_max, seq_by), z2_1b_true[,2], lty = 4)
# 
# small_df3 = quantile_mat[,c(1,4)]
# small_df3 = small_df3[order(small_df3[,2]),]
# upper_q = small_df3[975, 1]
# lower_q = small_df3[25, 1]
# 
# z2_1b_upp = prob_evolution(2, 3, 0, stacked_chains[upper_q, ])
# z2_1b_low = prob_evolution(2, 3, 0, stacked_chains[lower_q, ])
# lines(seq(seq_min, seq_max, seq_by), z2_1b_upp[,2], lwd=1, col = 'blue')
# lines(seq(seq_min, seq_max, seq_by), z2_1b_low[,2], lwd=1, col = 'red')
# 
# if(simulation) {
#     legend("bottomright", legend=c("Estimated", "True"), lty = c(1,4), 
#        bty="n", horiz = T)
# }
# 
# # Plot 4 -------------------------------------------------------------------
# plot( NULL, ylab='Probability', xlab='time', ylim=c(0,1), xlim=c(seq_min, seq_max), 
#       main=paste0("S2 -> S3, Z2 = 1"))
# 
# z2_2b = prob_evolution(2, 3, 1, vals)
# z2_2b_true = prob_evolution(2, 3, 1, true_vals)
# 
# lines(seq(seq_min, seq_max, seq_by), z2_2b[,2], lwd=2)
# if(simulation) lines(seq(seq_min, seq_max, seq_by), z2_2b_true[,2], lty = 4)
# 
# small_df4 = quantile_mat[,c(1,5)]
# small_df4 = small_df4[order(small_df4[,2]),]
# upper_q = small_df4[975, 1]
# lower_q = small_df4[25, 1]
# 
# z2_2b_upp = prob_evolution(2, 3, 1, stacked_chains[upper_q, ])
# z2_2b_low = prob_evolution(2, 3, 1, stacked_chains[lower_q, ])
# lines(seq(seq_min, seq_max, seq_by), z2_2b_upp[,2], lwd=1, col = 'blue')
# lines(seq(seq_min, seq_max, seq_by), z2_2b_low[,2], lwd=1, col = 'red')
# 
# if(simulation) {
#     legend("bottomright", legend=c("Estimated", "True"), lty = c(1,4), 
#        bty="n", horiz = T)
# }
# 
# dev.off()
