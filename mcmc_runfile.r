source("mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

trial_num = 3
simulation = T
no_s2_s3 = T

# x_i is the transition time from S1 -> S2
# t_i/s_i is the transition time from S2 -> S3
# In an AFT, they are functions of the covariates, but we will instead use
# the covariates to determine the transition rates which directly correspond to
# the transition times.

if(simulation) {
    init_par = c(  3,  0.5, 0.5, # beta(0, x), beta(1, x), beta(2, x)
                 1.2,  0.5, 0.5) # beta(0, t), beta(1, t), beta(2, t)
                    #  -6, -6,          # P(obs S1 | S2), P(obs S1 | S3)
                    #  -6, -6,          # P(obs S2 | S1), P(obs S2 | S3)
                    #  -6, -6,          # P(obs S3 | S1), P(obs S3 | S2)
                    #  -6, -6)          # P(init S2), P(init S3)
} else {
    init_par = c(-3, 0, 0,
                 -3, 0, 0)
}

par_index = list( beta=1:6) # , misclass=7:12, pi_logit=13:14)

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

if(simulation) {
    if(no_s2_s3) {
        load(paste0("Data/hmm_sim", ind, "_no.rda"))
        temp_data = as.matrix(rawData); rownames(temp_data) = NULL
    } else {
        load(paste0("Data/hmm_sim", ind, "_yes.rda"))
        temp_data = as.matrix(rawData); rownames(temp_data) = NULL
    }
} else {
    load(paste0("Data/bayestsm_dat", ind, ".rda"))
    temp_data = as.matrix(bayestsm_dat); rownames(temp_data) = NULL
}
id = temp_data[, "id"]
y  = temp_data[, "state"]
x  = temp_data[, c("Z.1", "Z.2"), drop=F]
t  = temp_data[, "time"]

steps = 100000
burnin = 5000

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin)

e_time = Sys.time() - s_time; print(e_time)

if(simulation) {
    if(no_s2_s3) {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_no.rda"))
    } else {
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_yes.rda"))
    }
} else {
    save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, "_", trial_num, "_btsm.rda"))
}
