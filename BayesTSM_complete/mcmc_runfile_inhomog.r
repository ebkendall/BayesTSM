source("mcmc_routine.r")

args = commandArgs(TRUE)

seed_num = as.numeric(args[1]) 
set.seed(seed_num)
print(seed_num)

ind = 1
trial_num = 1

simulation = T
no_s2_s3 = F
disc = T

if(simulation) {
    init_par = c(  3,  0.5, 0.5, 0.5, # beta(0, x), beta(1, x), beta(2, x), beta(2, x_time)
                 1.2,  0.5, 0.5, 0.5) # beta(0, t), beta(1, t), beta(2, t), beta(2, t_time)
} else {
    init_par = c(-3, 0, 0, 0,
                 -3, 0, 0, 0)
}

par_index = list( beta=1:8 )

# Defining the mean and variance for the uninformed Gaussian priors for the MH update
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

if(simulation) {
    if(no_s2_s3) {
        if(disc) {
            print("A")
            load(paste0("Data/hmm_sim_extended", ind, "_no_inhomog.rda"))
            temp_data = as.matrix(rawData); rownames(temp_data) = NULL
        } else {
            print("B")
            load(paste0("Data/hmm_sim", ind, "_no_inhomog.rda"))
            temp_data = as.matrix(rawData); rownames(temp_data) = NULL
        }
    } else {
        if(disc) {
            print("C")
            load(paste0("Data/hmm_sim_extended", ind, "_yes_inhomog.rda"))
            temp_data = as.matrix(rawData); rownames(temp_data) = NULL
        } else {
            print("D")
            load(paste0("Data/hmm_sim", ind, "_yes_inhomog.rda"))
            temp_data = as.matrix(rawData); rownames(temp_data) = NULL
        }
    }
} else {
    if(disc) {
        print("E")
        load(paste0("Data/bayestsm_dat_extended", ind, ".rda"))
        temp_data = as.matrix(bayestsm_dat); rownames(temp_data) = NULL
    } else {
        print("F")
        load(paste0("Data/bayestsm_dat", ind, ".rda"))
        temp_data = as.matrix(bayestsm_dat); rownames(temp_data) = NULL
    }
}

id = temp_data[, "id"]
y  = temp_data[, "state"]
x  = temp_data[, c("Z.1", "Z.2"), drop=F]
t  = temp_data[, "time"]

disc_t = rep(-1, length(t))
if(disc) disc_t = temp_data[,"disc_time"]

steps = 20000
burnin = 5000

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, disc, disc_t)

e_time = Sys.time() - s_time; print(e_time)

if(simulation) {
    if(no_s2_s3) {
        if(disc) {
            print("A")
            save(mcmc_out, file = paste0("Model_out/mcmc_out_extended_", seed_num, "_", trial_num, "_no_inhomog.rda"))
        } else {
            print("B")
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", seed_num, "_", trial_num, "_no_inhomog.rda"))
        }
    } else {
        if(disc) {
            print("C")
            save(mcmc_out, file = paste0("Model_out/mcmc_out_extended_", seed_num, "_", trial_num, "_yes_inhomog.rda"))
        } else {
            print("D")
            save(mcmc_out, file = paste0("Model_out/mcmc_out_", seed_num, "_", trial_num, "_yes_inhomog.rda"))
        }
    }
} else {
    if(disc) {
        print("E")
        save(mcmc_out, file = paste0("Model_out/mcmc_out_extended_", seed_num, "_", trial_num, "_btsm_inhomog.rda"))
    } else {
        print("F")
        save(mcmc_out, file = paste0("Model_out/mcmc_out_", seed_num, "_", trial_num, "_btsm_inhomog.rda"))
    }
}
