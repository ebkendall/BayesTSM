library(expm)
source("mcmc_routine.r")

args = commandArgs(TRUE)

seed_num = as.numeric(args[1]) 
set.seed(seed_num)
print(seed_num)

# *************************************************************************
# disc = TRUE  --> piece-wise time homogeneous solution
# disc = FALSE --> true continuous time solution
disc = FALSE
# *************************************************************************

init_par = c(-3, 0, 0, 0, # beta(0, x), beta(1, x), beta(2, x), beta(2, x_time)
             -3, 0, 0, 0) # beta(0, t), beta(1, t), beta(2, t), beta(2, t_time)

par_index = list( beta=1:8 )

# Defining the mean and variance for the uninformed Gaussian priors for the MH update
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

if(disc) {
    load(paste0("Data/bayestsm_dat_extended", seed_num, ".rda"))
    temp_data = as.matrix(bayestsm_dat); rownames(temp_data) = NULL
} else {
    load(paste0("Data/bayestsm_dat", seed_num, ".rda"))
    temp_data = as.matrix(bayestsm_dat); rownames(temp_data) = NULL
    # *************************************************************************
    # SCALING TIME SO ALL TIMES ARE < 1 
    temp_data[, "time"] = temp_data[, "time"] / 100
    # *************************************************************************
}

id = temp_data[, "id"]
y  = temp_data[, "state"]
x  = temp_data[, c("Z.1", "Z.2"), drop=F]
t  = temp_data[, "time"]

disc_t = rep(-1, length(t))
if(disc) disc_t = temp_data[,"disc_time"]

# *************************************************************************
# Feel free to change steps and burnin (I would not go above 10,000 for burnin)
# *************************************************************************
steps = 20000
burnin = 5000

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, disc, disc_t)

e_time = Sys.time() - s_time; print(e_time)

if(disc) {
    print("E")
    save(mcmc_out, file = paste0("Model_out/mcmc_out_extended_", seed_num, "_btsm_inhomog.rda"))
} else {
    print("F")
    save(mcmc_out, file = paste0("Model_out/mcmc_out_", seed_num, "_btsm_inhomog.rda"))
}

# ------------------------------------------------------------------------------
# Passage Probabilities --------------------------------------------------------
# ------------------------------------------------------------------------------

pdf_title = NULL

if(disc) {
    pdf_title = paste0('Plots/probEvol_',seed_num,'_extended_btsm.pdf')
} else {
    pdf_title = paste0('Plots/probEvol_',seed_num,'_btsm.pdf')
}

pdf(pdf_title)
par(mfrow=c(2, 1)) 

# ------------------------------------------------------------------------------
# Load the dataset we are investigating ----------------------------------------
# ------------------------------------------------------------------------------
if(disc) {
    load(paste0("Data/bayestsm_dat_extended", seed_num, ".rda"))
    myData = bayestsm_dat
} else {
    load(paste0("Data/bayestsm_dat", seed_num, ".rda"))
    myData = bayestsm_dat
    # *************************************************************************
    # SCALING TIME SO ALL TIMES ARE < 1 
    myData$time = myData$time / 100
    # *************************************************************************
}

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
        t_seq = seq(0, 0.50, length.out=100)
    } else {
        t_seq = seq(0, 0.50, length.out=100)
    }
    
    probEvo = rep(NA,length(t_seq))
    
    P = rep(0,3)
    P[state_from] = 1
    index = 1
    
    for(t in 2:length(t_seq)) {
        beta = par
        new_P = hmm_solver_P(par, c(0,t_seq[t]), c(z1, z2))
        P_t = P %*% new_P
        
        if(state_to == 2) {
            probEvo[index] = P_t[state_to] + P_t[state_to+1]   
        } else {
            probEvo[index] = P_t[state_to]
        }
        
        index = index + 1
    }
    
    return(probEvo)
}

prob_evolution_expm <- function(state_from, state_to, z1, z2, par) {
    
    # *************************************************************************
    # BECAUSE THIS USES THE PIECE-WISE TIME HOMOGENEOUS ASSUMPTION, THE t_seq 
    # USED NEEDS TO BE AT LEAST AT THE DISCRETIZATION LEVEL OF THE UPDATES FOR
    # THE TRANSITION RATE MATRIX 
    # (reference simulate_data_inhomog.r for more info on the discretization level)
    # *************************************************************************

    t_seq = NULL
    if(state_to == 2) {
        t_seq = seq(0, 100, by = 1)
    } else {
        t_seq = seq(0, 100, by = 1)
    }
    
    probEvo = rep(NA,length(t_seq))
    
    P = rep(0,3)
    P[state_from] = 1
    index = 1
    
    for(t in 2:length(t_seq)) {
        beta = par
        Q_mat = Q(c(z1, z2, t_seq[t-1]), par)
        new_P = expm((t_seq[t] - t_seq[t-1]) * Q_mat)
        P = P %*% new_P
        
        if(state_to == 2) {
            probEvo[index] = P[state_to] + P[state_to+1]   
        } else {
            probEvo[index] = P[state_to]
        }
        
        index = index + 1
    }
    
    return(probEvo)
}

if(disc) {
    t.x = seq(0, 100, by = 1)
    t.y = seq(0, 100, by = 1)
} else {
    t.x = seq(0, 0.50, length.out=100) # seq(0, 100, length.out=100)
    t.y = seq(0, 0.50, length.out=100) # seq(0, 20, length.out=100)
}

# create an empty matrix to store CDFs for each unique covariate combination
unique_ids = unique(myData[,'id'])
cdf.x = cdf.y = matrix(nrow = length(unique_ids), ncol = length(t.x)) 

par = colMeans(mcmc_out$chain)
for(i in 1:nrow(cdf.x)) {
    print(i)
    subDat = myData[myData[,'id'] == unique_ids[i], ]
    z_1 = subDat[1,'Z.1']
    z_2 = subDat[1,'Z.2']
    if(disc) {
        cdf.x[i,] = prob_evolution_expm(1, 2, z_1, z_2, par)
        cdf.y[i,] = prob_evolution_expm(2, 3, z_1, z_2, par)
    } else {
        cdf.x[i,] = prob_evolution_hmmsolver(1, 2, z_1, z_2, par)
        cdf.y[i,] = prob_evolution_hmmsolver(2, 3, z_1, z_2, par)
    }
}


###  State 1 to 2
plot(t.x,colMeans(cdf.x),type="l",col=1,lty=1,lwd=2,ylim=c(0,1),
     xlab="time since entry in state 1",ylab="CDF",
     main = "Transition from state 1 to 2")

###  State 2 to 3
plot(t.y,colMeans(cdf.y),type="l",col=1,lty=1,lwd=2,ylim=c(0,1),
     xlab="time since entry in state 2",ylab="CDF",
     main = "Transition from state 2 to 3")
dev.off()
