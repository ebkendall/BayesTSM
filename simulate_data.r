args <- commandArgs(TRUE)
no_s2_s3 = args[1] # takes value TRUE or FALSE


# Construct the transition rate matrix
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

true_vals = c(  3,  0.5, 0.5, # beta(0, x), beta(1, x), beta(2, x)
              1.2,  0.5, 0.5) # beta(0, t), beta(1, t), beta(2, t)
            #  -4, -4)          # P(obs S3 | S2), P(obs S2 | S3)

par_index = list( beta=1:6) #, misclass=7:10,  pi_logit=13:14)

N = 1000
# Choose the discretization of time.
dt <- 1/200

# errorMat_temp = matrix(c(1, 0, 0,
#                          0, 1, exp(true_vals[par_index$misclass][1]),
#                          0, exp(true_vals[par_index$misclass][2]), 1), ncol=3, byrow=TRUE)

# errorMat = errorMat_temp / rowSums(errorMat_temp)
errorMat = diag(3)

# Initial state probabilities
# initProbs_temp = c( 1, exp(true_vals[par_index$pi_logit][1]), exp(true_vals[par_index$pi_logit][2]))
# initProbs = initProbs_temp / sum(initProbs_temp)
initProbs = c(1,0,0)

# 10 replicates of the simulated data
for(ttt in 1:10) {
  set.seed(ttt)
  print(ttt)

  # Learning interobservation times
  inter_obs_time = NULL
  load(paste0('Data/bayestsm_dat', ttt, '.rda'))
  for(i in unique(bayestsm_dat$id)) {
    subject = bayestsm_dat[bayestsm_dat$id == i, , drop = F]
    inter_obs_time = c(inter_obs_time, diff(subject$time))
  }
  
  inter_obs_time = inter_obs_time / 200
  X_and_S = matrix(nrow = N, ncol = 2)
  colnames(X_and_S) = c("X", "S")
  
  rawData <- NULL
  NumObs_sim <- NULL
  for(i in 1:N){
    
    # z1 sample from N(0,1)
    z1 = rnorm(1)
    # z2 sample from Bernoulli(0.5)
    z2 = rbinom(1, 1, 0.5)

    # Sample for an initial state.
    trueState <- sample(1:3, size=1, prob=initProbs)

    # Sample the remaining states until death.
    time_main <- 0
    time1 <- 0
    s <- trueState
    while(s < 3){

      # Infinitesimal transition rates.
      qmat <- Q(c(z1, z2),true_vals[par_index$beta])

      # Possible next states.
      moveToStates <- which(qmat[s,] > 0)

      # Sample the wait times before transition to each of the next possible states.
      waitTimes <- rexp( n=length(moveToStates), rate= qmat[s,moveToStates])

      # If any of the wait times are smaller than dt, then transition to the state with the minimum wait time.
      min_waitTime <- min(waitTimes)
      if(min_waitTime < dt){  s <- moveToStates[ which(waitTimes == min_waitTime) ]  }

      time1 <- time1 + dt

      time_main <- c( time_main, time1)
      trueState <- c( trueState, s)
    }
    timeOfDeath <- tail(time_main,1)

    # Calculating the observed transition times
    X_and_S[i, ] = diff(c(0,time_main[which(diff(trueState) == 1) + 1]))
    
    # Sample inter-observation times from the cav data set.  Maximum of 20 time_main in study.
    visitTimes <- NULL
    time2 <- 0

    while(time2 < min( 0.15, timeOfDeath)){

      visitTimes <- c( visitTimes, time2)
      time2 <- time2 + sample( inter_obs_time, size=1)
    }

    # If first visit time occurs after death, then subject is NOT entered into the study.
    if( !is.null(visitTimes) ){

      # If death occured before the study ended, then record the time of death.
      if( timeOfDeath < 0.15 ){  visitTimes <- c( visitTimes, timeOfDeath) }

      state <- NULL
      for(k in 1:length(visitTimes)){  state <- c( state, tail( trueState[ time_main <= visitTimes[k] ], 1))  }
      
      temp_dat = cbind(state, visitTimes)
        
      if(no_s2_s3) {
          if((2 %in% state) & (3 %in% state)) {
              s2s3 = bayestsm_dat$state[bayestsm_dat$state != 1]
              if(runif(n=1, min = 0, max = 1) > mean(s2s3 == 3)) {
                  # state 2 is terminal state
                  terminal_index = min(which(state == 2))
                  temp_dat = temp_dat[1:terminal_index, , drop = F]
              } else {
                  # state 3 is terminal state
                  temp_dat = temp_dat[state != 2, , drop = F]
                  terminal_index = min(which(temp_dat[,"state"] == 3))
                  temp_dat = temp_dat[1:terminal_index, , drop = F]
              }
          }
      }
      
      if(nrow(temp_dat) > 1 & temp_dat[1,"state"] == 1) {
        n_i <- nrow(temp_dat)
        ptnum <- rep(i,n_i)
        time_main <- temp_dat[,"visitTimes"]
        true_state <- temp_dat[,"state"]
        rawData <- rbind( rawData, data.frame(ptnum,time_main,z1, z2,true_state) )

        NumObs_sim <- c( NumObs_sim, n_i)
      } else {
        print(paste0(i, ": not correct"))
      }
    }
  }
  
  colnames(rawData) <- c('id','time','Z.1', 'Z.2', 'state')

  temp_raw = t(table(rawData$id))
  print(summary(c(temp_raw)))
  
  temp_bayes = t(table(bayestsm_dat$id))
  print(summary(c(temp_bayes)))
  
  # # Add noise to the states.
  # for(i in 1:nrow(rawData)){	rawData$state[i] <- sample(1:4, size=1, prob=errorMat[rawData$state[i],])  }
  
  if(no_s2_s3) {
      save(rawData, file = paste0('Data/hmm_sim', ttt, '_no.rda'))
      save(X_and_S, file = paste0('Data/X_and_S', ttt, '_no.rda'))
  } else {
      save(rawData, file = paste0('Data/hmm_sim', ttt, '_yes.rda'))
      save(X_and_S, file = paste0('Data/X_and_S', ttt, '_yes.rda'))
  }
}
