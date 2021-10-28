##HMM-Forward-Backward algorithm

#Stationary probability calculation
stationaryProb <- function(tau){
  return(tau/sum(tau))
}

#Transition matrix is constructed
transProb <- function(tau){
  return(matrix(c(1-tau[1], tau[1], tau[2], 1-tau[2]), byrow=TRUE, nrow=2))
}


#Simulates some data for the chain
simulate <- function(TT, tau, mu, sigma2){
  #get transition matrix and stationary probabilities
  tm <- transProb(tau)
  st <- stationaryProb(tau)
  
  #hidden states are 1 for resting and 2 for active
  s <- numeric(TT)
  s[1] <- sample(1:2, 1, prob=st)
  for(i in 2:TT){
    s[i] <- sample(1:2, 1, prob=tm[s[i-1],])
  }
  
  #Data is simulated
  d <- rnorm(TT, mean=mu[s], sd=sqrt(sigma2));
  
  return(data.frame(s=s, d=d));
}

#Function for forward recursion
runForward <- function(data, tau, mu, sigma2){
  #Gets the transition matrix and stationary probabilities
  tm <- transProb(tau)
  st <- stationaryProb(tau)
  sigma <- sqrt(sigma2)
  
  #Storage for timepoints
  TT <- length(data)
  #alpha is the recursive variable, generated for each timepoint
  alpha <- matrix(data=0, nrow=TT, ncol=2)
  
  #Alpha 1 is created
  alpha[1,] <- dnorm(data[1], mu, sigma) * st 
  #Normalization
  s <- sum(alpha[1,])
  alpha[1,] <- alpha[1,] / s
  scale <- log(s)
  
  
  for(i in 2:TT){
    #Recursively calculating the future alphas. %*% is matrix multiplication
    #We get the previous alpha, times the transition matrix, times the probability of the current state (normal distribution) with the estimates for mu, sigma
    alpha[i,] <- (alpha[i-1,] %*% tm) * dnorm(data[i], mu, sigma)
    #More normalization
    s <- sum(alpha[i,])
    alpha[i,] <- alpha[i,] / s
    
    scale <- scale + log(s)    
  }
  
  return(list("alpha"= alpha, "scale"=scale));
}


#Functions for backwards recursion

runBackward <- function(data, tau, mu, sigma2, alpha){
  #Getting the transition matrix
  tm <- transProb(tau)
  sigma <- sqrt(sigma2)
  
  #Storage is prepared
  TT <- length(data)
  beta <- matrix(data=0, nrow=TT, ncol=2)
  gamma <- matrix(data=0, nrow=TT, ncol=2)
  
  #Initial beta is always 1
  beta[TT,] <- 1
  gamma[TT,] <- alpha[TT,] * beta[TT,]
  gamma[TT,] <- gamma[TT,] / sum(gamma[TT,])
  
  #Recursion starts
  for(i in (TT-1):1){
    #t(tm) is the transpose of the transition matrix. 
    beta[i,] <- (beta[i+1,] * dnorm(data[i+1], mu, sigma)) %*% t(tm)
    s <- sum(beta[i,])
    beta[i,] <- beta[i,] / s
    
    #then gamma is a normalized alpha*beta, and then normalized.
    gamma[i,] <- alpha[i,] * beta[i,]
    gamma[i,] <- gamma[i,] / sum(gamma[i,])
  }
  
  #return beta and gamma
  return(list("beta"=beta, "gamma"=gamma))
}



#Bayesian inference
estimateStatesBayesian <- function(data, tau, mu, sigma2){
  #Forward-backward recursion is run
  f <- runForward(data, tau, mu, sigma2)
  b <- runBackward(data, tau, mu, sigma2, f$alpha)
  
  #Backward object is returned
  return(b$gamma)
}

#Plotting is done here
plotStatePosteriors <- function(data, gamma, col=c('red', 'dodgerblue'), lwd=1.5){
  layout(matrix(1:2, ncol=1), heights=c(0.1,1))
  par(mar=c(0,4,0.5,1.0), oma=c(4,0,0,0), yaxs='i', xaxs='i', las=1, xpd=NA)
  
  #plot true states on top
  plot(0, type='n', xaxt='n', yaxt='n', xlab="", ylab="", xlim=c(1, length(data$s)), bty='n')
  change <- c(1, which(data$s[2:length(data$s)] != data$s[1:(length(data$s)-1)])+1, length(data$s)+1)
  for(i in 2:length(change)){
    rect(change[i-1]-0.5, par("usr")[3], change[i]-0.5, par("usr")[4], col=col[data$s[change[i-1]]], border=NA)
  }
  
  #plot posterior probabilities
  plot(gamma[,1], type='l', col=col[1], xlab="Time", ylab="Posterior Probability", lwd=lwd, ylim=c(-0.02,1.02))
  lines(gamma[,2], type='l', col=col[2], lwd=lwd)
}


#Everything is run

TT <- 100
tau <- c(0.05, 0.05)
##If these are set to be even, the algorithm can't distinguish the states
mu <- c(10, 140)
sigma2 <- 1000

data <- simulate(TT, tau, mu, sigma2)
gamma <- estimateStatesBayesian(data$d, tau, mu, sigma2)
plotStatePosteriors(data, gamma)


