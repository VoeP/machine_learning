#Analytical posterior

posterior <- function(lambda, data_obs, alpha, beta){
  alpha_prime <- alpha + length(data_obs)
  beta_prime <- beta + sum(data_obs)
  
  return(dgamma(lambda, alpha_prime, beta_prime))
}


##Functions for ABC

simulate <- function(m, lambda){
  N <- length(lambda)
  return(matrix(rexp(N*m, rate=lambda), nrow=N, byrow=FALSE))
}

#Summary statistic calculation
calcStats <- function(D){
  S <- data.frame("mean" = apply(D, 1, mean),
                  "var"  = apply(D, 1, var),
                  "median" = apply(D, 1, median),
                  "min"    = apply(D, 1, min),
                  "max"    = apply(D, 1, max),
                  "range"  = apply(D, 1, max) - apply(D, 1, min))
  return(S)
}

#Function for generating simulations
makeSims <- function(N, m, alpha, beta){
  #Parameters are drawn from the prior (gamma distribution is used)
  lambda <- rgamma(N, alpha, beta)
  
  #data is simulated using the normal distribution
  D <- simulate(m, lambda)
  
  #Using the simulated data, we calculate summary stats here
  S <- calcStats(D)
  
  #return alist with the lambdas as well as the summary statistics
  return(list(P=as.matrix(lambda), S=S))
}

#Functions to get the posterior stats
getABCPosteriorSamples <- function(Sims, sobs, tolerance){
  #match returns a vector of the positions of matches of its first arguments in its second
  kept_stats <- match(names(sobs), names(Sims$S))
  S <- Sims$S[,kept_stats, drop=FALSE]
  
  #standardisation of stats
  for(i in 1:ncol(S)){
    mean <- mean(S[,i])
    sd <- sd(S[,i])
    S[,i] <- (S[,i] - mean) / sd
    sobs[i] <- (sobs[i] - mean) / sd
  }
  
  #Euclidian distance is calculated, t() takes the matrix transpose
  x <- t(S) - as.numeric(sobs)
  distances <- sqrt(colSums(x*x))
  
  #Rejection
  distancesSorted <- sort(distances, index.return=T)
  #best simulations are chosen
  bestSims <- distancesSorted$ix[1:round(tolerance*nrow(Sims$P))]
  return(list(P=Sims$P[bestSims,], S=S[bestSims,]))
}


#Functions for plotting

#Plotting posterior distributions
plotPosteriors <- function(d, sobs, alpha, beta, Sims, delta, trueLambda, col=c("dodgerblue", "purple", "orange2", "red"), main=""){
  par(mar=c(4.2,4,1.2,0.2))
  
  #Calculating smoothed ABC posterior densities
  abcDensity <- list(length(delta))
  xlim <- rep(trueLambda, 2)
  ymax <- 0
  for(i in 1:length(delta)){
    abcSamples <- getABCPosteriorSamples(Sims, sobs, delta[i]) 
    abcDensity[[i]] <- density(abcSamples$P) 
    
    xlim <- range(c(abcDensity[[i]]$x, xlim))
    ymax <- max(ymax, abcDensity[[i]]$y)
  }
  
  #Preparation of grid
  lambda <- seq(xlim[1], xlim[2], length.out=1000)
  analyticalDensity <- posterior(lambda, d, alpha, beta)
  
  #Analytical posteriors are plotted
  ylim <- c(0, max(analyticalDensity, ymax))
  plot(lambda, analyticalDensity, xlim=xlim, ylim=ylim, type='l', xlab=expression(lambda), ylab="Posterior density", main=main)
  lines(lambda, dgamma(lambda, alpha, beta), lty=2)
  
  abline(v=trueLambda, lty=2)
  
  #ABC posteriors are plotted
  for(i in 1:length(delta)){
    lines(abcDensity[[i]], col=col[i])
  }
  
  legend <- c("Prior", "Analytical",  paste("ABC", delta))
  legend('topright', lwd=1, col=c('black', 'black', col), lty=c(2, 1, rep(1, length(delta))), legend=legend, bty='n')
}



#Functions are run

#True values
trueLambda <- 2
m <- 10

#Parameters for prior
alpha <- 5
beta <- 1

#ABC settings
N <- 100000
toleranceDelta <- c(0.5, 0.1, 0.05, 0.01)

#Simulations are generated
Sims <- makeSims(N, m, alpha, beta)

#Observations are generated
d <- simulate(m, trueLambda)
sobs <- calcStats(d)

# Effect of num simulations (smoothing)
par(mfrow=c(1,2))
plotPosteriors(d, sobs, alpha, beta, Sims, toleranceDelta, trueLambda, main=paste("N =", N))
Sims_short <- list(P=Sims$P[1:1000,, drop=FALSE], S=Sims$S[1:1000,])
plotPosteriors(d, sobs, alpha, beta, Sims_short, toleranceDelta, trueLambda, main="N = 1000")

# Effect fo summary statistics choice
par(mfrow=c(1,3))
plotPosteriors(d, sobs, alpha, beta, Sims_short, toleranceDelta, trueLambda, main="All statistics")
plotPosteriors(d, sobs[,1,drop=FALSE], alpha, beta, Sims, toleranceDelta, trueLambda, main="Only mean")
plotPosteriors(d, sobs[,3,drop=FALSE], alpha, beta, Sims, toleranceDelta, trueLambda, main="Only min")
