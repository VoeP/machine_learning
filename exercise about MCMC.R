##machine learning exercises: MCMC

#Function to simulate data
simData <- function(n, m, sigma_2=1, sigma_2_mu=1){
  mu ~ N(0,1)
  mu <- rnorm(n, 0, sqrt(sigma_2_mu))
  
  #simulate x_ij ~ N(mu_i, sigma_2)
  x <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:n)
    x[i,] <- rnorm(m, mu[i], sqrt(sigma_2))
  
  return(list(x=x, sigma_2=sigma_2, mu=mu))
}


#Functions for calculation of MLE

MLE_mu <- function(x){
  #needs to return vector of length equal to x
  m <- ncol(x)
  return(1/m * rowSums(x))
}

MLE_sigma2 <- function(x, mu){
  m <- ncol(x)
  n <- nrow(x)
  return((1/(n*m)) * sum(rowSums((x - mu)^2)))
}


##Functions to run MCMC

runMCMC <- function(x, mu_1, sigma2_1, propSd_Sigma2=0.1, propSd_mu=0.01, len=2000){
  n <- nrow(x)
  
  #matrices for storage
  mu <- matrix(0, nrow = len, ncol = n)
  mu[1,] <- mu_1
  sigma2 <- numeric(len)
  sigma2[1] <- sigma2_1
  
  #storage for normal density for each data point
  lh_old <- apply(x, 2, dnorm, mu[1], sqrt(sigma2[1]))
  
  #MCMC is run here
  for(t in 2:len){ 

    newMus <- rnorm(n, mean=mu[t-1,], sd=propSd_mu)
    
    #Hastings ratio
    lh_new <- apply(x, 2, dnorm, newMus, sqrt(sigma2[t-1]))
    #prod is for i
    h <- apply(lh_new / lh_old, 1, prod) 
    
    #get accepted values
    random <- runif(n)
    accepted <- which(random < h)
    mu[t, accepted] <- newMus[accepted]
    lh_old[accepted,]  <- lh_new[accepted,]
    
    #get rejected values
    rejected <- which(random >= h)
    mu[t,rejected] <- mu[t-1, rejected]
    

    #functions for updating sigma 2
    
    #propose sigmas
    newSigma2 <- abs(rnorm(1, mean=sigma2[t-1], sd=propSd_Sigma2))
    
    #calculation of Hastings ratio
    lh_new <- apply(x, 2, dnorm, as.numeric(mu[t,]), sqrt(newSigma2))
    #prod is for i and j
    h <- prod(lh_new / lh_old)
    
    #accept or reject
    if(runif(1) < h){ # accept
      sigma2[t] <- newSigma2
      lh_old <- lh_new
    } else { 
      sigma2[t] <- sigma2[t-1]
    }
  }
  
  return(list(sigma2=sigma2, mu=mu))
}


#Plotting the traces

plotOneTrace <- function(chain, trueValue, MLE, ylim, col='orange2', ylab="Parameter", writeLegend=FALSE){
  #This is the trace
  plot(chain, type='l', ylim=ylim, xlab="Iteration", ylab=ylab, col=col)
  #This is the true mu
  abline(h = trueValue, lty=2, col='black')
  #MLE estimate of mu
  abline(h = MLE, lty=3, col='firebrick')
  
  #acceptance rate calculation
  len <- length(chain)
  acc <- sum(chain[1:(len-1)] != chain[2:len])/(len-1)
  
  #designing plots
  y <- par("usr")[4] - 0.1 * diff(par("usr")[3:4])
  text(len, y, labels=paste(round(100*acc, digits=2), "%", sep=""), pos=2)
  if(writeLegend){
    legend("topleft", c("truth", "MLE"), col = c("black", "firebrick"), lty = c(2, 3), bty = "n", cex=0.8)
  }
}

plotTraces <- function(chain, whichMuToPlot, true_mu, true_sigma2, MLE_mu, MLE_sigma2, lim_mu=c(-1,1), lim_sigma2=c(0,0.5)){
  par(xaxs='i', yaxs='i', mfrow=c(2,2), mar=c(3.5,3.75,0.5,1), mgp=c(2.25,0.66, 0), las=1)
  for (i in whichMuToPlot){
    plotOneTrace(chain$mu[,i], true_mu[i], MLE_mu[i], ylim=lim_mu, ylab=bquote(Mean ~ mu[.(i)]))
  }

  plotOneTrace(chain$sigma2, true_sigma2, MLE_sigma2, ylim=lim_sigma2, ylab=expression(paste("Variance ", sigma^2)), col="dodgerblue", writeLegend=TRUE)
}



#Functions for plotting the posteriors
plotOnePosterior <- function(chain, trueValue, MLE, xlim, col='orange2', xlab="Parameter", writeLegend=FALSE){
  #plot posterior density
  dens <- density(chain)
  plot(dens, type='l', xlim=xlim, xlab=xlab, ylab="Posterior density", main="", col=col)
  #add line for true value
  abline(v = trueValue, lty=2, col='black')
  #add line for MLE
  abline(v = MLE, lty=3, col='firebrick')
  
  # write legend
  if(writeLegend){
    legend("topleft", c("truth", "MLE"), col = c("black", "firebrick"), lty = c(2, 3), bty = "n", cex=0.8)
  }
}

plotPosteriors <- function(chain, whichMuToPlot, true_mu, true_sigma2, MLE_mu, MLE_sigma2, lim_mu=c(-1,1), lim_sigma2=c(0,0.5)){
  #open plot
  par(xaxs='i', yaxs='i', mfrow=c(2,2), mar=c(3.5,3.75,0.5,1), mgp=c(2.25,0.66, 0), las=1)
  
  #plot posterior for three mu
  for (i in whichMuToPlot){
    print(i)
    plotOnePosterior(chain$mu[,i], true_mu[i], MLE_mu[i], xlim=lim_mu, xlab=bquote(Mean ~ mu[.(i)]))
  }
  
  #plot posterior for sigma2
  plotOnePosterior(chain$sigma2, true_sigma2, MLE_sigma2, xlim=lim_sigma2, xlab=expression(paste("Variance ", sigma^2)), col="dodgerblue", writeLegend=TRUE)
}

##Data simulation and plotting

n <- 1000
m <- 2
data <- simData(n, m, sigma_2=2, sigma_2_mu=3)


#MCMC parameters

#nr iterations
len=1000
#standard deviation for proposal kernel of mu
propSd_Sigma2=0.1 
#standard deviation for proposal kernel of sd
propSd_mu=2.2 

#MLE calculation
mu_hat <- MLE_mu(data$x)
sigma2_hat <- MLE_sigma2(data$x, mu_hat)
print(paste("true sigma2 = ", 1, " MLE sigma2 = ", sigma2_hat, sep = ""))

#MLE is starting point for MCMC
chain <- runMCMC(data$x, mu_hat, sigma2_hat, propSd_Sigma2=propSd_Sigma2, propSd_mu=propSd_mu, len=len)

#plotting
theseMu <- sample(length(data$mu), 3)
lim_mu <- range(data$mu)*1.5
lim_sigma2 <- c(0, 1.5*data$sigma_2)

plotTraces(chain, theseMu, data$mu, data$sigma_2, mu_hat, sigma2_hat, lim_mu, lim_sigma2)

#Marginal posterior densities are plotted
plotPosteriors(chain, theseMu, data$mu, data$sigma_2, mu_hat, sigma2_hat, lim_mu, lim_sigma2)