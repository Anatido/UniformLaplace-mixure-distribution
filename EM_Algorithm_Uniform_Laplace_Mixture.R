library(rmutil)
library(microbenchmark)


emOptimizerULM <- function(df, p0, sig0, maxiter = 100, tol = 1e-5){
  
  df <- sort(df)
  n <- length(df) 
  
  for (i in 1:n ){
    if (df[i] < 0){ 
      stop("Error: A data point is negative. Take absolute of data to make it postive")
    }
  }
  
  if (p0 < 0 || sig0 < 0){
    stop(" Error: Initial paramater should be positive ")
  }
  
  
  ## initialize dataframe to store estimates in subspaces
  estim <- array(NA_real_, dim = c(n+1,4), 
                 dimnames = list(NULL, NULL))
  
  ## When a < y1
  
  ahat = min(df)/2
  p_oo = 0
  sig_oo = mean(df)
  lo = -n*(log(sig_oo))-1
  
  estim[1,1] <- ahat
  estim[1,2] <- p_oo 
  estim[1,3] <- sig_oo
  estim[1,4] <- lo
  
  
  
  ## When k = n 
  p_n <- p0
  sig_n <- sig0
  
  ln <- sum(log(p_n/max(df) + ((1-p_n)/sig_n)*exp(-df/sig_n)))
  diff_n = 1
  iter_n = 0
  
  while (diff_n > tol & iter_n < maxiter){
    
    ## E-step
    expD_n <- (1/sig_n)*exp(-df/sig_n)
    num_n <- p_n/max(df)
    cexpect_n <- num_n/(num_n + (1-p_n)*expD_n)
    
    ## M-step: 
    p.new_n <- sum(cexpect_n)/n
    sig.new_n <- (sum(df*cexpect_n) - n*mean(df))/(sum(cexpect_n) - n)
    
    #oberseved likelihood
    likeli_n <- sum(log( p.new_n/max(df) + ((1-p.new_n)/sig.new_n)*exp(-df/sig.new_n)))
    
    ## checking convergence
    diff_n <- likeli_n - ln
    ln <- likeli_n
    p_n <- p.new_n
    sig_n <- sig.new_n
    
    
    iter_n = iter_n + 1
    
  }
  
  estim[n+1,2] <- p_n
  estim[n+1,3] <- sig_n
  estim[n+1,4] <- ln 
  
  
  ## subspace k = 1,..n-1 
  for (k in (1:n-1)){
    p<- p0
    sig <- sig0
    likel0 <- sum(log(p/max(df[1:k]) + ((1-p)/sig)*exp(-df[1:k]/sig))) - (n-k)*log(sig) + (n-k)*log(1-p)-sum(df[(k+1):n]/sig)
    
    diff = 1
    iter = 0
    
    while (diff > tol & iter < maxiter){
      
      ## E-step
      expD <- (1/sig)*exp(-df[1:k]/sig)
      num <- p/max(df[1:k])
      cexpect <- num/(num + (1-p)*expD)
      
      ## M-step: 
      p.new <- sum(cexpect)/n
      sig.new <- (sum(df[1:k]*cexpect) - n*mean(df))/(sum(cexpect) - n)
      
      #oberseved likelihood
      likeli <- sum(log(p.new/max(df[1:k]) + ((1-p.new)/sig.new)*exp(-df[c(1:k)]/sig.new))) - (n-k)*log(sig.new) + (n-k)*log(1-p.new)-sum(df[(k+1):n]/sig.new)
      
      ## checking convergence
      diff <- likeli - likel0
      likel0 <- likeli
      p <- p.new
      sig <- sig.new
      
      iter = iter + 1;
      
    }
    
    estim[2:(n+1),1] <- df
    estim[k+1,2] <- p
    estim[k+1,3] <- sig
    estim[k+1,4] <- likel0             
  }
  
  colnames(estim) <- c("a_hat", "p_hat", "sig_hat", "log-likelihood")
  optimalEm = subset(estim, estim[, 4] == max(estim[, 4]))[1,]
  optimalEm
}

# Example

sampData <- runifLaplace(a = 3, p = .4, sigma = 0.60, n = 500)
emOptimizer(sampData, p0 = 0.5, sig0 = mean(abs(sampData)))
