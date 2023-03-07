
emOptimizerUGM <- function(df, p0, sig0, maxiter = 100, tol = 1e-5){
  
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
  
  #initialize empty frame
  estim <- array(NA_real_, dim = c(n+1,4), 
                 dimnames = list(NULL, NULL))
  
  ## When a < y1
  ahat = min(df)/2
  p_oo = 0
  sig_oo = mean(df^2)
  lo = -n*log(sig_oo) + n*log((sqrt(2/pi))) - n/2
  
  estim[1,1] <- ahat
  estim[1,2] <- p_oo 
  estim[1,3] <- sig_oo
  estim[1,4] <- lo
  
  # When k = n
  p_n <- p0
  sig_n <- sig0
  
  ln <- sum(log(p_n/max(df) + sqrt(2/pi)*((1-p_n)/sig_n)*exp(-0.5*(df/sig_n)^2)))
  diff_n1 = 1
  diff_n2 = 1
  iter_n  = 0
  
  while (all(c(diff_n1, diff_n2) > tol)  & (iter_n < maxiter)){
    
    ## E-step
    expD_n <- (1/(sig_n))*sqrt(2/pi)*exp(-0.5*(df/sig_n)^2)
    num_n <- p_n/max(df)
    cexpect_n <- num_n/(num_n + (1-p_n)*expD_n)
    
    ## M-step: 
    p.new_n <- sum(cexpect_n)/n
    sig.new_n <- ((sum((df^2)*cexpect_n) - n*mean(df^2))/(sum(cexpect_n) - n))^0.5
   
    #oberseved likelihood
    likeli_n <- sum(log(p.new_n/max(df) + sqrt(2/pi)*((1-p.new_n)/sig.new_n)*exp(-0.5*(df/sig.new_n)^2)))
    
    ## checking convergence
    #diff_n <- abs(likeli_n - ln)
    diff_n1 <- abs(p.new_n - p_n)
    diff_n2 <- abs(sig.new_n - sig_n)
    
    #ln <- likeli_n
    p_n <- p.new_n
    sig_n <- sig.new_n
    
    #expD_est <- (1/sig)*exp(-df[1:k]/sig)
    #num_est <- p/max(df[1:k])
    #cexpect_est <- num/(num + (1-p)*expD)
    #E_est <- sum(cexpect_est*log(p) + cexpect_est*log(1/max(df[1:k])) + (1- cexpect_est)*log(1-p) + (1- cexpect_est)*log(expD_est)) + sum(log(1-p)+log((1/sig)*exp(-df[(k+1):n]/sig)))
    
    iter_n = iter_n + 1
  }
  
  
  estim[n+1,2] <- p_n
  estim[n+1,3] <- sig_n
  estim[n+1,4] <- likeli_n 
  
  
  ## k = 1,..n-1 
  for (k in (1:n-1)){
    p <- p0
    sig  <- sig0
    likel0 <- sum(log(p/max(df[1:k]) + sqrt(2/pi)*((1-p)/sig)*exp(-0.5*(df[1:k]/sig)^2))) - (n-k)*log(sig) + (n-k)*log((1-p)*sqrt(2/pi))-(1/(2*sig^2))*sum((df[(k+1):n])^2)
    
    diff1 =  1
    diff2 =  1
    iter = 0
    
    while (all(c(diff1, diff2) > tol) & iter < maxiter){
      
      ## E-step
      expD <- ((1/(sig))*sqrt(2/pi))*exp(-0.5*(df[1:k]/sig)^2)
      num <- p/max(df[1:k])
      cexpect <- num/(num + (1-p)*expD)
      
      ## M-step: 
      p.new <- sum(cexpect)/n
      sig.new <- ((sum((df[1:k]^2)*cexpect) - n*mean(df^2))/(sum(cexpect) - n))^0.5
      
      #oberseved likelihood
      likeli <- sum(log(p.new/max(df[1:k]) + sqrt(2/pi)*((1-p.new)/sig.new)*exp(-0.5*(df[1:k])^2/sig.new^2))) - (n-k)*log(sig.new) + (n-k)*log((1-p.new)*sqrt(2/pi))-(1/(2*sig.new^2))*sum((df[(k+1):n])^2)
      
      ## checking convergence
      #diff <- likeli - likel0
      #likel0 <- likeli
      
      diff1 <- abs(p.new - p)
      diff2 <- abs(sig.new - sig)
      p <- p.new
      sig <- sig.new
      
      #expD_est <- (1/sig)*exp(-df[1:k]/sig)
      #num_est <- p/max(df[1:k])
      #cexpect_est <- num/(num + (1-p)*expD)
      #E_est <- sum(cexpect_est*log(p) + cexpect_est*log(1/max(df[1:k])) + (1- cexpect_est)*log(1-p) + (1- cexpect_est)*log(expD_est)) + sum(log(1-p)+log((1/sig)*exp(-df[(k+1):n]/sig)))
      
      iter = iter + 1;
      # cat("Iter", iter, ": p =", p.new, ", sig =",sig.new, "diff = ", diff, "\n")
    }
    estim[2:(n+1),1] <- df
    estim[k+1,2] <- p
    estim[k+1,3] <- sig
    estim[k+1,4] <- likeli            
  }
  
  colnames(estim)<-c("a_hat", "p_hat", "sig_hat", "likelihood")
  return(subset(estim, estim[, 4]== max(estim[, 4])))
}

