mleOptimizer <- function(df){
  df <- sort(df)
  n <- length(df) 
  
  for (i in 1:n ){
    if (df[i] < 0){ 
      stop("Error: A data point is negative. Take absolute of data to make it postive")
    }
  }
  
  
  df <- sort(df)
  n <- length(df)
  hist <- array(NA_real_, dim = c(length(df),4), 
                dimnames = list(NULL, NULL))
  hist2 <- array(NA_real_, dim = c(1,4), 
                 dimnames = list(NULL, NULL))
  
  lo<- -n*log(mean(df))-n
  
  for (i in 1: length(df)){ 
    k <- i
    lik1 <- function(x){
      P <- x[1]
      sigma <- x[2]
      a <- sum(log(P/df[k] + (1-P)/sigma*exp(-df[1:k]/sigma)))
      b <- -(n-k)*log(sigma)+(n-k)*log(1-P)-sum(df[-(1:k)]/sigma)
      -(a+b)
    }
    opk<-optim(c(P= 0.2,sigma = 1),lik1)
    
    if(opk$par[1]<0){
      phat<-0
      sighat<-mean(df)
      lks<- -n*log(mean(df))-n
    }else{
      phat<-opk$par[1]
      sighat<-opk$par[2]
      lks<-sum(log(opk$par[1]/df[k] + (1-opk$par[1])/opk$par[2]*exp(-df[1:k]/opk$par[2])))-(n-k)*log(opk$par[2])+(n-k)*log(1-opk$par[1])-sum(df[-(1:k)]/opk$par[2])
    }
    
    hist[, 1] <- df
    hist[i,2] <- phat
    hist[i,3] <- sighat
    hist[i,4] <- lks
    hist2[1,1] <- min(df)/2
    hist2[1,2] <- 0
    hist2[1,3] <- mean(df)
    hist2[1,4] <- lo
  }
  
  mles = rbind(round(hist2,4), round(hist,4))
  colnames(mles) <- c("a_hat", "p_hat", "sig_hat", "log-likelihood")
  
  optimalmle <-  subset(mles, mles[, 4] == max(mles[, 4]))[1,]
  optimalmle
}

sampData <- runifLaplace(a = 3, p = .4, sigma = 0.60, n = 500)
mleOptimizer(abs(sampData))

