library(plotly)

# liklihood values
loglikelihood_values <- function(a, p, sigma, n0, k, k_optimal = FALSE, sigma_interval = 0.5, sigma_end = 5){
  p_values <- as.double(seq(0.0001, .99, 0.1))[-1]
  sig_values <- as.double(seq(0.0001, sigma_end, sigma_interval))[-1]
  grid_values <- array(NA_real_, dim = c(length(sig_values), length(p_values)))
  set.seed(56)
  ## data
  sim_data <- function(a, p, sigma, n0){
    ## Sample N random bernoulli Random variables
    B = rbinom(n0, 1, p)
    # Variable to store the samples from the mixture distribution
    rand.samples = rep(NA, n0)
    #Sampling from the mixture
    for(i in 1:n0){
      if(B[i]==1){
        rand.samples[i] = runif(1,0,a)
      }else{
        rand.samples[i] = rexp(1, sigma)
      }
    }
    return(rand.samples)
  }
  df <- sort(sim_data(a, p, sigma, n0))
  n <- length(df)
  if(k_optimal == FALSE){
    if(k < n & k != 1){
      for( i in 1: length(sig_values)){
        for( j in 1: length(p_values)){
          likel0 <- sum(log(p_values[j]/max(df[1:k]) + ((1-p_values[j])/sig_values[i])*exp(-df[1:k]/sig_values[i]))) - (n-k)*log(sig_values[i]) + (n-k)*log(1- p_values[j])- sum(df[(k+1):n]/sig_values[i])
          grid_values[i, j] <- likel0
        }
      }
    }else if( k >= n ){
      for( i in 1: length(sig_values)){
        for( j in 1: length(p_values)){
          likel0 <- sum(log(p_values[j]/max(df) + ((1-p_values[j])/sig_values[i])*exp(-df/sig_values[i])))
          grid_values[i, j] <- likel0
        }
      }
    }else{
      for( i in 1: length(sig_values)){
        for( j in 1: length(p_values)){
          likel0 = n*log(1-p_values[j])-n*log(sig_values[i]) - sum(df)/sig_values[i]
          grid_values[i, j] <- likel0
        }
      }
    }
  }else{
    k_optim <- max(which(a >= df))
    if(k_optim != -Inf & k_optim < n  ){
      for( i in 1: length(sig_values)){
        for( j in 1 : length(p_values)){
          likel0 <- sum(log(p_values[j]/max(df[1:k_optim]) + ((1 - p_values[j])/sig_values[i])*exp(-df[1:k_optim]/sig_values[i]))) -
            (n - k_optim)*log(sig_values[i]) + (n - k_optim)*log(1- p_values[j])- sum(df[(k_optim + 1):n]/sig_values[i])
          grid_values[i, j] <- likel0
        }
      }
    }else if(k_optim != -Inf & k_optim >= n){
      for( i in 1: length(sig_values)){
        for( j in 1: length(p_values)){
          likel0 <- sum(log(p_values[j]/max(df) + ((1- p_values[j])/sig_values[i])*exp(-df/sig_values[i])))
          grid_values[i, j] <- likel0
        }
      }
    }else{
      for( i in 1: length(sig_values)){
        for( j in 1: length(p_values)){
          likel0 = n*log(1-p_values[j])-n*log(sig_values[i]) - sum(df)/sig_values[i]
          grid_values[i, j] <- likel0
        }
      }
    }
  }
  list(pp = p_values, ss = sig_values, grid_values = grid_values)
}



# 3D plot
plot_3d_function <- function(a = 6,  p = .2, sigma = 1/5, n0 = 30,
                             k = 15, k_optimal = FALSE , sigma_interval = 0.5,
                             sigma_end = 8, scene){
  likelihood_values <- loglikelihood_values(a, p, sigma , n0 ,
                                            k, k_optimal, sigma_interval, sigma_end )
  tt <-loglikelihood_values(a,  p , sigma , n0 ,
                            k , k_optimal  , sigma_interval ,
                            sigma_end )
  if(k_optimal == FALSE){
    axx <- list(
      title = "p"
    )
    axy <- list(
      title = "sigma"
    )
    axz <- list(
      title = "log-likelihood"
    )
    plot_ly(z = tt$grid_values,x = tt$pp, y = tt$ss, scene = scene)%>% add_surface() %>%
      layout(scene = list(xaxis= axx,yaxis= axy,zaxis= axz))
    #plot_ly(z = likelihood_values$grid_values , x = likelihood_values$pp, y = likelihood_values$ss )  %>% add_surface() %>%
    # layout(scene = list(xaxis= axx,yaxis= axy,zaxis= axz))
  }else{
    axx <- list(
      title = "p"
    )
    axy <- list(
      title = "sigma"
    )
    axz <- list(
      title = "log-likelihood"
    )
    plot_ly(z = tt$grid_values,x = tt$pp, y = tt$ss, scene = scene)%>% add_surface() %>%
      layout(scene = list(xaxis= axx,yaxis= axy,zaxis= axz))
    #plot_ly(z = likelihood_values$grid_values, x =likelihood_values$pp, likelihood_values$ss) %>% add_surface() %>%
    # layout(scene = list(xaxis= axx,yaxis= axy,zaxis= axz))
  }
}
