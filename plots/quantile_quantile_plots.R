library(tibble)
library(ggplot2)

# Quantile function
# Laplace

qunifLaplace <- function(u, a = estimates[1], p = estimates[2], sig = estimates[3]){
  fa1 <- (1-p)/2 * exp(-a/sig)
  fa2 <- p + (1-p)/2 *(2-exp(-a/sig))
  
  if(u > 0 && u <= fa1){
    r <- sig*log(2*u/(1-p))
  }else if(u >= fa1 && u <= 0.5){
    g1 <- function(x) p*(x+a)/(2*a) + (1-p)/2 * exp(x/sig) - u
    r <-uniroot(g1, c(-a, 0))$root
    
  }else if(u >= 0.5 && u <= fa2){
    g2 <- function(x) p*(x+a)/(2*a) + (1-p)/2 * (2 - exp(-x/sig)) - u
    r <-uniroot(g2, c(0, a+300))$root
    
  }else{
    r <- -sig*log(2*(1-u)/(1-p))
  }
  return(r)
}


# normal
qunifnormal <- function(u, a = estimates_norm[1] , p = estimates_norm[2], sig = estimates_norm[3]){
  fa1 <- (1-p) * pnorm(-a, mean = 0, sd = sig)
  fa2 <- p + (1-p)*pnorm(a, mean = 0, sd = sig)
  
  if(u > 0 && u <= fa1){
    g1 <- function(x) (1-p)*pnorm(x, mean = 0, sd = sig)-u
    r <- uniroot(g1, c(-300, 100))$root
  }else if(u >= fa1 && u <= fa2){
    g2 <- function(x) p*(x+a)/(2*a) + (1-p)*pnorm(x, 0, sig) - u
    r <-uniroot(g2, c(-a-10, a+10))$root
  }else{
    g3 <- function(x) p + (1-p)*pnorm(x, 0, sig) - u
    r <-uniroot(g3, c(a, a+500))$root
  }
  return(r)
}



## dat  <- runifLaplace()
## dat <- runifnormal()

df <- tibble(
  x_sample =  dat, 
  q = (rank(x_sample ) - .5) /length(x_sample),
  sample_quant =  quantile(x_sample, probs = q),
  x_theoretical = sapply(q, qunifLaplace)
)

qqlap = ggplot(df) + 
  geom_point(aes(x = x_theoretical, y = sample_quant)) +
  labs(
    x = "Theoretical quantiles", 
    y = "Sample quantiles",
    title = "ULM"
  ) +
  geom_abline(intercept = 0, slope = 1, col = "red")+
  theme_bw()

df <- tibble(
  x_sample = dat, 
  q = (rank(x_sample ) - .5) /length(x_sample),
  sample_quant =  quantile(x_sample, probs = q),
  x_theoretical = sapply(q, qunifnormal)
)

qqnormal = ggplot(d) + 
  geom_point(aes(x = x_theoretical, y = sample_quant)) +
  labs(
    x = "Theoretical quantiles", 
    y = "Sample quantiles",
    title = "UGM"
  ) +
  geom_abline(intercept = 0, slope = 1, col = "red")+
  theme_bw()

cowplot::plot_grid(qqlap, qqnormal)