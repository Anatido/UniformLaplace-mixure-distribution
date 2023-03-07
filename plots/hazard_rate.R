library(tidyverse)
library(scales)
library(cowplot)

hazardrate <- function(y, sig = 5, a = 10, p = 0.1){
if(y<=a){
(1/sig)*((sig*p + a*(1-p)*exp(-y/sig))/(p*(a-y) + a*(1-p)*exp(-y/sig)))
}else{
1/sig
}
}



yval <- seq(0, 60, 0.1)
p_val = sapply(yval, hazardrate, sig = .5, a = 5, p = 0.5)

plot(yval, p_val, "l")
hazardval_df <- data.frame(yvalue = yval, value = sapply(yval, hazardrate, sig = 5, a = 5, p = 0.2), parray = rep("0.2", length(yval)))

for (i in p_val){
hazardval_df <- rbind(hazardval_df, cbind(data.frame(value = sapply(yval, hazardrate, sig = 5, a = 5, p = p_val[i])), data.frame(parray = rep(as.character(p_val[i]), length(yval)))))
}

param_df <- rbind(param_df, cbind(as.data.frame(iter_param(df = samDat, a = 0.5, k = 3, p = pval[i], sig = simVal[i], 
  maxiter = maxiteration)$estimates), data.frame(Initialization = 
rep(paste("p = ",as.character(pval[i]), "", "sigma = ", as.character(simVal[i])), maxiteration))))

yval1 <- seq(0, 4, 0.001)
hazardval1 <- data.frame(yvalue = yval1, value = sapply(yval1, hazardrate, sig = 5, a = 2, p = 0.30), p = rep("0.2", length(yval1)))
hazardval1_0.5 <- data.frame(yvalue = yval1, value = sapply(yval1, hazardrate, sig = 2, a = 2, p = 0.5), p = rep("0.5", length(yval1)))
hazardval1_0.8 <- data.frame(yvalue = yval1, value = sapply(yval1, hazardrate, sig = 2, a = 2, p = 0.8), p = rep("0.8", length(yval1)))

hazardval_df1 <- rbind(hazardval1_0.2, hazardval1_0.5, hazardval1_0.8)

H1 <- ggplot(hazardval1 ,aes(x = yvalue, y = value))+
  geom_point(size = 0.1)+ 
  geom_point(aes(x = 2, y = max(value) ),size = 3)+
  geom_point(aes(x = 2, y = 1/5 ),size = 3, shape = 1)+
  xlab("y") + ylab("hazard rate")+ theme_bw() +
  theme(legend.position = c(.9, .7)) 


yval2 <- seq(0, 15, 0.001)
hazardval2 <- data.frame(yvalue = yval2, value = sapply(yval2, hazardrate, sig = 10, a = 10 , p = 0.50), p = rep("0.2", length(yval2)))
hazardval2_0.5 <- data.frame(yvalue = yval2, value = sapply(yval2, hazardrate, sig = 10, a = 5, p = 0.5), p = rep("0.5", length(yval2)))
hazardval2_0.8 <- data.frame(yvalue = yval2, value = sapply(yval2, hazardrate, sig = 10, a = 5, p = 0.8), p = rep("0.8", length(yval2)))

hazardval_df2 <- rbind(hazardval2_0.2, hazardval2_0.5, hazardval2_0.8)



H2 <- ggplot(hazardval2 ,aes(x = yvalue, y = value))+
  geom_point(size = 0.1)+ 
  geom_point(aes(x = 10, y = max(value) ),size = 3)+
  geom_point(aes(x = 10, y = 1/10 ),size = 3, shape = 1)+
  xlab("y") + ylab("hazard rate")+ theme_bw() +
  theme(legend.position = c(.9, .7)) 

yval3 <- seq(0, 8, 0.0001)
hazardval3_0.2 <- data.frame(yvalue = yval3, value = sapply(yval3, hazardrate, sig = 1, a = 5 , p = 0.1), p = rep("0.2", length(yval3)))
hazardval3_0.5 <- data.frame(yvalue = yval3, value = sapply(yval3, hazardrate, sig = 0.5, a = .5, p = 0.5), p = rep("0.5", length(yval3)))
hazardval3_0.8 <- data.frame(yvalue = yval3, value = sapply(yval3, hazardrate, sig = 0.5, a = .5, p = 0.8), p = rep("0.8", length(yval3)))

hazardval_df3 <- rbind(hazardval3_0.2, hazardval3_0.5, hazardval3_0.8)



H3 <- ggplot(hazardval3_0.2 ,aes(x = yvalue, y = value))+
  geom_point(size = 0.1)+ 
  geom_point(aes(x = 5, y = max(value) ),size = 3)+
  geom_point(aes(x = 5, y = 1/1 ),size = 3, shape = 1)+
    xlab("y") + ylab("hazard rate")+ theme_bw() +
  theme(legend.position = c(.9, .7)) 

ggplot(hazardval3_0.2 ,aes(x = yvalue, y = value))+
  geom_point(size = 0.1)+ 
  geom_point(aes(x = 10, y = max(value) ),size = 3)+
  geom_point(aes(x = 10, y = 1/2 ),size = 3, shape = 1)+
  xlab("y") + ylab("hazard rate")+ theme_bw() +
  theme(legend.position = c(.9, .7)) 
plot_grid(H1, H2, H3, nrow = 1)

hazardrate(1, a = 10, p = 0.2, sig = 0.1)

yval3 <- seq(1.5, 8, 0.005)
ggplot(hazardval_df3 ,aes(x = yvalue, y = value, color = p)) + geom_point()
  geom_line(aes(linetype = p), size  = 1) + 
  scale_x_continuous(breaks= pretty_breaks()) +
  xlab("y") + ylab("hazard rate")+ theme_bw() +
  theme(legend.position = c(.9, .7))
