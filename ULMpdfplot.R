library(ggplot2)
library(rmutil)


#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")

## Uniform-Laplace Mixture##

p  <- 0.5

l1 <- function(x) p*dunif(x,-5,5) +(1-p)*dlaplace(x,0,1)
l2 <- function(x) p*dunif(x,-6,6)+(1-p)*dlaplace(x,0,1)
l3 <- function(x) p*dunif(x,-12,12)+(1-p)*dlaplace(x,0,1)
l4 <- function(x) p*dunif(x,-18,18)+(1-p)*dlaplace(x,0,1)

ggplot2::ggplot(data.frame(x = c(-20, 20)), aes(x)) +theme_bw()+
              stat_function(fun=l1, aes(linetype = "l1" ))+
              stat_function(fun=l2,aes(  linetype = "l2"))+
              stat_function(fun = l3,aes( linetype = "l3"))+
              stat_function(fun = l4,aes( linetype = "l4"))+
              scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash"),name="a",label = c("a = 5", "a = 6", "a = 12","a = 18"))+
               xlab("x")+ylab("pdf")+
              ggtitle(expression(paste("ULM ( ", a, " ,", sigma, " = 1,",p,"= 0.5)")))+
              theme(plot.title = element_text(size = 10,hjust = 0.5))+
              guides(fill = guide_legend(override.aes = list(color = NA)),color = FALSE)
            
            

##Uniform exponential Mixture##

g1 <- function(x) p*dunif(x,0,5) +(1-p)*dexp(x,1)
g2 <- function(x) p*dunif(x,0,6)+(1-p)*dexp(x,1)
g3 <- function(x) p*dunif(x,0,12)+(1-p)*dexp(x,1)
g4 <- function(x) p*dunif(x,0,18)+(1-p)*dexp(x,1)

 ggplot(data.frame(x=c(0, 22)), aes(x)) +theme_bw()+
  stat_function(fun=g1, aes(linetype = "g1"))+
  stat_function(fun=g2,aes(linetype = "g2"))+
  stat_function(fun = g3,aes( linetype = "g3"))+
  stat_function(fun = g4,aes( linetype = "g4"))+
  scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash"), name=" a",label = c("a = 5", "a = 6", "a = 12", "a = 18"))+
  xlab("y")+ylab("pdf")+ggtitle(expression(paste("UEM ( ", a, ", ", sigma,
                                                 " = 1,", p," = 0.5)")))+
  theme(plot.title = element_text(size = 10, hjust = 0.5))
