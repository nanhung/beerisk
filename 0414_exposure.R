if(!require(fitdistrplus)) {
  install.packages("fitdistrplus"); 
  require(fitdistrplus)
} 

require(gridExtra)
library(ggplot2)
library(mc2d)
library(EnvStats)

dat <- read.table("exposure.dat", head = T)
EU_non0 <- subset(dat, adj>0.1)
fitln_EU <- fitdist(EU_non0$adj,"lnorm",method="mle")
boot.conc <- bootdist(fitln_EU, niter = 1001)
summary(boot.conc)

## fitdistrplus
meanlog_EU <- fitln_EU$est[1]
sdlog_EU <- fitln_EU$est[2]
summary(fitln_EU)
plot(fitln_EU)

Boot_EU <- bootdist(fitln_EU, bootmethod="param", niter=1001)

## mc2d ----
# Specification of uncertainty and variability dimensions
ndvar(1001)
ndunc(1001)

# Uncertainty on the distribution parameters
Mean_conso_EU <- mcdata(Boot_EU$estim$meanlog, type="U")
Sd_conso_EU <- mcdata(Boot_EU$estim$sdlog, type="U")

conso1_EU <- mcstoc(rlnorm, type="VU", meanlog= Mean_conso_EU, sdlog= Sd_conso_EU, 
                    rtrunc=TRUE, linf=0, lsup=10000)
conso0 <- mcdata(0, type="V")
plot(conso1_EU)

## data manipulate
pzero_eu <- sum(dat$adj < 0.2)/length(dat$adj)
v_EU <- mcprobtree(c(pzero_eu,1-pzero_eu), list("0"=conso0,"1"=conso1_EU), type = "VU")
v_EU<-as.numeric(v_EU)
dose<-data.frame(rep(deparse(substitute(v_EU)), length(v_EU)), v_EU)
colnames(dose) <- c("name","dose")

#plot exposure dose
data <- c(unname(quantile(conso1_EU, c(0.6))), unname(quantile(conso1_EU, c(0.7))), unname(quantile(conso1_EU, c(0.8))), unname(quantile(conso1_EU, c(0.95))))
percentile <- c("60 Percentile", "70 Percentile", "80 Percentile", "95 Percentile")
df1 <- data.frame(data, percentile)

Ex <- ggplot(dose, aes(x=dose, fill=name))+
  geom_histogram(aes(y=..density..), alpha=0.8, colour="darkgrey", fill="grey", 
                 position="identity", bins = 50)+
  scale_x_log10(breaks=c(1,10,100,1000,10000), limits = c(0.01, 100000))+
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))+
  ylab("Density")+xlab(bquote('Dose ('*mu*g~L^-1*')'))+theme(legend.position = "none")+
  geom_density(fill="white", alpha = 0.05)+
  ggtitle("Simulated exposure dose")+
  geom_vline(df1, mapping=aes(xintercept=data), color="red") 

# 161220
bee.pop <- function(t, state, parms) {
  f = state [1]
  B = state [2]
  H = state [3]
  F = state [4]
  
  with(as.list(c(parms)), 
       {
         # compute time-dependent parameter
         c1 <- c0*(0.448-0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
         l1 <- l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365))
         m1 <- m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365))
         
         # compute toxicity
         a <- 3941.3673; n <- 2.0033; c <- 0.3476 
         md <- c*d^n/(a^n+d^n)
         
         #Function for survival of brood
         S <- (f^2/(f^2+b^2))*(H/(H+v))
         #recruitment function
         R <- amin-s*F/(H+F)+amax*(b^2/(b^2+f^2))
         
         # compute derivatives
         df <- c1*F-ga*(H+F)-gb*B
         dB <- l1*S-fi*B
         dH <- fi*B-R*H
         dF <- R*H-(m1+md)*F
         
         # combine results
         list(c(df, dB, dH, dF))
       })
}

# set parameters
parms <- function(x){
  c(c0=0.1, l0=2000, amin=0.25, amax=0.25, s=0.75,m0=0.1, 
    ga=0.007, gb=0.018, fi=0.1, b=500, v=5000, d=unname(quantile(conso1_EU, c(x))))
}


# set states
state = c(f = 3000, B = 2000, H = 2000, F = 1000)

# set output times
times = seq(0, 1000, by = 1)

# ode output
out.5 <- ode(y = state, times = times, func = bee.pop, parms = parms(.5))
out.75 <- ode(y = state, times = times, func = bee.pop, parms = parms(.75))
out.9 <- ode(y = state, times = times, func = bee.pop, parms = parms(.9))
out.975 <- ode(y = state, times = times, func = bee.pop, parms = parms(.975))

## Makes the data in a more ggplot-friendly format
stacker <- function(df){
  df.stack <- stack(df[, -1])
  df.stack$time <- rep(seq_len(nrow(df)), length(table(df.stack$ind)))
  names(df.stack)[1:2] <- c("population", "compartments")
  df.stack$compartments <- factor(df.stack$compartments, 
                                  levels=c("f", "B","H","F"), ordered=TRUE)
  return(df.stack)
  
}

gg.dynamic.plot <- function(x,y){
  # Data manipulate
  guess <- stacker(as.data.frame(ode(y = state, times = times, func = bee.pop, parms = parms(x))))
  
  # Base plot
  ggplot(guess, aes(x=time, y=population, 
                    group=compartments, 
                    color=compartments)) + 
    geom_line(aes(colour = compartments), size=1.5) + 
    xlab(" ") + ylab(" ") +
    theme(legend.position="none") + theme(legend.title=element_blank()) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
    scale_y_continuous(lim=c(0, y))
}

E60 <- gg.dynamic.plot(.60, 80000)+ggtitle("Exposure dose at 60 perentile") + ylab("Food Storagrs and Bee Population")
E70 <- gg.dynamic.plot(.70, 80000)+ggtitle("Exposure dose at 70 perentile") 
E80 <- gg.dynamic.plot(.80, 80000)+ggtitle("Exposure dose at 80 perentile") + xlab("Time (day)") + ylab("Food Storagrs and Bee Population")
E95 <- gg.dynamic.plot(.95, 20000)+ggtitle("Exposure dose at 95 perentile") + xlab("Time (day)")

x11(8, 11)
grid.arrange(Ex, arrangeGrob(E60, E70,ncol=2),arrangeGrob(E80, E95,ncol=2), ncol=1)
