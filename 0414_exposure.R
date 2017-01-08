if(!require(fitdistrplus)) {
  install.packages("fitdistrplus"); 
  require(fitdistrplus)
} 

require(gridExtra)
require(dplyr)
library(ggplot2)
library(mc2d)
library(EnvStats)
require(deSolve)

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
ndvar(101)
ndunc(1001)

# Uncertainty on the distribution parameters
Mean_conso_EU <- mcdata(Boot_EU$estim$meanlog, type="U")
Sd_conso_EU <- mcdata(Boot_EU$estim$sdlog, type="U")

conso1_EU <- mcstoc(rlnorm, type="VU", meanlog= Mean_conso_EU, sdlog= Sd_conso_EU, 
                    rtrunc=TRUE, linf=0, lsup=10000)
conso0 <- mcdata(0, type="V")
plot(conso1_EU)
conso1_EU_val <- unclass(conso1_EU)

## data manipulate
pzero_eu <- sum(dat$adj < 0.2)/length(dat$adj)
v_EU <- mcprobtree(c(pzero_eu,1-pzero_eu), list("0"=conso0,"1"=conso1_EU), type = "VU")
v_EU<-as.numeric(v_EU)
dose<-data.frame(rep(deparse(substitute(v_EU)), length(v_EU)), v_EU)
colnames(dose) <- c("name","dose")

#plot exposure dose
data <- c(unname(quantile(conso1_EU_val, c(0.5))), unname(quantile(conso1_EU_val, c(0.75))), unname(quantile(conso1_EU_val, c(0.9))), unname(quantile(conso1_EU_val, c(0.95))))
percentile <- c("50 Percentile", "75 Percentile", "90 Percentile", "95 Percentile")
df1 <- data.frame(data, percentile)

Ex <- ggplot(dose, aes(x=dose, fill=name))+
  geom_histogram(aes(y=..density..), alpha=0.8, colour="darkgrey", fill="grey", 
                 position="identity", bins = 50)+
  scale_x_log10(breaks=c(1,10,100,1000,10000), limits = c(0.01, 100000))+
  theme(legend.title=element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  ylab("Probability")+xlab(bquote('Dose ('*mu*g~L^-1*')'))+theme(legend.position = "none")+
  geom_density(fill="white", alpha = 0.05)+
  ggtitle("Simulated exposure dose")+
  theme_bw()+
  geom_vline(df1, mapping=aes(xintercept=data), color= c("green","lightyellow3","orange","red")) 

# 0107 Risk
slp <- mcstoc(rnorm, type="VU", 4.672e-05, 1.804e-06)
int <- mcstoc(rnorm, type="VU", -6.049e-04, 3.400e-03)
Imi <- mcstoc(rlnorm, type="VU", meanlog= Mean_conso_EU, sdlog= Sd_conso_EU,
              rtrunc=TRUE, linf=0, lsup=10000)
md <- slp*Imi

# T = 10 day
Mor_1 <- (1-exp(-md*10)/1)*100
plot(Mor_1, xlab="%Mortality", ylab="Cumulative probability")
mormodel <- mc(slp, md, Imi, Mor_1)
tor<-tornado(mormodel)
plot(tor)

# T = 30 day
Mor_2 <- (1-exp(-md*30)/1)*100
plot(Mor_2, xlab="%Mortality", ylab="Cumulative probability")
mormodel <- mc(slp, md, Imi, Mor_2)
tor<-tornado(mormodel)
plot(tor)


pre.ecdf <- ecdf(Mor_1[,1,])
curve(1-pre.ecdf(x), col="red", xlim=c(0,100), ylim=c(0,1),
      xlab="10-day% mortality", ylab = "Exceedence risk")

for (i in 2:1001){
  post.ecdf <- ecdf(Mor_1[,i,])
  curve(1-post.ecdf(x), col="blue", add=TRUE)
}

pre.ecdf <- ecdf(Mor_2[,1,])
curve(1-pre.ecdf(x), col="red", xlim=c(0,100), ylim=c(0,1),
      xlab="30-day% mortality", ylab = "Exceedence risk")

for (i in 2:1001){
  post.ecdf <- ecdf(Mor_2[,i,])
  curve(1-post.ecdf(x), col="darkblue", add=TRUE)
}

# Exceedance risk
library(reshape2)
library(plyr)
melt1 <- suppressMessages(melt(Mor_1[,,]))
melt2 <- suppressMessages(melt(Mor_2[,,]))
melt1$Var1<-NULL
melt2$Var1<-NULL
names(melt1) <- c('v', 'x')
names(melt2) <- c('v', 'x')
ecdf1 <- ddply(melt1, c("v"), mutate, ecdf = ecdf(x)(unique(x))*length(x))
ecdf2 <- ddply(melt2, c("v"), mutate, ecdf = ecdf(x)(unique(x))*length(x))
ecdf_1 <- ddply(ecdf1 , "v", mutate, 
                ecdf =scale(ecdf,center=min(ecdf),scale=diff(range(ecdf))))
ecdf_2 <- ddply(ecdf2 , "v", mutate, 
                ecdf =scale(ecdf,center=min(ecdf),scale=diff(range(ecdf))))
R1 <- ggplot(ecdf_1, aes(x,1-ecdf, color = v)) + geom_point() + theme_bw() +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title=element_text(size=12))+
  xlab('10-day %mortality') + ylab('Exceedance risk') + geom_hline(yintercept = c(0, 1), linetype = "dashed", color = 'black')+
  scale_x_continuous(lim=c(0, 100))+
  geom_hline(yintercept = c(0.5,0.25,0.1,0.05), color= c("green","lightyellow3","orange","red")) 
R2 <- ggplot(ecdf_2, aes(x,1-ecdf, color = v)) + geom_point() + theme_bw() +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title=element_text(size=12))+
  xlab('30-day %mortality') + ylab('Exceedance risk') + geom_hline(yintercept = c(0, 1), linetype = "dashed", color = 'black')+
  scale_x_continuous(lim=c(0, 100))+
  geom_hline(yintercept = c(0.5,0.25,0.1,0.05), color= c("green","lightyellow3","orange","red")) 


#
probs <- c(0.5, 0.75, 0.9, 0.95)
quant1 <- as.data.frame(t(apply(Mor_1[,,], 2, quantile, probs = probs)))
quant2 <- as.data.frame(t(apply(Mor_2[,,], 2, quantile, probs = probs)))
cols <- c('green', 'lightyellow3', 'orange', 'red')
quant.melt1 <- suppressMessages(melt(quant1))
quant.melt2 <- suppressMessages(melt(quant2))
names(quant.melt1) <- c('quantile', 'x')
names(quant.melt2) <- c('quantile', 'x')
ecdf2.1 <- ddply(quant.melt1, c("quantile"), mutate, ecdf = ecdf(x)(unique(x))*length(x))
ecdf2.2 <- ddply(quant.melt2, c("quantile"), mutate, ecdf = ecdf(x)(unique(x))*length(x))
ecdf_2.1 <- ddply(ecdf2.1 , "quantile", mutate, 
                ecdf =scale(ecdf,center=min(ecdf),scale=diff(range(ecdf))))
ecdf_2.2 <- ddply(ecdf2.2 , "quantile", mutate, 
                ecdf =scale(ecdf,center=min(ecdf),scale=diff(range(ecdf))))
R3 <- ggplot(ecdf_2.1, aes(x,ecdf, color = quantile)) + geom_step(size=1) + theme_bw() +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title=element_text(size=12))+
  xlab('10-day %mortality') + ylab('Cumulative probability') + geom_hline(yintercept = c(0, 1), linetype = "dashed", color = 'black')+
  scale_x_continuous(lim=c(0, 100)) + labs(colour = "Exposue quantile") + scale_colour_manual(values = cols)
R4 <- ggplot(ecdf_2.2, aes(x,ecdf, color = quantile)) + geom_step(size=1) + theme_bw() +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title=element_text(size=12))+
  xlab('30-day %mortality') + ylab('Cumulative probability') + geom_hline(yintercept = c(0, 1), linetype = "dashed", color = 'black')+
  scale_x_continuous(lim=c(0, 100)) + labs(colour = "Exposue quantile") + scale_colour_manual(values = cols)

#
x11(8, 11)
grid.arrange(Ex, arrangeGrob(R1, R2,ncol=2), arrangeGrob(R3, R4,ncol=2), ncol=1)

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
         ylag <- ifelse((t - 0) <= 0, 2000, lagvalue(t - 0)[2])
         df <- c1*F-ga*(H+F)-gb*B
         dB <- l1*S-fi*B
         dH <- fi*ylag-R*H
         dF <- R*H-(m1+md)*F
         
         # combine results
         list(c(df, dB, dH, dF))
       })
}


# set parameters
parms <- function(x){
  c(c0=0.1, l0=2000, amin=0.25, amax=0.25, s=0.75,m0=0.1, 
    ga=0.007, gb=0.018, fi=0.1, b=500, v=5000, d=unname(quantile(conso1_EU_val, c(x))))
}

# set states
state = c(f = 3000, B = 2000, H = 2000, F = 1000)

# set output times
times = seq(0, 1000, by = 1)

# ode output
out.5 <- dede(y = state, times = times, func = bee.pop, parms = parms(.5))
out.75 <- dede(y = state, times = times, func = bee.pop, parms = parms(.75))
out.9 <- dede(y = state, times = times, func = bee.pop, parms = parms(.9))
out.95 <- dede(y = state, times = times, func = bee.pop, parms = parms(.95))

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
  guess <- stacker(as.data.frame(dede(y = state, times = times, func = bee.pop, parms = parms(x))))
  
  # Base plot
  ggplot(guess, aes(x=time, y=population, 
                    group=compartments, 
                    color=compartments)) + 
    geom_line(aes(colour = compartments), size=1.5) + 
    xlab(" ") + ylab(" ") +
    theme(legend.title=element_blank()) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
    scale_y_continuous(lim=c(0, y))
}

E50 <- gg.dynamic.plot(.50, 90000)+ggtitle("Exposure dose at 50 perentile") + ylab("Food Storagrs and Bee Population") + theme(legend.position="none")
E75 <- gg.dynamic.plot(.75, 90000)+ggtitle("Exposure dose at 75 perentile") + theme(legend.position="none")
E90 <- gg.dynamic.plot(.90, 90000)+ggtitle("Exposure dose at 90 perentile") + xlab("Time (day)") + ylab("Food Storagrs and Bee Population") + theme(legend.position="none")
E95 <- gg.dynamic.plot(.95, 20000)+ggtitle("Exposure dose at 95 perentile") + xlab("Time (day)") + theme(legend.position=c(.8, .7))

x11(8, 11)
grid.arrange(Ex, arrangeGrob(E50, E75,ncol=2),arrangeGrob(E90, E95,ncol=2), ncol=1)







# Find quantile
quantile(md, c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95))

# Parameters
l0 <- 2000; v<-5000; amax<-0.25; amin<-0.25; s<-0.75; m0<-0.14; b<-500 
ra <- 0.007; rb<-0.018; fi<-0.11; c0 <- 0.15

t <- 1:365

#
m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.5))
c1 <- c0*(0.448+0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
R50 <- data.frame(t=t, m1=m1, c1=c1)
R50 <- R50 %>%
  mutate (r = m1/(fi*ra/rb)/(c1/ra-1))
Percentile<- rep(0.5, 365)
R50<- cbind(R50,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.6))
c1 <- c0*(0.448+0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
R60 <- data.frame(t=t, m1=m1, c1=c1)
R60 <- R60 %>%
  mutate (r = m1/(fi*ra/rb)/(c1/ra-1))
          Percentile<- rep(0.6, 365)
R60<- cbind(R60,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.7))
c1 <- c0*(0.448+0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
R70 <- data.frame(t=t, m1=m1, c1=c1)
R70 <- R70 %>%
  mutate (r = m1/(fi*ra/rb)/(c1/ra-1))
Percentile<- rep(0.7, 365)
R70<- cbind(R70,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.8))
c1 <- c0*(0.448+0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
R80 <- data.frame(t=t, m1=m1, c1=c1)
R80 <- R80 %>%
  mutate (r = m1/(fi*ra/rb)/(c1/ra-1))
Percentile<- rep(0.8, 365)
R80<- cbind(R80,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.9))
c1 <- c0*(0.448+0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
R90 <- data.frame(t=t, m1=m1, c1=c1)
R90 <- R90 %>%
  mutate (r = m1/(fi*ra/rb)/(c1/ra-1))
Percentile<- rep(0.9, 365)
R90<- cbind(R90,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.95))
c1 <- c0*(0.448+0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
R95 <- data.frame(t=t, m1=m1, c1=c1)
R95 <- R95 %>%
  mutate (r = m1/(fi*ra/rb)/(c1/ra-1))
Percentile<- rep(0.95, 365)
R95<- cbind(R95,Percentile)

data1<-rbind(R95, R90, R80, R70, R60, R50)

# Plot
R3 <- ggplot(data1, aes(x=t, y=r, color=Percentile)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  geom_line()+
  scale_y_log10(limits = c(0.05,100))+
  geom_hline(aes(yintercept = 1), linetype = "dashed", colour = "red")+
  theme(legend.position=c(0.9,0.2)) + labs(colour = "Exposue")+
  ylab("Population extinction risk")+xlab("Annual day")+
  scale_color_gradient(low="blue", high="red")

x11(8, 11)
grid.arrange(grid.arrange(R1, R2, ncol=2), R3, ncol=1)

# Sensitivity
library(sensitivity)
library(randtoolbox)
risk.fun <- function(X){
  ra <- 0.007; rb<-0.018; fi<-0.11
  c <- 0.1 * X[, 1]; m <- 0.14 * X[, 2]; md <- X[, 3]
  (m+md)/(fi*ra/rb)/(c/ra-1)
}

x <- delsa(model=risk.fun,
           par.ranges=list(c=c(0.1, 1),m=c(0.29, 1), md=c(0.01,0.347)),
           samples=c(10,10,10),method="grid")

# Summary of sensitivity indices of each parameter across parameter space
print(x)
library(ggplot2)
library(reshape2)
library(tidyr) # long_par
x11(8, 11)
plot(x)

# Parameter value vs. Output---
par <- data.frame(x$X0)
colnames(par) <- c("c", "m", "md")
long_par <- par %>% gather(par, val)
out <- x$y[1:3000]
df <- cbind(long_par, out)
df$par_f <- factor(df$par, levels=c("c", "m", "md"))

#
delsafirst <- data.frame(x$delsafirst)
colnames(delsafirst) <- c("c", "m", "md")
long_df <- delsafirst %>% gather(delsafirst, Sensitivity)
df <- cbind(df, long_df)

# Plot
S1 <- ggplot(df,aes(val, out, color=Sensitivity))+
  geom_point(size = 2)+
  xlab("Parameter value")+
  ylab("Population extinction risk")+
  ggtitle("")+
  facet_grid(. ~ par_f, scale = "free")+theme(legend.position="top")+
  scale_y_log10(limits = c(0.05,100))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  geom_hline(aes(yintercept = 1), linetype = "dashed", colour = "red")+
  scale_color_gradient(low="pink", high="red")
S2 <- ggplot(long_df, aes(delsafirst, Sensitivity))+ geom_boxplot()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  xlab("Parameter")+
  ylab("DELSA first-order sensitivity")

x11(8, 11)
grid.arrange(S1, S2, ncol=1)

