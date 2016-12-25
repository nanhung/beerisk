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
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))+
  ylab("Probability")+xlab(bquote('Dose ('*mu*g~L^-1*')'))+theme(legend.position = "none")+
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

# 1223 Risk
EC50<-mcstoc(rnorm, type="VU", 5543.9879, 330.4317)
n <-mcstoc(rnorm, type="VU", 2.1720, 0.3314)
Imidacloprid <- mcstoc(rlnorm, type="VU", meanlog= Mean_conso_EU, sdlog= Sd_conso_EU, 
                    rtrunc=TRUE, linf=0, lsup=10000)
Missing_Ratio <- 100*Imidacloprid ^n/(EC50^n+Imidacloprid ^n)
plot(Missing_Ratio, xlab="Missing ratio", ylab="Cumulative probability")
missmodel <- mc(EC50, n, Imidacloprid, Missing_Ratio)
tor<-tornado(missmodel)
plot(tor)

#
Emax<-mcstoc(rnorm, type="VU", 3.476e-01, 5.594e-02)
EC50<-mcstoc(rnorm, type="VU", 3.941e+03, 7.193e+02)
n <-mcstoc(rnorm, type="VU", 2.003, 5.594e-02)
md <- Emax*Imidacloprid^n/(EC50^n+Imidacloprid^n)
m_d <- 1/0.1- 1/(0.1+Emax*Imidacloprid^n/(EC50^n+Imidacloprid^n))
plot(m_d, xlab="Decreasing lifespan", ylab="Cumulative probability")
mormodel <- mc(Emax, EC50, n, Imidacloprid, m_d)
tor<-tornado(mormodel)
plot(tor)

# Exceedance risk
library(reshape2)
library(plyr)
probs <- c(0.50, 0.75, 0.9, 0.95)
quant1 <- as.data.frame(t(apply(Missing_Ratio[,,], 1, quantile, probs = probs)))
quant2 <- as.data.frame(t(apply(m_d[,,], 1, quantile, probs = probs)))
cols <- c('green', 'lightyellow3', 'orange', 'red')
quant.melt1 <- suppressMessages(melt(quant1))
quant.melt2 <- suppressMessages(melt(quant2))
names(quant.melt1) <- c('quantile', 'x')
names(quant.melt2) <- c('quantile', 'x')
ecdf1 <- ddply(quant.melt1, c("quantile"), mutate, ecdf = ecdf(x)(unique(x))*length(x))
ecdf2 <- ddply(quant.melt2, c("quantile"), mutate, ecdf = ecdf(x)(unique(x))*length(x))
ecdf_1 <- ddply(ecdf1 , "quantile", mutate, 
                ecdf =scale(ecdf,center=min(ecdf),scale=diff(range(ecdf))))
ecdf_2 <- ddply(ecdf2 , "quantile", mutate, 
                ecdf =scale(ecdf,center=min(ecdf),scale=diff(range(ecdf))))
R1 <- ggplot(ecdf_1, aes(x,1-ecdf, colour = quantile)) + geom_step() +
  xlab('Missing ratio (%)') + ylab('Exxceedance risk') + geom_hline(yintercept = c(0, 1), linetype = "dashed", color = 'black') +
  theme(legend.position=c(0.85,0.8)) +  labs(colour = "Exposue quantile") + scale_colour_manual(values = cols)
R2 <- ggplot(ecdf_2, aes(x,1-ecdf, colour = quantile)) + geom_step() +
  xlab('Decreasing lifespan (day)') + ylab('Exxceedance risk') + geom_hline(yintercept = c(0, 1), linetype = "dashed", color = 'black') +
  theme(legend.position="none") +  labs(colour = "Exposue quantile") + scale_colour_manual(values = cols)

x11(8, 11)
grid.arrange(R1, R2, ncol=1)


# Find quantile
quantile(md, c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95))

# Parameters
l0 <- 2000; v<-5000; amax<-0.25; amin<-0.25; s<-0.75; m0<-0.1; b<-500 
ra <- 0.007; rb<-0.018; fi<-0.11

t <- 1:365

#
m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.5))
mp <- 1-(quantile(as.numeric(Missing_Ratio), c(0.5)))
Q <- data.frame(t=t, m1=m1, q=(-(amin-s-m1)+((amin-s-m1)^2+4*amin*m1)^0.5)/(2*amin))
L <- Q %>% 
  mutate (l1=l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365)))
R50 <- L %>%
  mutate (r = m1*(q+v)/(l1*q))
Percentile<- rep(0.5, 365)
R50<- cbind(R50,Percentile)


m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.6))
mp <- 1-(quantile(as.numeric(Missing_Ratio), c(0.6)))
Q <- data.frame(t=t, m1=m1, q=(-(amin-s-m1)+((amin-s-m1)^2+4*amin*m1)^0.5)/(2*amin))
L <- Q %>% 
  mutate (l1=l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365)))
R60 <- L %>%
  mutate (r = m1*(q+v)/(l1*q))
Percentile<- rep(0.6, 365)
R60<- cbind(R60,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.7))
mp <- 1-(quantile(as.numeric(Missing_Ratio), c(0.7)))
Q <- data.frame(t=t, m1=m1, q=(-(amin-s-m1)+((amin-s-m1)^2+4*amin*m1)^0.5)/(2*amin))
L <- Q %>% 
  mutate (l1=l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365)))
R70 <- L %>%
  mutate (r = m1*(q+v)/(l1*q))
Percentile<- rep(0.7, 365)
R70<- cbind(R70,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.8))
mp <- 1-(quantile(as.numeric(Missing_Ratio), c(0.7)))
Q <- data.frame(t=t, m1=m1, q=(-(amin-s-m1)+((amin-s-m1)^2+4*amin*m1)^0.5)/(2*amin))
L <- Q %>% 
  mutate (l1=l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365)))
R80 <- L %>%
  mutate (r = m1*(q+v)/(l1*q))
Percentile<- rep(0.8, 365)
R80<- cbind(R80,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.9))
mp <- 1-(quantile(as.numeric(Missing_Ratio), c(0.9)))
Q <- data.frame(t=t, m1=m1, q=(-(amin-s-m1)+((amin-s-m1)^2+4*amin*m1)^0.5)/(2*amin))
L <- Q %>% 
  mutate (l1=l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365)))
R90 <- L %>%
  mutate (r = m1*(q+v)/(l1*q))
Percentile<- rep(0.9, 365)
R90<- cbind(R90,Percentile)

m1 <-(m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365)))+quantile(as.numeric(md), c(0.95))
mp <- quantile(as.numeric(Missing_Ratio), c(0.95))
Q <- data.frame(t=t, m1=m1, q=(-(amin-s-m1)+((amin-s-m1)^2+4*amin*m1)^0.5)/(2*amin))
L <- Q %>% 
  mutate (l1=l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365)))
R95 <- L %>%
  mutate (r = m1*(q+v)/(l1*q))
Percentile<- rep(0.95, 365)
R95<- cbind(R95,Percentile)

data1<-rbind(R95, R90, R80, R70, R60, R50)

# Plot
ggplot(data1, aes(x=t, y=r, color=Percentile))+
  geom_line((aes(colour = Percentile)))+
  scale_y_log10(limits = c(0.05,100))+
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))+
  geom_hline(aes(yintercept = 1), linetype = "dashed", colour = "red")+
  theme(legend.position=c(0.8,0.8))+
  ylab("Risk")+xlab("Time (day)")
