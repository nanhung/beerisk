if(!require(fitdistrplus)) {
  install.packages("fitdistrplus"); 
  require(fitdistrplus)
} 


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
ggplot(dose, aes(x=dose, fill=name))+
  geom_histogram(aes(y=..density..), alpha=0.8, colour="darkgrey", fill="grey", 
                 position="identity", bins = 50)+
  scale_x_log10(breaks=c(1,10,100,1000,10000), limits = c(0.01, 100000))+
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))+
  ylab("Density")+xlab(bquote('Dose ('*mu*g~L^-1*')'))+theme(legend.position = "none")+
  geom_density(fill="white", alpha = 0.05)


