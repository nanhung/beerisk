library(fitdistrplus)
library(ggplot2)
library(mc2d)

dat <- read.table("exposure.dat", head = T)
EU_non0 <- subset(dat, adj>0.1)
fitln_EU <- fitdist(EU_non0$adj,"lnorm",method="mle")

## fitdistrplus
meanlog_EU <- fitln_EU$est[1]
sdlog_EU <- fitln_EU$est[2]
summary(fitln_EU)
plot(fitln_EU)

## mc2d
ndvar(10001)
ndunc(1001)
mcstoc(rempiricalD,type="V",values=dat$adj)
Boot_EU <- bootdist(fitln_EU, bootmethod="param", niter=ndunc())
Mean_conso_EU <- mcdata(Boot_EU$estim$meanlog, type="U")
Sd_conso_EU <- mcdata(Boot_EU$estim$sdlog, type="U")
conso1_EU <- mcstoc(rlnorm, type="VU", meanlog= Mean_conso_EU, sdlog= Sd_conso_EU)
conso0 <- mcdata(0,type="V")
conso1_EU <- mcstoc(rlnorm, type="V", meanlog=meanlog_EU, sdlog=sdlog_EU)

## data manipulate
pzero_eu <- sum(dat$adj < 0.2)/length(dat$adj)
v_EU <- mcprobtree(c(pzero_eu,1-pzero_eu), list("0"=conso0,"1"=conso1_EU), type = "V")
v_EU<-as.numeric(v_EU)
dose<-data.frame(rep(deparse(substitute(v_EU)), length(v_EU)), v_EU)
colnames(dose) <- c("name","dose")

#plot exposure dose
ggplot(dose, aes(x=dose, fill=name))+
  geom_histogram(aes(y=..density..), alpha=0.8, colour="darkgrey", fill="grey", position="identity")+
  scale_x_log10(breaks=c(1,10,100,1000,10000))+
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))+
  ylab("Density \n")+xlab(bquote('Dose ('*mu*g~L^-1*')'))+theme(legend.position = "none")+
  geom_density(fill="white", alpha = 0.05)


## No used code
## transfer to md
md <- exp(0.0007767*v_EU-4.487)

## parameters
L<-2000; v<-5000; amax<-0.25; amin<-0.25; s<-0.75; m0<-0.1; b<-500 
ga<-0.007; gb<-0.018;c<-0.1; fi<-1/9

c_jan<-0.03; c_feb<-0.03; c_mar<-0.48; c_apr<-0.93; c_may<-0.97; c_jun<-1; c_jul<-0.6; c_aug<-0.19; c_sep<-0.32; c_oct<-0.45; c_nov<-0.25; c_dec<-0.1
m_jan<-0.29; m_feb<-0.29; m_mar<-0.54; m_apr<-0.79; m_may<-0.68; m_jun<-0.57; m_jul<-0.78; m_aug<-1; m_sep<-0.83; m_oct<-0.67; m_nov<-0.29; m_dec<-0.29

## risk
risk<-(m0)/(fi*ga/gb*(c/ga-1))
risk_jan<-(m0*m_jan+md)/(fi*ga/gb*(c*c_jan/ga-1))
risk_feb<-(m0*m_feb+md)/(fi*ga/gb*(c*c_feb/ga-1))
risk_mar<-(m0*m_mar+md)/(fi*ga/gb*(c*c_mar/ga-1))
risk_apr<-(m0*m_apr+md)/(fi*ga/gb*(c*c_apr/ga-1))
risk_may<-(m0*m_may+md)/(fi*ga/gb*(c*c_may/ga-1))
risk_jun<-(m0*m_jun+md)/(fi*ga/gb*(c*c_jun/ga-1))
risk_jul<-(m0*m_jul+md)/(fi*ga/gb*(c*c_jul/ga-1))
risk_aug<-(m0*m_aug+md)/(fi*ga/gb*(c*c_aug/ga-1))
risk_sep<-(m0*m_sep+md)/(fi*ga/gb*(c*c_sep/ga-1))
risk_oct<-(m0*m_oct+md)/(fi*ga/gb*(c*c_oct/ga-1))
risk_nov<-(m0*m_nov+md)/(fi*ga/gb*(c*c_nov/ga-1))
risk_dec<-(m0*m_dec+md)/(fi*ga/gb*(c*c_dec/ga-1))

## Confidence interval
median(risk_jan); median(risk_feb); median(risk_mar); median(risk_apr);median(risk_may); median(risk_jun); median(risk_jul); median(risk_aug); median(risk_sep); median(risk_oct); median(risk_nov); median(risk_dec)
quantile(risk_jan, c(0.025, 0.975)); quantile(risk_feb, c(0.025, 0.975)); quantile(risk_mar, c(0.025, 0.975)); quantile(risk_apr, c(0.025, 0.975)); quantile(risk_may, c(0.025, 0.975)); quantile(risk_jun, c(0.025, 0.975)); quantile(risk_jul, c(0.025, 0.975)); quantile(risk_aug, c(0.025, 0.975)); quantile(risk_sep, c(0.025, 0.975)); quantile(risk_oct, c(0.025, 0.975)); quantile(risk_nov, c(0.025, 0.975)); quantile(risk_dec, c(0.025, 0.975))

## Vector Combined (http://stackoverflow.com/questions/14620972/how-to-combine-two-vectors-into-a-data-frame)
A<-data.frame(rep(deparse(substitute(risk_jan)), length(risk_jan)), risk_mar)
colnames(A) <- c("Season","Risk")
B<-data.frame(rep(deparse(substitute(risk_feb)), length(risk_feb)), risk_apr)
colnames(B) <- c("Season","Risk")
C<-data.frame(rep(deparse(substitute(risk_mar)), length(risk_mar)), risk_may)
colnames(C) <- c("Season","Risk")
D<-data.frame(rep(deparse(substitute(risk_apr)), length(risk_apr)), risk_jun)
colnames(D) <- c("Season","Risk")
E<-data.frame(rep(deparse(substitute(risk_may)), length(risk_may)), risk_jul)
colnames(E) <- c("Season","Risk")
G<-data.frame(rep(deparse(substitute(risk_jun)), length(risk_jun)), risk_aug)
colnames(G) <- c("Season","Risk")
H<-data.frame(rep(deparse(substitute(risk_jul)), length(risk_jul)), risk_sep)
colnames(H) <- c("Season","Risk")
I<-data.frame(rep(deparse(substitute(risk_aug)), length(risk_aug)), risk_oct)
colnames(I) <- c("Season","Risk")
J<-data.frame(rep(deparse(substitute(risk_sep)), length(risk_sep)), risk_nov)
colnames(J) <- c("Season","Risk")
K<-data.frame(rep(deparse(substitute(risk_oct)), length(risk_oct)), risk_nov)
colnames(K) <- c("Season","Risk")
L<-data.frame(rep(deparse(substitute(risk_nov)), length(risk_nov)), risk_nov)
colnames(L) <- c("Season","Risk")
M<-data.frame(rep(deparse(substitute(risk_dec)), length(risk_dec)), risk_nov)
colnames(M) <- c("Season","Risk")

Data<-rbind(A,B,C,D,E,G,H,I,J,K,L,M)

##ggplot
ggplot(Data, aes(x=Season, y=Risk, fill=Season)) +
  geom_violin(binaxis='y', stackdir='center', stackratio=1.5, dotsize=2, binwidth=0.02) +
  scale_fill_brewer(palette="Blues", guide=FALSE) +
  stat_summary(fun.y=mean, geom="point", size=5, color="red") +
  ylab("Risk")+xlab("Month")+
  scale_y_log10(breaks=c(1,10,100))+
  scale_x_discrete(breaks=c("risk_mar", "risk_apr", "risk_may"),
                   labels=c("March", "April", "May"))+
  theme(axis.title.x = element_text(face="bold", colour="#990000", size=18),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(face="bold", colour="#990000", size=18),
        axis.text.y  = element_text(vjust=0.5, size=16))
  


