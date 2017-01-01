## Loading required package: ggplot2, methods, gridExtra
require(ggplot2)
require(gridExtra)
library(data.table) # rbindlist
require(investr)

# 1230
x<-c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
y0 <- c(100, 100, 100, 100, 100, 100, 100, 100, 97, 94, 90)
y025 <- c(100, 97, 97, 94, 94, 94, 94, 88, 88, 88, 88)
y05 <- c(100, 94, 94, 90, 90, 90, 90, 90, 86, 86, 80)
y1 <- c(100, 97, 94, 88, 80, 78, 78, 70, 70, 68, 68)
y2 <- c(100, 90, 78, 70, 70, 70, 70, 62, 52, 42, 40)
y4 <- c(100, 92, 72, 60, 50, 40, 30, 30, 18, 0.1, 0.1)

df0 <- cbind(x, y0, 0); df025 <- cbind(x, y025, 250); df05 <- cbind(x, y05, 500);
df1 <- cbind(x, y1, 1000); df2 <- cbind(x, y2, 2000); df4 <- cbind(x, y4, 4000);
df0<-as.data.frame(df0);df025<-as.data.frame(df025);df05<-as.data.frame(df05)
df1<-as.data.frame(df1);df2<-as.data.frame(df2);df4<-as.data.frame(df4)
colnames(df0)<-colnames(df025)<-colnames(df05)<-colnames(df1)<-colnames(df2)<-colnames(df4)<-c("day", "sur", "dos")
Dat1.list <- list(df0,df025,df05,df1,df2,df4)
Dat1 <- rbindlist(Dat1.list)

exp.model0 <- lm(log(y0)~ 0+x, offset=rep(4.60517, length(x)))
exp.model025 <- lm(log(y025)~ 0+x, offset=rep(4.60517, length(x)))
exp.model05 <- lm(log(y05)~ 0+x, offset=rep(4.60517, length(x)))
exp.model1 <- lm(log(y1)~ 0+x, offset=rep(4.60517, length(x)))
exp.model2 <- lm(log(y2)~ 0+x, offset=rep(4.60517, length(x)))
exp.model4 <- lm(log(y4)~ 0+x, offset=rep(4.60517, length(x)))

p1 <-ggplot(Dat1, aes(day, sur))+ scale_y_continuous(lim=c(0, 100))+
  stat_function(fun=function(x)100*exp(-0.004816*x), geom="line", color="grey") +
  stat_function(fun=function(x)100*exp(-0.014419*x), geom="line", color="grey") +
  stat_function(fun=function(x)100*exp(-0.019719*x), geom="line", color="grey") +
  stat_function(fun=function(x)100*exp(-0.043743*x), geom="line", color="grey") +
  stat_function(fun=function(x)100*exp(-0.084599*x), geom="line", color="grey") +
  stat_function(fun=function(x)100*exp(-0.191161*x), geom="line", color="grey")+ geom_point(aes(shape = factor(dos)))+
  ylab("Survival (%)")+xlab(bquote('Day'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  labs(shape="Dose")


mor.x<-c(0, 250, 500, 1000, 2000, 4000)
mor.y<-c(0.004816, 0.014419, 0.019719, 0.043743, 0.084599, 0.191161)

Dat2 <- data.frame(mor.x, mor.y2)
mod2 <- lm(mor.y~ mor.x+0.004816)
summary(mod2)

day <- c(1,32,60,91,121,152,182,213,244,274,305,335,366,1,32,60,91,121,152,182,213,244,274,305,335,366,1,32,60,91,121,152,182,213,244,274,305,335,366)
rate_ratio<-c(0,0.1,0.63,0.8,0.98,0.99,1,1,0.34,0.45,0.56,0.28,0,0.29,0.29,0.29,0.54,0.79,0.68,0.57,0.78,1,0.83,0.67,0.29,0.29,0.1,0.03,0.03,0.48,0.93,0.97,1.00,0.60,0.19,0.32,0.45,0.25,0.1)
parameter<-c("l","l","l","l","l","l","l","l","l","l","l","l","l","m","m","m","m","m","m","m","m","m","m","m","m","m","c","c","c","c","c","c","c","c","c","c","c","c","c")
Data3 <- data.frame(parameter, day, rate_ratio)


## plot base + points
p2 <- ggplot(Dat2, aes(mor.x, mor.y)) + geom_point() + stat_smooth(method=lm, formula = y~x+0, color="blue")+
  ylab("Mortality rate (per day)")+xlab(bquote('Dose ('*mu*g~L^-1*')'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
p3 <- ggplot(Data3, aes(y=rate_ratio, x=day)) + facet_grid(. ~ parameter) + stat_smooth(method=lm, formula = y ~ sin(2 * pi * (x - 1)/365) + cos(2 * pi * (x - 1)/365), color="red") + geom_point()+
  ylab("Rate ratio")+xlab("Time (day)")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))

##Use grid.arrange(), from the gridExtra package
x11(8,11)
grid.arrange(p1,p2,p3, ncol=1)





## Creat dataframe
#dose<- c(0.1,1,2,3,4,5,6,7,8,9,10)
#acute_mortality_rate <- c(0.01474,0.11612,0.18873,0.23891,0.27584,0.30422,0.32675,0.34508,0.36031,0.37315,0.38414)
#Data1 <- data.frame(dose,acute_mortality_rate)
#Data1 (dose, missing_porb)
x <- c(0, 600, 1600, 3000, 4000, 6000)
y <- c(0, 0, 0, 22.6, 36.4, 51.6)
Data1 <- data.frame(x, y)
Data1
mod1 <- nls(y ~ 100*x^b/(a^b+x^b), data = Data1, start=list(a=5540,b=2.2))
summary(mod1)
RSS.p<-sum(residuals(mod1)^2)
TSS<-sum((y-mean(y))^2)
r.squared<-1-(RSS.p/TSS)
r.squared

# test 95% CI
plotFit(mod1,interval="prediction",ylim=c(0,60),pch=19,col.pred='light grey',shade=T,
        xlab=expression('Dose ('*mu*g~L^-1*')'), ylab= "Missing ratio (%)")

#1230 missing
# Dose = 1600
x<-c(0, 1, 2, 3, 4, 5)
y1<-c(1000,774,599,464,359,278)
y2<-c(1000,636,404,257,163,103)
y3<-c(1000,484,234,113,54,26)
exp.model1 <- lm(log(y1)~ x)
exp.model2 <- lm(log(y2)~ x)
exp.model3 <- lm(log(y3)~ x)

x<-c(0, 600, 1600, 3000, 4000, 6000)
y<-c(0, 0, 0, 0.256, 0.454, 0.730)