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

#
x5 <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
y5 <- c(100, 92, 72, 60, 50, 40, 30, 30, 18)
exp.model5 <- lm(log(y5)~ 0+x5, offset=rep(4.60517, length(x5)))
#

exp.model0 <- lm(log(y0)~ 0+x, offset=rep(4.60517, length(x)))
exp.model025 <- lm(log(y025)~ 0+x, offset=rep(4.60517, length(x)))
exp.model05 <- lm(log(y05)~ 0+x, offset=rep(4.60517, length(x)))
exp.model1 <- lm(log(y1)~ 0+x, offset=rep(4.60517, length(x)))
exp.model2 <- lm(log(y2)~ 0+x, offset=rep(4.60517, length(x)))
exp.model4 <- lm(log(y4)~ 0+x, offset=rep(4.60517, length(x)))

p1 <-ggplot(Dat1, aes(day, sur))+ geom_point(aes(color = factor(dos)))+ scale_y_continuous(lim=c(0, 100))+
  stat_function(fun=function(x)100*exp(-0.004816*x), geom="line", color="gray10") +
  stat_function(fun=function(x)100*exp(-0.014419*x), geom="line", color="gray10") +
  stat_function(fun=function(x)100*exp(-0.019719*x), geom="line", color="gray10") +
  stat_function(fun=function(x)100*exp(-0.043743*x), geom="line", color="gray10") +
  stat_function(fun=function(x)100*exp(-0.084599*x), geom="line", color="gray10") +
  stat_function(fun=function(x)100*exp(-0.191161*x), geom="line", color="gray10")+ 
  ylab("%Survival")+xlab(bquote('Day'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  labs(color="Dose")+theme_bw()

#
P_y0 <- exp(predict(lm(log(y0)~ 0+x, offset=rep(4.60517, length(x)))))
P_y025 <- exp(predict(lm(log(y025)~ 0+x, offset=rep(4.60517, length(x)))))
P_y05 <- exp(predict(lm(log(y05)~ 0+x, offset=rep(4.60517, length(x)))))
P_y1 <- exp(predict(lm(log(y1)~ 0+x, offset=rep(4.60517, length(x)))))
P_y2 <- exp(predict(lm(log(y2)~ 0+x, offset=rep(4.60517, length(x)))))
P_y4 <- exp(predict(lm(log(y5)~ 0+x5, offset=rep(4.60517, length(x5)))))

P2_df <- cbind(c(y0,y025,y05,y1,y2,y5) ,c(P_y0,P_y025,P_y05,P_y1,P_y2,P_y4))
P2_df <- as.data.frame(P2_df)
dose <- c(rep("0",11),rep("250",11),rep("500",11),rep("1000",11),rep("2000",11),rep("4000",9))
rownames(P2_df) <- c(1:64)
P2_df <- cbind(dose, P2_df)
colnames(P2_df) <- c("Dose", "Observed", "Predicted")
P2_df[, 1] <- factor(P2_df[, 1], levels = c("0","250","500","1000","2000","4000"))

## plot base + points
p2<- ggplot(P2_df, aes(Predicted, Observed)) + geom_point(aes(color = Dose))+theme_bw()+
  scale_x_continuous(lim=c(0, 100))+scale_y_continuous(lim=c(0, 100)) + 
  geom_abline(intercept = 0, cex = 1) +
  theme(legend.position = "none", axis.text=element_text(size=12), axis.title=element_text(size=12))+
  xlab("Predicted %survial")+ylab('Observed %survial')


#
h_lif <- c(log(2)/(-(coef(summary(exp.model4))[, 1])),
           log(2)/(-(coef(summary(exp.model2))[, 1])),
           log(2)/(-(coef(summary(exp.model1))[, 1])),
           log(2)/(-(coef(summary(exp.model05))[, 1])))
h_dos <- c(4000,2000,1000,500)
h_df <- data.frame(h_lif,h_dos)
p3 <- ggplot(h_df, aes(h_lif, h_dos)) + geom_point() + scale_y_continuous(lim=c(0, 4000))+
  ylab(bquote('ED50 ('*mu*g~L^-1*')'))+xlab(bquote('Day'))+ 
  stat_function(fun=function(x)4457.9-1164*log(x), geom="line", color="gray10", cex = 1)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  labs(color="Dose")+theme_bw()



mor.x<-c(0, 250, 500, 1000, 2000, 4000)
mor.y<-c(0.004816, 0.014419, 0.019719, 0.043743, 0.084599, 0.191161)
mor.y.sd<-c(coef(summary(exp.model0))[, 2],
             coef(summary(exp.model025))[, 2],
             coef(summary(exp.model05))[, 2],
             coef(summary(exp.model1))[, 2],
             coef(summary(exp.model2))[, 2],
             coef(summary(exp.model5))[, 2])

Dat3 <- data.frame(mor.x, mor.y, mor.y.sd)
mod3 <- lm(mor.y~ mor.x)
summary(mod3)

## plot base + points
p4 <- ggplot(Dat3, aes(mor.x, mor.y)) + geom_point() + stat_smooth(method=lm, formula = y~x, color="gray10")+
  ylab("Mortality rate (per day)")+xlab(bquote('Dose ('*mu*g~L^-1*')'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  geom_errorbar(aes(x = mor.x, ymin = mor.y-mor.y.sd, ymax =  mor.y+mor.y.sd, width = .5))+
  theme_bw()

##Use grid.arrange(), from the gridExtra package
x11(8.5,11)
grid.arrange(p1, p2, arrangeGrob(p3, p4,ncol=2), ncol=1)







# Seasonality
day <- c(1,32,60,91,121,152,182,213,244,274,305,335,366,1,32,60,91,121,152,182,213,244,274,305,335,366,1,32,60,91,121,152,182,213,244,274,305,335,366)
rate_ratio<-c(0,0.1,0.63,0.8,0.98,0.99,1,1,0.34,0.45,0.56,0.28,0,0.29,0.29,0.29,0.54,0.79,0.68,0.57,0.78,1,0.83,0.67,0.29,0.29,0.1,0.03,0.03,0.48,0.93,0.97,1.00,0.60,0.19,0.32,0.45,0.25,0.1)
parameter<-c("l","l","l","l","l","l","l","l","l","l","l","l","l","m","m","m","m","m","m","m","m","m","m","m","m","m","c","c","c","c","c","c","c","c","c","c","c","c","c")
Data4 <- data.frame(parameter, day, rate_ratio)

#
x11(11,8.5)
ggplot(Data4, aes(y=rate_ratio, x=day)) + facet_grid(. ~ parameter) + 
  stat_smooth(method=lm, formula = y ~ sin(2 * pi * (x - 1)/365) + cos(2 * pi * (x - 1)/365), color="red") + geom_point()+
  ylab("Rate ratio")+xlab("Time (day)")+ 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
