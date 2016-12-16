## Loading required package: ggplot2, methods, gridExtra
require(ggplot2)
require(gridExtra)

## Creat dataframe
#dose<- c(0.1,1,2,3,4,5,6,7,8,9,10)
#acute_mortality_rate <- c(0.01474,0.11612,0.18873,0.23891,0.27584,0.30422,0.32675,0.34508,0.36031,0.37315,0.38414)
#Data1 <- data.frame(dose,acute_mortality_rate)
#Data1
dose <- c(0, 600, 1600, 3000, 4000, 6000)
missing_porb <- c(0, 0, 0, 22.6, 36.4, 51.6)
Data1 <- data.frame(dose, missing_porb)
Data1

dose <- c(0,1,10,100,1000,1500,2000,2500,3000,3500,4000)
mortality_rate <- c(0,0.00002,0.00023,0.00233,0.02614,0.04231,0.06162,0.08557,0.11712,0.16348,0.24169)
Data2 <- data.frame(dose,mortality_rate)
Data2

day <- c(1,32,60,91,121,152,182,213,244,274,305,335,366,1,32,60,91,121,152,182,213,244,274,305,335,366,1,32,60,91,121,152,182,213,244,274,305,335,366)
rate_ratio<-c(0,0.1,0.63,0.8,0.98,0.99,1,1,0.34,0.45,0.56,0.28,0,0.29,0.29,0.29,0.54,0.79,0.68,0.57,0.78,1,0.83,0.67,0.29,0.29,0.1,0.03,0.03,0.48,0.93,0.97,1.00,0.60,0.19,0.32,0.45,0.25,0.1)
parameter<-c("l","l","l","l","l","l","l","l","l","l","l","l","l","m","m","m","m","m","m","m","m","m","m","m","m","m","c","c","c","c","c","c","c","c","c","c","c","c","c")
Data3 <- data.frame(parameter, day, rate_ratio)
Data3

#
mod <- nls(y ~ a*(x^c)/((b^c)+(x^c)), data = Data1, start=list(a=50,b=3000,c=1), trace=T, algorithm = "port")

## plot base + points
p1 <- ggplot(Data1, aes(dose, missing_porb)) + geom_point()
p2 <- ggplot(Data2, aes(dose, mortality_rate)) + geom_point()
p3 <- ggplot(Data3, aes(y=rate_ratio, x=day)) + facet_grid(. ~ parameter)

## Automatically create fitted function
p11<-p1 + stat_smooth(method = "nls", formula = y ~ a*(x^c)/((b^c)+(x^c)), data = Data1, method.args=list(start= c(a=60,b=3000,c=1)), se=FALSE)+
  ylab("Missing bee (%)")+xlab(bquote('Dose ('*mu*g~L^-1*')'))
p22<-p2 + stat_smooth(method = "nls", formula = y ~ a*(x^c)/((b^c)+(x^c)), data = Data2, method.args=list(start= c(a=0.2,b=3000,c=1)), se=FALSE)+
  ylab("Mortality rate (per day)")+xlab(bquote('Dose ('*mu*g~L^-1*')'))
p33<-p3 + stat_smooth(method=lm, formula = y ~ sin(2 * pi * (x - 1)/365) + cos(2 * pi * (x - 1)/365), color="red") + geom_point()+
  ylab("Rate ratio")+xlab("Time (day)")

##Use grid.arrange(), from the gridExtra package
grid.arrange(p11, p22,p33, ncol=1)
