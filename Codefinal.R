####### PROJET DE SERIES TEMPORELLES ########
#############################################

rm(list=ls())

## Importation packages : 
require(zoo)
require(tseries)
require(fUnitRoots)
require(stargazer)
require(urca)
require(FitAR)
require(caschrono)
require(tidyverse)
require(pracma)

# Importation des donnees :
###########################
data=read.csv("/Users/jeremymarck//Desktop/Projetst/data1.csv",sep=";",header=T)
colnames(data)=c('Dates','IPI','Codes')

# Mise en place des dates :
###########################
data_sources=as.character(data$Date)
ipi_source=data$IPI
start=1990+1/12
end=2019+2/12
dates_sequence=seq(from=start,to=end,by=1/12) # d??finit un s??quencage des dates
dates=as.yearmon(dates_sequence)
ipi=zoo(ipi_source,order.by = dates)
dev.off()

# Plot de la serie brute :
##########################
dev.off()
plot(ipi,main="Serie brute",ylab="IPI",xlab="Dates")
## Deux remarques : 
# - presence d'un trend baissier puis haussier
# - presence d'une saisonnalite
acf(ipi)
pacf(ipi)

## Recuperation du log-ipi : 
###########################
lipi=log(ipi)
plot(lipi,main="Serie log-transformee",ylab="Log-ipi",xlab="Dates")
acf(lipi,lag.max=40,main="ACF Log-IPI")
pacf(lipi,lag.max=40)
# Saisonnalite d'ordre 12

## Differenciation saisonniere :
################################
lipi12=diff(lipi,12)
plot(lipi12)
acf(lipi12)
pacf(lipi12) # Saisonnalite donc on differencie


## Differenciation d'ordre 1: car sinon on n'avait pas quelque chose de stationnaire
#############################
dlipi12=diff(lipi12,1)
plot(dlipi12,xlab='Dates')
acf(dlipi12) # lagmax de 13 pour le test adf


## Tests de stationnarisation :
###############################
test=ur.df(y=dlipi12,type="trend",lags=13)
summary(test)
kpss.test(x=dlipi12,null='Trend')
stargazer(a)# Apparement c'est bon 
# Ici l'hypothese nulle est la stationnarite et on accepte l'hypothse nulle
help(kpss.test)

pp.test(dlipi12)
stargazer(a)# Stationnaire
plot(dlipi12)
adf.test(dlipi12)
stargazer(a)# Stationnaire
# Pour les deux derniers tests, on rejette H0 qui est l'hypoth??se de racine unit?? 
# donc on rejette la non-stationnarit??. 
# Ok, on a un mod??le stationnaire 

## Choix d'un mod??le : 
######################
statio=dlipi12
par(mfrow=c(1,2))
acf(statio,lag.max=24,main='')
pacf(statio,lag.max=24,main='')

# On retient un SARIMA (2,1,2)(2,1,2)
#####################################
hessian=FALSE
modele1=arima(x=lipi,order=c(2,1,1),seasonal=list(order=c(1,1,2),period=12),method='ML')
t_stat(modele1)
modele1$aic
# AR(1) pas signif et sar(1) non plus --> a enlever, on commence par sar1
modele2=arima(x=lipi,order=c(1,1,1),seasonal=list(order=c(1,1,2),period=12),method='ML')
t_stat(modele2)
modele2$aic
# 
modele3=arima(x=lipi,order=c(1,1,1),seasonal=list(order=c(0,1,2),period=12),method='ML')
t_stat(modele3)
modele3$aic
# AR(1) a enlever
modele4=arima(x=lipi,order=c(0,1,1),seasonal=list(order=c(0,1,2),period=12),method='ML')
t_stat(modele4)
modele4$aic
#
modele5=arima(x=lipi,order=c(0,1,2),seasonal=list(order=c(0,1,2),period=12),method='ML')
t_stat(modele5)
modele5$aic
stargazer((a))

## Mod??le 5 retenu soit un SARIMA(0,1,2)(0,1,2)_{12}
####################################################

#### Blancheur des r??sidus : ###
################################

## Graphique et histogramme : 
modele=modele5
residus=modele$residuals
dev.off()
plot(residus,main="Residus du SARIMA retenu",xlab="Dates")
hist(residus,breaks=20,freq=F,main="Histogramme residus SARIMA retenu",col='lightblue',border='blue')
lines(density(residus),lwd=2)
lines(density(rnorm(1000000,mean(residus),sd(residus))),col='red')

## Graphique des p-valeurs : 
pval = c()
for(i in 1:50){
  pval[i]=Box.test(x=residus,lag=i,type="Ljung-Box")$p.value
}
vec = rep(0.05,50)
c=(rep(2,10))
plot(pval,xlab="Retards",col="red")
lines(vec,col="blue",lty=2)
legend(38,0.8,legend=c("p-valeurs", "y=0.05"),col=c("red", "blue"),lty=1:2,cex=0.7,text.font=4, bg='lightblue')
help(plot)
pval[11] # Point ambigu --> elle vaut 0.054, on accepte la blancheur des r??sidus :)


#### PREDICTION ####
####################

### Prediction de X_{T+1} et X_{T+2} :
x1=predict(modele,n.ahead = 2)$pred[1]
x2=predict(modele,n.ahead = 2)$pred[2]

### Variance estimee des residus du modele :
sigma2=var(residuals(modele))

### Parametre psi1 intervenant dans la variance de la prevision de X_{T+2} : 
psi1=modele$coef[4]
psi1

### Matrice de variance-covariance : 
v2=sigma2*(1+(1+psi1)^2)
cov=sigma2*(1+psi1)
sigma=matrix(c(sigma2,cov,cov,v2),nrow=2,ncol=2)

# R??gion de confiance : 
library(ellipse)
plot(ellipse(sigma,centre=c(x1,x2),type="l"),main="Region de confiance a 95% pour (X_{T+1},X_{T+2})",)

