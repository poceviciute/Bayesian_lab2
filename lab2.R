#lab 2

# Question 1

#a)

temps<-read.table("TempLinkoping.dat",header=TRUE)
#response: temp, covariate: time
temp <- temps$temp
time <- temps$time
temps_lm <-lm(temp ~ time + I(time^2))

#hyperparameters
mu0<-c(-10,90,-80)
omega0<-diag(3)
inv_omega0 <- solve(omega0)
nu0<-5
sigma20<-9


#b)
library(mvtnorm)

sim_X <- rchisq(50,nu0)
sigma2 <- nu0*sigma20/sim_X
betas <- data.frame(beta0 = c(), beta1 =c(), beta2= c())

for(i in 1:length(sigma2)){
  tp <- rmvnorm(n = 1, mean=mu0, sigma=sigma2[i]*inv_omega0)
  betas <- rbind(betas,tp)
}
colnames(betas)<-c("beta0", "beta1", "beta2")

X <- as.matrix(data.frame(intercept=rep(1, length(time)), time, time^2))
y <- apply(as.matrix(betas),1,function(a){X%*%a})

#plot
# Creating 50 shades of gray
cl <- grey.colors(50)
#plotting the results
plot(x=time, y=y[,1], col=cl[1], type="l", ylim=c(-15, 25), ylab="Computed temperature", xlab="Time", main="Regression curves")
for(i in 2:ncol(y)){
 lines(time, y[,i], col=cl[i]) 
}

#c)
n<-60
beta_hat<-solve(t(X)%*%X)%*%t(X)%*%temp
mun <- solve(t(X)%*%X+omega0)%*%(t(X)%*%X%*%beta_hat+omega0%*%mu0)
omegan <-t(X)%*%X+omega0
inv_omegan <- solve(omegan)
nun <- nu0+nrow(temps)
sigma2n <- (nu0*sigma20+(t(temp)%*%temp+t(mu0)%*%omega0%*%mu0-t(mun)%*%omegan%*%mun))/nun

sim_X2 <- rchisq(n,nun)
sigma2_post <- nun*as.vector(sigma2n)/sim_X2

betas_post <- data.frame(beta0 = c(), beta1 =c(), beta2= c())

for(i in 1:length(sigma2_post)){
  tp <- rmvnorm(n = 1, mean=mun, sigma=sigma2_post[i]*inv_omegan)
  betas_post <- rbind(betas_post,tp)
}
colnames(betas_post)<-c("beta0", "beta1", "beta2")

X <- as.matrix(data.frame(intercept=rep(1, length(time)), time, time^2))
y_post <- apply(as.matrix(betas_post),1,function(a){X%*%a})
y_mean <- rowMeans(y_post)

perc = 0.05*n
lowertail<-c()
uppertail<-c()
for(i in 1:nrow(y_post)){
  lowertail[i] <- y_post[i,order(y_post[i,], decreasing = FALSE)[perc+1]]
  uppertail[i] <- y_post[i,order(y_post[i,], decreasing = FALSE)[n-perc-1]]
}

plot(x=time, y=temp, pch=20,ylab="Temperature", xlab="Time", main="Posterior mean")
lines(time, y_mean, col="red") 
lines(time, lowertail, col="blue")
lines(time, uppertail, col="blue")

#ASK IN LAB ABOUT WHY INTERVAL IS NARROW

#d)
x_tilde<-c()
for(i in 1:n){
 x_tilde[i]<- -betas_post[i,2]/(2*betas_post[i,3])
}

#time where post temp is maximal

#e)

# Question 2

#a)

#b)

#c)