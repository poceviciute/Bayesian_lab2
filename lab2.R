#lab 2

# Question 1

#a)

temps<-read.table("TempLinkoping.dat",header=TRUE)
#response: temp, covariate: time
temp <- temps$temp
time <- temps$time
temps_lm <-lm(temp ~ time + I(time^2))

#hyperparameters
mu0<-c(-10,90,-80) #based on lm
omega0<-diag(3)
inv_omega0 <- solve(omega0)
nu0<-5
sigma20<-9

n1 <- 50
#b)
library(mvtnorm)

sim_X <- rchisq(n1,nu0)
sigma2 <- nu0*sigma20/sim_X

# betas <- matrix(nrow=n1, ncol=3)
# for(i in 1:n1){
#   betas[i,] <- rmvnorm(n = 1, mean=mu0, sigma=sigma2[i]*inv_omega0)
# }

betas_mat1 <- t(apply(as.matrix(sigma2),1,function(x){rmvnorm(n=1,mean=mu0,sigma=x*inv_omega0)}))
betas <- as.data.frame(betas_mat1)
colnames(betas)<-c("beta0", "beta1", "beta2")

X <- as.matrix(data.frame(intercept=rep(1, length(time)), time, time^2))
y <- apply(betas_mat1,1,function(a){X%*%a})

#plotting the results

# Creating 50 shades of gray
cl <- grey.colors(n1)
plot(x=time, y=y[,1], col=cl[1], type="l", ylim=c(-15, 25), ylab="Computed temperature", xlab="Time", main="Regression curves")
for(i in 2:ncol(y)){
 lines(time, y[,i], col=cl[i]) 
}

#c)
n<-60
beta_hat<-solve(t(X)%*%X)%*%t(X)%*%temp #same as found with lm
mun <- solve(t(X)%*%X+omega0)%*%(t(X)%*%X%*%beta_hat+omega0%*%mu0)
omegan <-t(X)%*%X+omega0
inv_omegan <- solve(omegan)
nun <- nu0+nrow(temps)
sigma2n <- (nu0*sigma20+(t(temp)%*%temp+t(mu0)%*%omega0%*%mu0-t(mun)%*%omegan%*%mun))/nun

sim_X2 <- rchisq(n,nun)
sigma2_post <- nun*as.vector(sigma2n)/sim_X2

# betas_post <- matrix(nrow = n, ncol=3)
# for(i in 1:n){
#   betas_post[i,] <- rmvnorm(n = 1, mean=mun, sigma=sigma2_post[i]*inv_omegan)
# }

betas_mat<-t(apply(as.matrix(sigma2_post), 1, function(x){rmvnorm(n=1, mean=mun, sigma=x*inv_omegan)}))
betas_post <- as.data.frame(betas_mat)
colnames(betas_post)<-c("beta0", "beta1", "beta2")

y_post <- apply(betas_mat,1,function(a){X%*%a})
y_mean <- rowMeans(y_post)

perc = 0.05*n
# lowertail<-vector(length=nrow(y_post))
# uppertail<-vector(length=nrow(y_post))
# for(i in 1:nrow(y_post)){
#   lowertail[i] <- y_post[i,order(y_post[i,], decreasing = FALSE)[perc+1]]
#   uppertail[i] <- y_post[i,order(y_post[i,], decreasing = FALSE)[n-perc-1]]
# }

low <- apply(y_post,1,function(x){x[order(x, decreasing = FALSE)[perc+1]]})
upp <- apply(y_post,1,function(x){x[order(x, decreasing = FALSE)[n-perc-1]]})

plot(x=time, y=temp, pch=20, ylab="Temperature", xlab="Time", main="Posterior mean")
lines(time, y_mean, col="red") 
lines(time, low, col="blue")
lines(time, upp, col="blue")

#ASK IN LAB ABOUT WHY INTERVAL IS NARROW

#d)

# x_tilde<-vector(length=n)
# for(i in 1:n){
#  x_tilde[i]<- -betas_post$beta1[i]/(2*betas_post$beta2[i])
# }
x_tilde <- apply(betas_mat,1,function(x){-x[2]/(2*x[3])})
#time where post temp is maximal

x_tilde_mat <- data.frame(intercept=rep(1,n), x_tilde, x_tilde^2)


#e)

# Question 2

#a)
work<-read.table("WomenWork.dat",header=TRUE)
glmModel <- glm(Work ~ 0 + ., data = work, family = binomial)

#b)


#c)