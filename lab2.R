#lab 2

# Question 1

#a)

temps<-read.table("TempLinkoping.dat",header=TRUE)
#response: temp, covariate: time
temp <- temps$temp
time <- temps$time
temps_lm <-lm(temp ~ time + I(time^2))

#hyperparameters
mu0<-c(-4,90,-80) #based on lm
omega0<-diag(3)
inv_omega0 <- solve(omega0)
nu0<-5
sigma20<-9


#b)
library(mvtnorm)
n1 <- 50
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
n<-200
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
# x_tilde_mat <- data.frame(intercept=rep(1,n), x_tilde, x_tilde^2)

x_tilde <- apply(betas_mat,1,function(x){-x[2]/(2*x[3])})
hist(x_tilde)
#time where post temp is maximal

#e)


# Question 2

#a)
work<-read.table("WomenWork.dat",header=TRUE)
glmModel <- glm(Work ~ 0 + ., data = work, family = binomial)

#b)

y_work <- work$Work
X_work <- work[,-work$Work]
nparam <- ncol(X_work)

#Prior
mu <- as.vector(rep(0, nparam))
tau2<-100
sigma2 <- tau2*diag(nparam)

logpost_logistic <- function(betas, y, X, mu, sigma2){
  npara <- length(betas)
  X<-as.matrix(X)
  lin_pred <- X%*%betas
  loglik <- sum(lin_pred*y-log(1+exp(lin_pred)))
  if(abs(loglik)==Inf){
    loglik <- -20000
  }
  log_prior <- dmvnorm(betas, mean=mu, sigma2, log=TRUE)
  log_post<-loglik+log_prior
  return(log_post)
}

betas_init <- as.vector(rep(0,nparam))
opt_results <- optim(betas_init, logpost_logistic, gr=NULL, y_work, X_work, mu, sigma2, method=c("BFGS"),control=list(fnscale=-1), hessian=TRUE)
#opt_results2 <- optim(betas_init, logpost_logistic, gr=NULL, y_work, X_work, mu, sigma2, method=c("CG"),control=list(fnscale=-1), hessian=TRUE)
#CG gave smaller value of fn
beta_tilde <- opt_results$par
beta_hessian <- -1*opt_results$hessian 
inv_hessian <- solve(beta_hessian) # Posterior covariance matrix is -inv(Hessian)

sim_beta<-rmvnorm(n=100, mean=beta_tilde, sigma=inv_hessian)
colnames(sim_beta) <- colnames(X_work)
sim_NSmallChild <- sim_beta[,"NSmallChild"]
perc2 = 0.025*length(sim_NSmallChild)
lower <- sim_NSmallChild[order(sim_NSmallChild, decreasing = FALSE)[perc2+1]]
upper <- sim_NSmallChild[order(sim_NSmallChild, decreasing = FALSE)[length(sim_NSmallChild)-perc2-1]]

hist(sim_NSmallChild, freq = FALSE, breaks = 20, col = "grey70", xlab="beta for NSmallChild", main="Histogram of simulated NSmallChild beta")
lines(c(lower, upper), c(0.1,0.1), col="grey20", lwd=6)
#legend(x = -0.7, y=1.35, c("95% Equal tail"), col=c("grey20"), lwd = 3)


#c)

log_pred <- function(husband, edu, exper, age, smallchild, bigchild, n){
  sim_b<-rmvnorm(n=n, mean=beta_tilde, sigma=inv_hessian)
  x <- c(1, husband, edu, exper, (exper/10)^2, age, smallchild, bigchild)
  prob_pred <- exp(sim_b%*%x)/(1+exp(sim_b%*%x))
  u <- runif(n = n)
  y_pred <- ifelse(prob_pred>u, 1, 0)
  y_pred
}

woman <- log_pred(10, 8, 10, 40, 1, 1, 1000)
hist(woman, freq=FALSE, col = "grey70")
