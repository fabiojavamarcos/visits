## ----setup, include=FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(dev = 'pdf')
knitr::opts_knit$set(global.par = TRUE)
par(mar= c(2, 2,.5,.5))


## ---- echo=TRUE, out.height = '50%', out.width = '50%', warning=FALSE-------------------------------------------------------
library(coda)
library(knitr)
library(kableExtra)
library(faraway)
library(mvtnorm)

suppressPackageStartupMessages(
  library(rstan,quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE))
data(prostate, package="faraway")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(8675309)

standx<- function(x) { 
  rangex<- range(x)
 (x - rangex[1]) / diff(rangex) 
}


## ---- echo=TRUE, out.height = '50%', out.width = '50%'----------------------------------------------------------------------
zprostate <- sapply(prostate[, 1:8], function(x) standx(x) - mean(standx(x)))
zprostate <- data.frame(lpsa = prostate$lpsa, zprostate)
colnames(zprostate) <- c('lpsa', paste0("z.", colnames(prostate[, -9])))
head(zprostate)


## ---- echo=TRUE, out.height = '50%', out.width = '50%'----------------------------------------------------------------------
summary(zprostate.lm <- lm(lpsa ~ ., data = zprostate))


## ---- echo=TRUE, out.height = '50%', out.width = '50%'----------------------------------------------------------------------
X <- model.matrix(zprostate.lm)
y <- zprostate.lm$model[,1]
n <- nrow(X)
p <- ncol(X)

XtX<- t(X) %*% X
bhat<- coef(zprostate.lm)
ybar<- mean(zprostate$lpsa)
nu0 <- ybar
yhat<- X %*%bhat
epshat<-y-yhat
sigmahat<-sum(epshat^2)/(n-p)


## ---- echo=TRUE, out.height = '50%', out.width = '50%'----------------------------------------------------------------------
sig02<- 0; ## prior scale squared
m0<- rep(0,p) ## prior beta mean
V<- (log(10)/1.96)^2
Sigma0inv<- diag(c(0,rep(1/V,p-1)))
nuhat<- nu0 + n
M <- 10000 ## MCMC iterations

beta_all <- list()
sigma2_all <- matrix(rep(NA, (M+1) * 3), ncol = 3)
sigma2_all[1, 1]<- sigmahat^2/4 ## initial sigma20
sigma2_all[1, 2]<- sigmahat^2  ## initial sigma20
sigma2_all[1, 3]<- 4 * sigmahat^2 ## initial sigma20
sXtX <- solve(XtX)


## ---- echo=TRUE, out.height = '50%', out.width = '50%'----------------------------------------------------------------------
for (s in 1:3){
  sigma2 <- as.numeric(sigma2_all[, s])
  beta <- matrix(NA,nrow=p,ncol=M+1)
  for (i in 1:M){
    W <- solve(sigma2[i] * Sigma0inv + XtX)%*%XtX 
    mhat <- (diag(p) - W) %*% m0 + W %*% bhat #
    Sighat <- W %*% sXtX 
    beta[, i+1] <- as.vector(rmvnorm(n=1,mean=mhat, sigma=sigma2[i]*Sighat))
    sigma2hat <- (nu0*sig02 + sum((y - X%*%beta[,i+1])^2))/nuhat
    sigma2[i+1] <- nuhat * sigma2hat / rchisq(n=1, df=nuhat)
  }
  beta_all[[s]] <- beta
  beta_all[[s]] <- beta
  sigma2_all[, s] <- sigma2
}
sigma2.1 <- sigma2_all[, 1]
sigma2.2 <- sigma2_all[, 2]
sigma2.3 <- sigma2_all[, 3]
beta1 <- beta_all[[1]]
beta2 <- beta_all[[2]]
beta3 <- beta_all[[3]]


## ---- echo=TRUE, eval=FALSE-------------------------------------------------------------------------------------------------
## par(mfrow=c(2, 3))
## for(j in 1:9){
##   plot(1:M, beta1[j,-1], type="l", xlab="", ylab="", main = colnames(zprostate)[j])
##   lines(1:M, beta2[j,-1], col = 2, lty = 2)
##   lines(1:M, beta3[j,-1], col = 4, lty = 2)
##   legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
##   legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
##          bty = 'n', lty = 1:3, col = c(1:2, 4))
## }; for(j in 1:9){
##   plot(density(na.omit(beta1[j,])), main = NA)
##   lines(density(na.omit(beta2[j,])), col = 2, lty = 2)
##   lines(density(na.omit(beta3[j,])), col = 4, lty = 4)
##   legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
##          bty = 'n', lty = 1:3, col = c(1:2, 4))
## }


## ---- echo=FALSE,  out.height = '80%', out.width = '80%', eval=TRUE---------------------------------------------------------
par(mfrow=c(2, 3),mar = rep(2, 4))
for(j in 1:3){
  plot(1:M, beta1[j,-1], type="l", xlab="", ylab="", main = colnames(zprostate)[j])
  lines(1:M, beta2[j,-1], col = 2, lty = 2)
  lines(1:M, beta3[j,-1], col = 4, lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}; for(j in 1:3){
  plot(density(na.omit(beta1[j,])), main = NA)
  lines(density(na.omit(beta2[j,])), col = 2, lty = 2)
  lines(density(na.omit(beta3[j,])), col = 4, lty = 4)
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}


## ---- echo=FALSE, out.height = '80%', out.width = '80%', eval=T-------------------------------------------------------------
par(mfrow=c(2, 3))
for(j in 4:6){
  plot(0:M, beta1[j,], type="l", xlab="", ylab="", main = colnames(zprostate)[j])
  lines(0:M, beta2[j,], col = 2, lty = 2)
  lines(0:M, beta3[j,], col = 4, lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}; for(j in 4:6){
  plot(density(na.omit(beta1[j,])), main = NA)
  lines(density(na.omit(beta2[j,])), col = 2, lty = 2)
  lines(density(na.omit(beta3[j,])), col = 4, lty = 4)
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}


## ---- echo=FALSE, out.height = '70%', out.width = '70%',eval=T--------------------------------------------------------------
par(mfrow=c(2, 3),mar = rep(2, 4))
for(j in 7:9){
  plot(0:M, beta1[j,], type="l", xlab="", ylab="", main = colnames(zprostate)[j])
  lines(0:M, beta2[j,], col = 2, lty = 2)
  lines(0:M, beta3[j,], col = 4, lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}; for(j in 7:9){
  plot(density(na.omit(beta1[j,])), main = NA)
  lines(density(na.omit(beta2[j,])), col = 2, lty = 2)
  lines(density(na.omit(beta3[j,])), col = 4, lty = 4)
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval=T---------------------------------------------------------------
stats <- c('mean', 'med', 'sd', 'x2.5%', 'x97.5%')
niceTables <- rbind(
  data.frame(grp = 'lm', stats = c('mean', 'sd'), 
             matrix(t(summary(zprostate.lm)$coefficients[, c('Estimate', 'Std. Error')]), nrow = 2)),
  data.frame(grp = 'Gibbs', stats, apply(
    cbind(beta1[, -1], beta2[, -1], beta3[, -1]),
          1, function(x) {
    c(mean(x), median(x), sd(x), quantile(x, c(.025, .975)))})))
colnames(niceTables) <- c('gr', 'stats', colnames(zprostate))
#head(niceTables)


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval=T---------------------------------------------------------------
kable(xtabs(lpsa ~ gr + stats, data = niceTables))
kable(xtabs(z.lcavol ~ gr + stats, data = niceTables))
kable(xtabs(z.lweight ~ gr + stats, data = niceTables))


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval=T---------------------------------------------------------------
kable(xtabs(z.age ~ gr + stats, data = niceTables))
kable(xtabs(z.lbph ~ gr + stats, data = niceTables))
kable(xtabs(z.svi ~ gr + stats, data = niceTables))


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval=T---------------------------------------------------------------
kable(xtabs(z.lcp ~ gr + stats, data = niceTables))
kable(xtabs(z.pgg45 ~ gr + stats, data = niceTables))
kable(xtabs(z.pgg45 ~ gr + stats, data = niceTables))


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval=T---------------------------------------------------------------
### writeLines(readLines("./Stan/zprostate.stan"))
# functions {
#   // nothing for this example
# }

# data{
#   int<lower=1> N; 
#   int<lower=1> p; 
#   matrix[N,p] X; 
#   vector[N] y;
#   vector[p] m0;
#   matrix[p,p] prec0;
#   real<lower=0> nu0;
#   real<lower=0> sig20;
# }
 
# transformed data{
#   // nothing for this example
#}


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval=T---------------------------------------------------------------
# parameters{
#   vector[p] beta;
#   real<lower=0> sigma2;
# }
 
# transformed parameters{
#   real<lower=0> sigma=sqrt(sigma2);
# }
 
# model{
#   sigma2 ~ scaled_inv_chi_square(nu0, sig20);
#   y ~ normal(X*beta, sigma);
#   beta ~ multi_normal_prec(m0, prec0);
# }

# generated quantities{
#   // nothing for this example
# }


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval= TRUE, cache=TRUE-----------------------------------------------
zprostate.stanc <- stanc(file="./Stan/zprostate.stan")
zprostate.stanmod <- stan_model(stanc_ret=zprostate.stanc)

zprostate.data <- list(
  N=dim(X)[1],
  p=dim(X)[2],
  X=X,
  y=zprostate.lm$model[,1],
  m0=rep(0, dim(X)[2]),
  prec0=diag(rep(1/V, ncol(X))),
  nu0=1/1000,
  sig20 = 1)


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval= TRUE-----------------------------------------------------------
set.seed(5551212)
binits1 <- rmvnorm(n=1,mean=coef(zprostate.lm), sigma=vcov(zprostate.lm)/4)
binits2 <- rmvnorm(n=1,mean=coef(zprostate.lm), sigma=vcov(zprostate.lm))
binits3 <- rmvnorm(n=1,mean=coef(zprostate.lm), sigma=vcov(zprostate.lm)*4)
attach(zprostate.data)
sigma2inits1 <- sum((y - X%*%binits1[1,])^2)/(N-p)
sigma2inits2 <- sum((y - X%*%binits2[1,])^2)/(N-p)
sigma2inits3 <- sum((y - X%*%binits3[1,])^2)/(N-p)
detach(zprostate.data)


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval= TRUE-----------------------------------------------------------
zprostate.init <- list(
  list(beta=binits1[1,], sigma2 = sigma2inits1),
  list(beta=binits2[1,], sigma2 = sigma2inits2),
  list(beta=binits3[1,], sigma2 = sigma2inits3)
)

zprostate.fit <- rstan::sampling(zprostate.stanmod,
                                  seed=8675309, warmup = 5000,
                                  data=zprostate.data,
                                  pars=c("beta","sigma2"),
                                  chains = 3, iter=10000,
                                  init = zprostate.init,
                                  refresh = 1000)


## ---- echo=F, out.height = '50%', out.width = '50%',eval=F------------------------------------------------------------------
## #load(file="./Stan/zprostate.fit.RData")


## ---- echo=TRUE, out.height = '50%', out.width = '50%',eval=T---------------------------------------------------------------
(zprostate.sum<- summary(zprostate.fit))
mcmc <- As.mcmc.list(zprostate.fit)

nchain(mcmc)
niter(mcmc)
nvar(mcmc)
varnames(mcmc)


## ---- echo=TRUE, eval=F-----------------------------------------------------------------------------------------------------
## for(j in 1:3){
##   plot(as.numeric(mcmc[[1]][, j]), type="l", xlab="",
##        ylab="", main = colnames(mcmc)[j])
##   lines( mcmc[[2]][, j], col = alpha(2, .5), lty = 2)
##   lines( mcmc[[3]][, j], col = alpha(4, .5), lty = 2)
##   legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
## }; for(j in 1:3){
##   plot(density(mcmc[[1]][, j]), main = NA)
##   lines(density(mcmc[[2]][, j]), col = 2, lty = 2)
##   lines(density(mcmc[[3]][, j]), col = 4, lty = 4)
##   legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
##          bty = 'n', lty = 1:3, col = c(1:2, 4))
## }


## ---- echo=FALSE, out.height = '80%', out.width = '80%', eval=T-------------------------------------------------------------
for(j in 1:3){
  plot(as.numeric(mcmc[[1]][, j]), type="l", xlab="", ylab="", main = names(zprostate.fit)[j])
  lines( mcmc[[2]][, j], col = alpha(2, .5), lty = 2)
  lines( mcmc[[3]][, j], col = alpha(4, .5), lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
}; for(j in 1:3){
  plot(density(mcmc[[1]][, j]), main = NA)
  lines(density(mcmc[[2]][, j]), col = 2, lty = 2)
  lines(density(mcmc[[3]][, j]), col = 4, lty = 4)
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}


## ---- echo=FALSE,  out.height = '80%', out.width = '80%', eval=T------------------------------------------------------------
for(j in 4:6){
  plot(as.numeric(mcmc[[1]][, j]), type="l", 
       xlab="", ylab="", main = names(zprostate.fit)[j])
  lines( mcmc[[2]][, j], col = alpha(2, .5), lty = 2)
  lines( mcmc[[3]][, j], col = alpha(4, .5), lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
}; for(j in 4:6){
  plot(density(mcmc[[1]][, j]), main = NA)
  lines(density(mcmc[[2]][, j]), col = 2, lty = 2)
  lines(density(mcmc[[3]][, j]), col = 4, lty = 4)
  legend('topright',  bty = 'n', lty = 1:3, col = c(1:2, 4),
         legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)))
}


## ---- echo=FALSE,  out.height = '80%', out.width = '80%', eval=T------------------------------------------------------------
for(j in 7:9){
  plot(as.numeric(mcmc[[1]][, j]), type="l", xlab="", ylab="",
       main = names(zprostate.fit)[j])
  lines( mcmc[[2]][, j], col = alpha(2, .5), lty = 2)
  lines( mcmc[[3]][, j], col = alpha(4, .5), lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
}; for(j in 7:9){
  plot(density(mcmc[[1]][, j]), main = NA)
  lines(density(mcmc[[2]][, j]), col = 2, lty = 2)
  lines(density(mcmc[[3]][, j]), col = 4, lty = 4)
  legend('topright', bty = 'n', lty = 1:3, col = c(1:2, 4),
         legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)), )
}


## ---- echo=FALSE,  out.height = '80%', out.width = '80%',eval=T-------------------------------------------------------------
par(mfrow = c(2,2))
for(j in 10:11){
  plot(as.numeric(mcmc[[1]][, j]), type="l", xlab="", ylab="", 
       main = names(zprostate.fit)[j])
  lines( mcmc[[2]][, j], col = alpha(2, .5), lty = 2)
  lines( mcmc[[3]][, j], col = alpha(4, .5), lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
}; for(j in 10:11){
  plot(density(mcmc[[1]][, j]), main = NA)
  lines(density(mcmc[[2]][, j]), col = 2, lty = 2)
  lines(density(mcmc[[3]][, j]), col = 4, lty = 4)
  legend('topright', bty = 'n', lty = 1:3, col = c(1:2, 4),
         legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)))
}


## ---- echo=TRUE, out.height = '70%', out.width = '70%',eval=T---------------------------------------------------------------
p2lo <- apply(t(cbind(beta1[2:9,-c(1:5001)], beta2[2:9,-c(1:5001)], 
                      beta3[2:9,-c(1:5001)])), 2, quantile, probs=c(0.025)) 
p2hi <- apply(t(cbind(beta1[2:9,-c(1:5001)], beta2[2:9,-c(1:5001)],
                      beta3[2:9,-c(1:5001)])), 2, quantile, probs=c(0.975))

p3lo <- apply(rbind(data.frame(mcmc[[1]]), data.frame(mcmc[[2]]),
                    data.frame(mcmc[[3]])), 2, quantile, probs=c(0.025))[2:9]
p3hi <- apply(rbind(data.frame(mcmc[[1]]), data.frame(mcmc[[2]]),
                    data.frame(mcmc[[3]])), 2, quantile, probs=c(0.975))[2:9]


## ---- echo=TRUE, eval=F-----------------------------------------------------------------------------------------------------
## par(mar=c(5,7,1,1)+.1, mfrow = c(1, 1))
## plot(p2lo, p2hi,xlim=c(-2.2,4.2),type="n",xlab="95% Credible Interval",
##      ylab="",axes="F",ylim=c(0,10))
## box(); axis(1)
## axis(2,at=seq(1,8),las=1,
##      labels=c("log can vol","log weight","age","log bph","svi", "log cap pen","gleason","pgg45"))
## 
## for (i in 1:8){
##   lines(y=c(i,i),x=c(p2lo[i],p2hi[i]),lty=2)
##   lines(y=c(i-.2,i-.2),x=c(p3lo[i],p3hi[i]),lty=3)
## }
## abline(v=0,lty=4)
## legend("topright",legend=c("Inform. Prior Gibbs", "Inform. Prior HMC"),lty=2:3,bty="n")


## ---- echo=F, out.height = '80%', out.width = '80%',eval=T------------------------------------------------------------------
par(mar=c(5,7,1,1)+.1, mfrow = c(1, 1))
plot(p2lo, p2hi,xlim=c(-2.2,4.2),type="n",xlab="95% Credible Interval", 
     ylab="",axes="F",ylim=c(0,10))
box(); axis(1)
axis(2,at=seq(1,8), las=1,
     labels=c("log can vol","log weight","age","log bph","svi", "log cap pen","gleason","pgg45"))

for (i in 1:8){ 
  lines(y=c(i,i),x=c(p2lo[i],p2hi[i]),lty=2)
  lines(y=c(i-.2,i-.2),x=c(p3lo[i],p3hi[i]),lty=3)  
}
abline(v=0,lty=4) 
legend("topright",legend=c("Inform. Prior Gibbs", "Inform. Prior HMC"), 
       lty=2:3,bty="n")

