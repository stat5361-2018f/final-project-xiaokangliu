library(spatstat)
library(grplasso)
library(MASS)
library(Matrix)

# read in the data
setwd("/Users/xiaokangliu/Desktop/statistical computing/hw/dslab-templates-master/")
load("bci.full7.rdata")
dat.tree <- dat.tree0 <- bci.full7
dat.cov <- read.table("bci-soil-20x20.txt", header = TRUE)

# add elevation and slope
ele <- as.data.frame(as.im(bei.extra$elev, W=owin(c(0,1000), c(0,500)),
                           dimyx = c(25, 50)))$value
slo <- as.data.frame(as.im(bei.extra$grad, W=owin(c(0,1000), c(0,500)),
                           dimyx = c(25, 50)))$value
dat.cov <- cbind(dat.cov,ele,slo)

# only the trees alive
dat.tree <- subset(dat.tree,status==c('A','P'))

# get covariates for all processes, M=1250
cov.mat0 <- cbind(dat.cov$x, dat.cov$y, scale(dat.cov[,-c(1,2)]))
cov.mat <- matrix(0, nrow = 1250, ncol = 15)
for (i in 1:50){
  cov.mat[c(((i-1)*25+1):(i*25)),] <- cov.mat0[rev(((i-1)*25+1):(i*25)),-c(1,2)]
}

p <- 15

############################################################
#####################     species      #### ################
############################################################

# select one species of tree
EachPPP <- function(spec.name, p, m){
  dat.spec <- subset(dat.tree,sp==spec.name)
  dat.spec <- dat.spec[,c(6,7,8)] # only keep the location information
  
  
  # construct ppp object
  x.spec <- dat.spec$gx
  y.spec <- dat.spec$gy
  ppp.spec <- ppp(x.spec, y.spec, c(0, 1000), c(0, 500))
  
  
  # quadrat count
  spec_ni <-  quadratcount.ppp(ppp.spec, nx = 50, ny = 25) # dim is 25 \times 50
  
  
  # plotting:
  #plot(ppp.spec, pch="+")
  #plot(spec_ni, add=TRUE, col="red", cex=1.5, lty=2)
  
  
  # get required input for Poisson GLM
  a <- 1000*500/50/25
  spec_ni_vec <- as.vector(spec_ni) # col by col, x: small->large
                                    # Note: y: large->small
  spec_ni_all <- spec_ni_vec+1 # obs + dummy points
  spec_delta <- (spec_ni_vec > 0) # delta_i
  spec_vi <- a/spec_ni_all 
  
  spec_ni_work <- c(spec_ni_vec, rep(1,1250))
  spec_delta_work <- c(spec_delta, rep(0, 1250))
  spec_vi_work <- rep(spec_vi, 2)
  
  cov.mat_work <- rbind(cov.mat, cov.mat)
  
  # fit Poisson GLM to get \tilde \beta^{EE}
  fit1 <- glm(spec_delta_work~cov.mat_work+offset(log(spec_vi_work)),
              weights = spec_ni_work,
              family = poisson)
  spec_beta_EE <- fit1$coefficients
  
  
  # get w(s) with m=10
  spec_lambda <- exp(cbind(rep(1,1250),cov.mat)%*%spec_beta_EE)
  spec_lambda_mat <- matrix(spec_lambda, nrow = 25, ncol = 50, byrow = FALSE)
  spec_lambda_im <- im(spec_lambda_mat, 
                         xcol=seq(0,1000,length=50), yrow=seq(0,500,length=25))
  #m <- 10
  spec_K <- Kinhom(ppp.spec, lambda = spec_lambda_im, r = c(0, m))$border[2]
  spec_f <- spec_K-pi*m^2
  
  spec_w <- 1/(1+spec_lambda*spec_f)
  
  # fit Poisson GLM to get \tilde \beta^{WEE}
  w <- rep(spec_w,2)
  fit2 <- glm(spec_delta_work~cov.mat_work+offset(log(spec_vi_work)),
              weights = spec_ni_work*w,
              family = poisson)
  spec_beta_WEE <- fit2$coefficients
  
  
  # set spec_beta_WEE as initial value
  init <- spec_beta_WEE
  
  A2 <- matrix(0, nrow = p+1, ncol = p+1)
  muls <- exp(cbind(rep(1,1250),cov.mat)%*%init)*a*spec_w
  Y0 <- rep(0, p+1)
  for (i in 1:1250){
    A2 <- A2+c(1,cov.mat[i,])%*%t(c(1,cov.mat[i,]))*muls[i]
    Y0 <- Y0+spec_ni_vec[i]*spec_w[i]*c(1,cov.mat[i,])-
      c(1,cov.mat[i,])*muls[i]
  }
  X <- A <- chol(A2)
  Y <- t(ginv(A))%*%Y0 + A%*%init
  
  # keep the first initial value in all the iterations
  # in each iteration, only updates w(s) and use new beta as initial values
  
  #Ynew <- Y-X[,1]*init[1]
  #Xnew <- X[,-1]
  #beta_init <- init[-1]
  
  Ynew <- Y
  Xnew <- X
  beta_init <- init
  
  return(list(Y=Ynew, X=Xnew, Beta_init=beta_init, Beta_0=init[1],
              ppp.spec=ppp.spec))
}


############################################################
#####################     genus       ######################
############################################################

# for genus: Eugenia
fit.spec1 <- EachPPP('eugeco', p = 15, m = 10)
fit.spec2 <- EachPPP('eugega', p = 15, m = 10)
fit.spec3 <- EachPPP('eugene', p = 15, m = 10)
fit.spec4 <- EachPPP('eugeoe', p = 15, m = 10)

# initial beta
grl_init <- c(fit.spec1$Beta_init,fit.spec2$Beta_init,fit.spec3$Beta_init,
              fit.spec4$Beta_init)
grl_X <- bdiag(fit.spec1$X,fit.spec2$X,fit.spec3$X,fit.spec4$X)
grl_Y <- c(fit.spec1$Y,fit.spec2$Y,fit.spec3$Y,fit.spec4$Y)

#index <- c(NA, rep(c(1:15),4))
index <- c(NA,rep(c(NA,1:15),4))
grl_init <- c(0, grl_init)
grl_X <- as.matrix(cbind(rep(1,dim(grl_X)[1]),grl_X))

lambda = exp(seq(log(20), log(0.005), length = 50))
fitlasso <- grplasso(x = grl_X, y = grl_Y, index = index, 
         coef.init = grl_init, model = LinReg(), lambda = lambda,
         control = grpl.control(trace = 0))
grp.est <- fitlasso$coefficients[-1,]

grps <- list()
#for (i in 1:p){
#  grps[[i]] <- grp.est[(1:4-1)*15+i,]
#}
for (i in 2:(p+1)){
  grps[[i-1]] <- grp.est[(1:4-1)*16+i,]
}
g1 <- colSums(grps[[1]]^2)
plot(g1)

par(mfrow = c(3,3))
for (i in 1:9){
  g1 <- colSums(grps[[i]]^2)
  plot(g1, ylab = "magnitude", ylim = c(-0.1, 0.7), type = "l")
}

par(mfrow = c(2,3))
for (i in 10:15){
  g1 <- colSums(grps[[i]]^2)
  plot(g1, ylab = "magnitude", ylim = c(-0.1, 0.7), type = "l")
}


# for genus: Inga
fit.spec12 <- EachPPP('ingas1', p = 15, m = 10)
fit.spec22 <- EachPPP('ingago', p = 15, m = 10)
fit.spec32 <- EachPPP('ingama', p = 15, m = 10)
fit.spec42 <- EachPPP('ingaqu', p = 15, m = 10)
fit.spec52 <- EachPPP('ingasa', p = 15, m = 10)
fit.spec62 <- EachPPP('ingath', p = 15, m = 10)
fit.spec72 <- EachPPP('ingaum', p = 15, m = 10)


# initial beta
grl_init2 <- c(fit.spec12$Beta_init,fit.spec22$Beta_init,fit.spec32$Beta_init,
              fit.spec42$Beta_init,fit.spec52$Beta_init,fit.spec62$Beta_init,
              fit.spec72$Beta_init)
grl_X2 <- bdiag(fit.spec12$X,fit.spec22$X,fit.spec32$X,fit.spec42$X,
                fit.spec52$X,fit.spec62$X,fit.spec72$X)
grl_Y2 <- c(fit.spec12$Y,fit.spec22$Y,fit.spec32$Y,fit.spec42$Y,
            fit.spec52$Y,fit.spec62$Y,fit.spec72$Y)

#index2 <- c(NA, rep(c(1:15),7))
index2 <- c(NA, rep(c(NA,1:15),7))
grl_init2 <- c(0, grl_init2)
grl_X2 <- as.matrix(cbind(rep(1,dim(grl_X2)[1]),grl_X2))

lambda = exp(seq(log(20), log(0.005), length = 50))
fitlasso2 <- grplasso(x = grl_X2, y = grl_Y2, index = index2, 
                     coef.init = grl_init2, model = LinReg(), lambda = lambda)
grp.est2 <- fitlasso2$coefficients[-1,]

grps2 <- list()

for (i in 2:(p+1)){
  grps2[[i-1]] <- grp.est2[(1:7-1)*16+i,]
}

par(mfrow = c(4,4))
for (i in 1:15){
  g1 <- colSums(grps2[[i]]^2)
  plot(g1, ylab = "magnitude", ylim = c(-0.1, 7), type = "l")
}





# for genus: Ocotea
fit.spec13 <- EachPPP('ocotce', p = 15, m = 10)
fit.spec23 <- EachPPP('ocotob', p = 15, m = 10)
fit.spec33 <- EachPPP('ocotpu', p = 15, m = 10)
fit.spec43 <- EachPPP('ocotwh', p = 15, m = 10)

# initial beta
grl_init3 <- c(fit.spec13$Beta_init,fit.spec23$Beta_init,fit.spec33$Beta_init,
              fit.spec43$Beta_init)
grl_X3 <- bdiag(fit.spec13$X,fit.spec23$X,fit.spec33$X,fit.spec43$X)
grl_Y3 <- c(fit.spec13$Y,fit.spec23$Y,fit.spec33$Y,fit.spec43$Y)

#index <- c(NA, rep(c(1:15),4))
index <- c(NA,rep(c(NA,1:15),4))
grl_init3 <- c(0, grl_init3)
grl_X3 <- as.matrix(cbind(rep(1,dim(grl_X3)[1]),grl_X3))

lambda = exp(seq(log(20), log(0.005), length = 50))
fitlasso3 <- grplasso(x = grl_X3, y = grl_Y3, index = index, 
                     coef.init = grl_init3, model = LinReg(), lambda = lambda,
                     control = grpl.control(trace = 0))
grp.est3 <- fitlasso3$coefficients[-1,]

grps3 <- list()
#for (i in 1:p){
#  grps[[i]] <- grp.est[(1:4-1)*15+i,]
#}
for (i in 2:(p+1)){
  grps3[[i-1]] <- grp.est3[(1:4-1)*16+i,]
}
g1 <- colSums(grps3[[1]]^2)
plot(g1)

par(mfrow = c(3,3))
for (i in 1:9){
  g1 <- colSums(grps3[[i]]^2)
  plot(g1, ylab = "magnitude", ylim = c(-0.1, 2), type = "l")
}

par(mfrow = c(2,3))
for (i in 10:15){
  g1 <- colSums(grps3[[i]]^2)
  plot(g1, ylab = "magnitude", ylim = c(-0.1, 2), type = "l")
}


