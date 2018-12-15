library(spatstat)
library(grplasso)
library(MASS)
library(Matrix)


EachPPPsimu <- function(spec.name, p, m){
  # quadrat count
  ppp.spec <- spec.name
  spec_ni <-  quadratcount.ppp(ppp.spec, nx = 50, ny = 25) # dim is 25 \times 50
  
  
  # plotting:
  #plot(ppp.spec, pch="+")
  #plot(spec_ni, add=TRUE, col="red", cex=1.5, lty=2)
  
  # get required input for Poisson GLM
  a <- 1000*500/50/25
  spec_ni_vec <- as.vector(spec_ni) # row by row, row is for x start from 0-20
  # Note: y is in a reverse order
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
  # in each iteration, only updates W(s) and use new beta as initial values
  
  #Ynew <- Y-X[,1]*init[1]
  #Xnew <- X[,-1]
  #beta_init <- init[-1]
  
  Ynew <- Y
  Xnew <- X
  beta_init <- init
  
  
  return(list(Y=Ynew, X=Xnew, Beta_init=beta_init, Beta_0=init[1],
              ppp.spec=ppp.spec))
}


set.seed(123)

# set working path
setwd("/Users/xiaokangliu/Desktop/statistical computing/hw/dslab-templates-master/")
dat.cov <- read.table("bci-soil-20x20.txt", header = TRUE)

# elevation, slope and 13 irrelavent covariates
dat.noise <- matrix(rnorm(13*1250, 0, 1), nrow = 1250, ncol = 13)

# add elevation and slope
ele <- as.data.frame(as.im(bei.extra$elev, W=owin(c(0,1000), c(0,500)),
                           dimyx = c(25, 50)))$value
slo <- as.data.frame(as.im(bei.extra$grad, W=owin(c(0,1000), c(0,500)),
                           dimyx = c(25, 50)))$value
dat.comb <- cbind(dat.noise,ele,slo)

# get covariates for all processes, M=1250
cov.mat0 <- cbind(dat.cov$x, dat.cov$y, scale(dat.comb))
cov.mat <- matrix(0, nrow = 1250, ncol = 15)
for (i in 1:50){
  cov.mat[c(((i-1)*25+1):(i*25)),] <- cov.mat0[rev(((i-1)*25+1):(i*25)),-c(1,2)]
}

p <- 15

# add collinearity into the design matrix
col.mat <- matrix(0, nrow = p, ncol = p)
for (i in 1:p){
  for (j in 1:p){
    col.mat[i,j] <- col.mat[j,i] <- 0.5^abs(i-j)
  }
}
col.mat[14,15] <- col.mat[15,14] <- 0

spec.col <- eigen(col.mat)
hfcol.mat <- spec.col$vectors%*%diag(sqrt(spec.col$values))%*%t(spec.col$vectors)

cov.mat <- cov.mat%*%hfcol.mat


############################################################
##################     generate ppp      ###################
############################################################

nsim <- 200
k.par <- 5e-5 # 1e-4 or 5e-4
w.par <- 30
q.par <- 5 # number of point process

# generate a pixel image for mu
beta.sim <- matrix(0, nrow = q.par, ncol = p)

# sim1
beta.sim[,15] <- c(4,4.1,4.9,4.5,4.7)-3 
beta.sim[,13] <- c(2.3,2.2,2.2,2.6,2.51)-1

# sim2 
beta.sim[,14] <- c(4,4.1,4.9,4.5,4.7)-3 
beta.sim[,13] <- c(2.3,2.2,2.2,2.6,2.51)-1

# sim3
beta.sim[,1] <- c(4,4.1,4.9,4.5,4.7)-3 
beta.sim[,10] <- c(2.3,2.2,2.2,2.6,2.51)-1

# sim4
beta.sim[,14] <- c(4,4.1,4.9,4.5,4.7)-3 

# sim5
beta.sim[,15] <- c(4,4.1,4.9,4.5,4.7)-3 



res <- matrix(0, nrow = nsim, ncol = 15)
for (l in 1:nsim){
  beta.10 <- 0.2
  mu_par01 <- exp(cov.mat%*%beta.sim[1,] + beta.10)
  mu_par01_mat <- matrix(mu_par01, nrow = 25, ncol = 50, byrow = FALSE)
  mu_par01_im <- im(mu_par01_mat, 
                    xcol=seq(0,1000,length=50), yrow=seq(0,500,length=25))
  #plot(mu_par01_im)
  spe1.sim <- rThomas(kappa = k.par, scale = w.par, mu = mu_par01_im, 
                      win = owin(c(0,1000),c(0,500)), nsim = 1)
  #spe1.sim
  #plot(spe1.sim)
  
  
  beta.20 <- 0.2
  mu_par02 <- exp(cov.mat%*%beta.sim[2,] + beta.20)
  mu_par02_mat <- matrix(mu_par02, nrow = 25, ncol = 50, byrow = FALSE)
  mu_par02_im <- im(mu_par02_mat, 
                    xcol=seq(0,1000,length=50), yrow=seq(0,500,length=25))
  #plot(mu_par02_im)
  spe2.sim <- rThomas(kappa = k.par, scale = w.par, mu = mu_par02_im, 
                      win = owin(c(0,1000),c(0,500)), nsim = 1)
  #spe2.sim
  #plot(spe2.sim)
  
  beta.30 <- 0.02
  mu_par03 <- exp(cov.mat%*%beta.sim[3,] + beta.30)
  mu_par03_mat <- matrix(mu_par03, nrow = 25, ncol = 50, byrow = FALSE)
  mu_par03_im <- im(mu_par03_mat, 
                    xcol=seq(0,1000,length=50), yrow=seq(0,500,length=25))
  #plot(mu_par03_im)
  spe3.sim <- rThomas(kappa = k.par, scale = w.par, mu = mu_par03_im, 
                      win = owin(c(0,1000),c(0,500)), nsim = 1)
  #spe3.sim
  #plot(spe3.sim)
  
  
  beta.40 <- 0.2
  mu_par04 <- exp(cov.mat%*%beta.sim[4,] + beta.40)
  mu_par04_mat <- matrix(mu_par04, nrow = 25, ncol = 50, byrow = FALSE)
  mu_par04_im <- im(mu_par04_mat, 
                    xcol=seq(0,1000,length=50), yrow=seq(0,500,length=25))
  #plot(mu_par04_im)
  spe4.sim <- rThomas(kappa = k.par, scale = w.par, mu = mu_par04_im, 
                      win = owin(c(0,1000),c(0,500)), nsim = 1)
  #spe4.sim
  #plot(spe4.sim)
  
  
  beta.50 <- 0.2
  mu_par05 <- exp(cov.mat%*%beta.sim[3,] + beta.50)
  mu_par05_mat <- matrix(mu_par05, nrow = 25, ncol = 50, byrow = FALSE)
  mu_par05_im <- im(mu_par05_mat, 
                    xcol=seq(0,1000,length=50), yrow=seq(0,500,length=25))
  #plot(mu_par05_im)
  spe5.sim <- rThomas(kappa = k.par, scale = w.par, mu = mu_par05_im, 
                      win = owin(c(0,1000),c(0,500)), nsim = 1)
  #spe5.sim
  #plot(spe5.sim)
  
  ############################################################
  ##################     several ppp       ###################
  ############################################################
  
  sim.spec1 <- EachPPPsimu(spe1.sim, p = 15, m = 10)
  sim.spec2 <- EachPPPsimu(spe2.sim, p = 15, m = 10)
  sim.spec3 <- EachPPPsimu(spe3.sim, p = 15, m = 10)
  sim.spec4 <- EachPPPsimu(spe4.sim, p = 15, m = 10)
  sim.spec5 <- EachPPPsimu(spe5.sim, p = 15, m = 10)
  
  
  # initial beta
  sim_grl_init <- c(sim.spec1$Beta_init,sim.spec2$Beta_init,sim.spec3$Beta_init,
                    sim.spec4$Beta_init,sim.spec5$Beta_init)
  sim_grl_X <- bdiag(sim.spec1$X,sim.spec2$X,sim.spec3$X,sim.spec4$X,sim.spec5$X)
  sim_grl_Y <- c(sim.spec1$Y,sim.spec2$Y,sim.spec3$Y,sim.spec4$Y,sim.spec5$Y)
  
  index <- c(NA, rep(c(NA,1:15),5))
  sim_grl_init <- c(0, sim_grl_init)
  sim_grl_X <- as.matrix(cbind(rep(1,dim(sim_grl_X)[1]),sim_grl_X))
  
  lambda = exp(seq(log(20), log(0.005), length = 50))
  sim_fitlasso <- grplasso(x = sim_grl_X, y = sim_grl_Y, index = index, 
                           coef.init = sim_grl_init, model = LinReg(), lambda = lambda)
  sim_grp.est <- sim_fitlasso$coefficients[-1,]
  
  sim_grps <- list()
  magni <- rep(0,p)
  for (i in 2:(p+1)){
    sim_grps[[i-1]] <- sim_grp.est[(1:5-1)*16+i,]
    magni[i-1] <- sum((sim_grps[[i-1]][,30])^2)
  }
  
  res[l,] <- magni
}

save(nsim, k.par, w.par, q.par, beta.sim, res, ord, file=paste("sim1.rda",sep=""))



