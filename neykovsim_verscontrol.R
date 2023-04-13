# Starting to code the Neykov estimation procedure
# First just need to import code for Dantzig selector and CLIME
# Step 1: Simulate data with longitudinal stucture (I should already have this code)
# Step 2: Estimate betas with Dantzig selector - Lambda = AKsqrt(log(d)/n), 
## A > 2sqrt(2) and K = subgaussian norm of epsilon
# Step 3: Estimate v via CLIME - Lambda' >= ||v*||_1 * 4 * A_x * K^2_X * sqrt(log(d)/n)
## A_x*sqrt(log(d)/n) leq 1 and K_X = subgaussian norm of row of X
# Step 4: Solve estimating equation for each single coefficient

j  <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(j)

args = commandArgs(TRUE)
n <- as.numeric(args[1])
pn <- as.numeric(args[2])
s <- as.numeric(args[3])
b <- as.numeric(args[4])
rho <- as.numeric(args[5])
#lambda_prime <- as.numeric(args[6])
m <- as.numeric(args[6])

## Step 1: Simulate data with longitudinal stucture
# required R packages and extra codes
library(mvtnorm)
library(gee)
library(geepack)
#library(fastclime)
library(PGEE)
library(tidyverse)
library(cvTools)
library(flare)
library(hdme)
library(clime)
library(matchingR)
library(matrixNormal)
library(expm)
source("/n/home08/nhuey/neykov/v7/quadform.R")
source("/n/home08/nhuey/neykov/v7/helper_functions.R")
source("/n/home08/nhuey/neykov/v7/helper_functions2.R")
source("/n/home08/nhuey/neykov/v7/loss_functions.R")
source("/n/home08/nhuey/neykov/v7/project_gradient.R")
source("/n/home08/nhuey/neykov/v7/gee_gds.R")
source("/n/home08/nhuey/neykov/climetuningv7/clime_cv_longit.R")
# True beta
beta.true <- c(0.5,rep(b,s),rep(0,pn-s-1))
coverage <- ci.length <- neykov.var <- vector(length = pn)
  
  # vector if subject ids
  id.vect <- rep(1:n, each = m) 
  
  # covariance matrix of (pn) number of continuous covariates 
  X.sigma.column <- matrix(0,(pn),(pn))
  {
    for (i in 1:(pn))
      X.sigma.column[i,] <- 0.5^(abs((1:(pn))-i))  
  }
  
 X.sigma.row <- matrix(0, m, m)
  { 
    for (i in 1:m)
      X.sigma.row[i,] <- 0.5^(abs((1:m)-i))
  }

  M <- matrix(0, nrow = m, ncol = pn)
  x.mat <- x.mat.2 <- matrix(NA, nrow = n*m, ncol = pn)
  # generate matrix of covariates (time varying)
  for(i in 1:n){
    x.mat[((i-1)*m+1):(m*i),] <- matrixNormal::rmatnorm(M = M, U = X.sigma.row, V = X.sigma.column)
    x.mat.2[((i-1)*m+1):(m*i),] <- matrixNormal::rmatnorm(M = M, U = X.sigma.row, V = X.sigma.column)
  }    
  
  # Normalize X's according to assumptions (for independence working correlation)
  ## Do I still need this?
  # for(i in 1:n){
  #   x.i <- x.mat[((m*i)-(m-1)):(m*i),]
  #   for(k in 1:pn){
  #     if(sum((x.i[,k])^2) > 1){ x.i[,k] <- (x.i[,k])/(sum((x.i[,k])^2))}
  #   }
  # }
  
  # Center and unit scale design matrix
  x.mat <- apply(x.mat, 2, scale)
  x.mat.2 <- apply(x.mat.2, 2, scale)

  # true values
  sigma2 <- 1
  if(m == 1){R <- matrix(rho,m,m)+rep(1-rho,m)}else{
    R <- matrix(rho,m,m)+diag(rep(1-rho,m))}
  
  # covariance matrix of error
  SIGMA <- sigma2*R
  error <- rmvnorm(n, mean = rep(0,m),SIGMA)
  error.2 <- rmvnorm(n, mean = rep(0,m), SIGMA)

  # generate longitudinal data with continuous outcomes
  y.temp <- x.mat%*%beta.true
  y.temp.2 <- x.mat.2%*%beta.true
  y.vect <- y.temp+as.vector(t(error))
  y.vect.2 <- y.temp.2+as.vector(t(error.2))
  #y.vect <- scale(y.vect)
  
  mydata <- data.frame(id.vect,y.vect,x.mat)
  mydata.2 <- data.frame(id.vect, y.vect.2, x.mat.2) 
  colnames(mydata) <- c("id","y",paste("x",1:length(beta.true),sep = ""))
  colnames(mydata.2) <- c("id","y",paste("x",1:length(beta.true),sep = ""))
  
  # Step 2: Estimate betas with Dantzig selector - Lambda = AKsqrt(log(d)/n), 
  ## A > 2sqrt(2) and K = subgaussian norm of epsilon
  
  A <- 2*sqrt(2)
  K <- 1
  lambda <- A*K*sqrt(log(pn)/(n*m))
  invisible(capture.output(R_hat <- generate.working.correlation(mydata.2, structure = "exchangeable")))
  invisible(capture.output(R_hat.2 <- generate.working.correlation(mydata, structure = "exchangeable")))
  
  #tuning <- list(tuning.psi = c(3.443689, 4.685061))
  
  #dantzig.model <- dantzig.cv(x.mat, y.vect, lambda, nlambda = 100)
  #beta_hat <- dantzig.model$beta.cv
  
  ## Add ID column to x.mat
  x.mat<- cbind(id.vect,x.mat)
  x.mat.2 <- cbind(id.vect, x.mat.2)
  cv_fit <- cv_gds_longit(x.mat, y.vect, family = "gaussian", no_lambda = 50, n_folds = 10, R_hat = diag(m))
  cv_fit.2 <- cv_gds_longit(x.mat.2, y.vect.2, family = "gaussian", no_lambda = 50, n_folds = 10, R_hat = diag(m))
  # Remove ID column
  x.mat <- x.mat[,-1]
  x.mat.2 <- x.mat.2[,-1]
  fit <- gds(x.mat, y.vect, lambda = cv_fit$lambda_min, family = "gaussian")
  fit.2 <- gds(x.mat.2, y.vect.2, lambda = cv_fit.2$lambda_min, family = "gaussian")
  #fit <- gds(x.mat, y.vect, lambda = lambda, family = "gaussian")
  beta_hat <- fit$beta
  beta_hat.2 <- fit.2$beta
  
  #dantzig_run=slim(X=x.mat,Y=y.vect,method="dantzig",nlambda = 10,verbose=TRUE)
  # beta.cv <-  (dantzig_run$beta)
  # 
  # sum(abs(beta_hat))
  # sum(abs(beta.cv))
  # 
  # (1/n)*max(abs(t(x.mat)%*%(x.mat%*%beta_hat - y.vect)))
  # (1/n)*max(abs(t(x.mat)%*%(x.mat%*%beta.cv - y.vect)))
  # (1/n)*max(abs(t(x.mat)%*%(x.mat%*%fit$beta - y.vect)))
  
  # Step 3: Estimate v via CLIME - Lambda' >= ||v*||_1 * 4 * A_x * K^2_x * sqrt(log(d)/n)
  ## A_x*sqrt(log(d)/n) leq 1 and K_X = subgaussian norm of row of X
  
  
  #Sigma_n <- (1/(n*m)) * (crossprod(x.mat))
  
  #A_x <- 0.01/(sqrt(log(pn)/(n*m)))
  #K_x <- 1

  
  # Fast CLIME algo. I don't think the same lambda is used for each column
  # Phi_clime <- fastclime(Sigma_n, lambda.min = lambda_prime, nlambda = 100)
  # Phi_c <- fastclime.selector(Phi_clime$lambdamtx, Phi_clime$icovlist, lambda_prime)
  # Phi <- Phi_c$icov
  
  resids <- y.vect.2 - (x.mat.2 %*% beta_hat.2)
  V <- (resids)^2
  d.V <- data.frame(time = rep(1:m,n),V = V)
  A_i.neghalf <- diag(((d.V %>% group_by(time) %>% summarise(mean = mean(V)))$mean)^(-0.5))
  #A_i.neghalf <- diag(m)
  R_hat.neg0.5 <- sqrtm(solve(R_hat))
  x.mat.forclime <- x.mat
  
  resids.2 <- y.vect - (x.mat %*% beta_hat)
  V.2 <- (resids.2)^2
  d.V.2 <- data.frame(time = rep(1:m,n),V = V.2)
  A_i.neghalf.2 <- diag(((d.V.2 %>% group_by(time) %>% summarise(mean = mean(V)))$mean)^(-0.5))
  #A_i.neghalf <- diag(m)
  R_hat.neg0.5.2 <- sqrtm(solve(R_hat.2))
  x.mat.forclime.2 <- x.mat.2
  
  # Transform design matrix for appropriate Hessian
  resid_matrix <- matrix(data = 0, nrow = m, ncol = m)
  
  for(i in 1:n){
    ep.hat.i <- resids[((i*m)-(m-1)):((i*m))]
    resid_matrix <- resid_matrix + outer(ep.hat.i, ep.hat.i)
  }
  resids.piece <- ((1/n)*resid_matrix)
  
  x.piece <- matrix(data = 0, nrow = pn, ncol = pn)
  for(i in 1:n){
    x.mat.i <- x.mat[(((i*m)-(m-1)):((i*m))),]
    #A_i.neghalf <- diag((V[(((i*m)-(m-1)):((i*m)))])^(-0.5))
    x.transform <- R_hat.neg0.5 %*% A_i.neghalf 
    x.mat.forclime[((i*m)-(m-1)):((i*m)), ] <- x.transform %*% x.mat.i
    A.piece <- A_i.neghalf %*% solve(R_hat) %*% A_i.neghalf
    x.temp <- t(x.mat.i) %*% A.piece %*% resids.piece %*% A.piece %*% x.mat.i
    x.piece <- x.piece + x.temp
  }
  x.piece <- (1/n)*x.piece
  
  # May want to use GLasso for continuous data eventually
  # Phi.glasso <- CVglasso(x.mat.forclime, cores = 2)
  # Phi <- Phi.glasso$Omega
  tuning.params <- seq(0.01, 0.5, length.out = 50)
  #Phi_clime.opt <- clime(x.mat.forclime, lambda = lambda_prime, standardize = FALSE)
  re.clime <- clime(x.mat.forclime, lambda = tuning.params, standardize = FALSE)
  re.cv <- clime.cv.longit(re.clime, loss = "tracel2", fold=10, m=m)
  re.clime.opt <- clime(x.mat.forclime, re.cv$lambdaopt)
  Phi <- (re.clime.opt$Omegalist)[[1]]
  
  neykov_beta <- vector(length = pn)
  neykov_variance <- vector(length = pn)
  lower <- upper <- vector(length = pn)
  
    # Step 4: Solve estimating equation for each single coefficient
  for(i in 1:pn){
    S_hat.theta.gamma_hat <- function(x){ beta_hat_i <- beta_hat
    beta_hat_i[i] <- x
    x.for.solve <- matrix(0, nrow = n*m, ncol = pn)
    for(j in 1:n){
      x.mat.i <- x.mat[(((j*m)-(m-1)):((j*m))),]
     # A_i.neghalf <- diag((V[(((j*m)-(m-1)):((j*m)))])^(-0.5))
     # A_i.neghalf <- diag(m)
      #print(dim(x.for.solve[(((j*m)-(m-1)):((j*m))),]))
      x.for.solve[(((j*m)-(m-1)):((j*m))),] <- (A_i.neghalf %*% solve(R_hat) %*% A_i.neghalf %*% x.mat.i)
    }
    return(t(Phi[,i]) %*% t(x.for.solve) %*%((x.mat%*%beta_hat_i)-y.vect))} 
    neykov_solve <- uniroot(S_hat.theta.gamma_hat, c(-1000, 1000))
    neykov_beta[i] <- neykov_solve$root
    # May need to modify for m = 1
      neykov_variance[i] <- (1/m) * (quadform(x.piece, Phi[,i])) 
    
    lower[i] <- neykov_beta[i] - (qnorm(.975)*sqrt(neykov_variance[i]/(n*m)))
    upper[i] <- neykov_beta[i] + (qnorm(.975)*sqrt(neykov_variance[i]/(n*m)))
  }
  
  
  
  ######################################
  resid_matrix.2 <- matrix(data = 0, nrow = m, ncol = m)
  
  for(i in 1:n){
    ep.hat.i.2 <- resids.2[((i*m)-(m-1)):((i*m))]
    resid_matrix.2 <- resid_matrix.2 + outer(ep.hat.i.2, ep.hat.i.2)
  }
  resids.piece.2 <- ((1/n)*resid_matrix.2)
  
  x.piece.2 <- matrix(data = 0, nrow = pn, ncol = pn)
  for(i in 1:n){
    x.mat.i.2 <- x.mat.2[(((i*m)-(m-1)):((i*m))),]
    #A_i.neghalf <- diag((V[(((i*m)-(m-1)):((i*m)))])^(-0.5))
    x.transform.2 <- R_hat.neg0.5.2 %*% A_i.neghalf.2 
    x.mat.forclime.2[((i*m)-(m-1)):((i*m)), ] <- x.transform.2 %*% x.mat.i.2
    A.piece.2 <- A_i.neghalf.2 %*% solve(R_hat.2) %*% A_i.neghalf.2
    x.temp.2 <- t(x.mat.i.2) %*% A.piece.2 %*% resids.piece.2 %*% A.piece.2 %*% x.mat.i.2
    x.piece.2 <- x.piece.2 + x.temp.2
  }
  x.piece.2 <- (1/n)*x.piece.2
  
  # May want to use GLasso for continuous data eventually
  # Phi.glasso <- CVglasso(x.mat.forclime, cores = 2)
  # Phi <- Phi.glasso$Omega
  tuning.params <- seq(0.01, 0.5, length.out = 50)
  #Phi_clime.opt <- clime(x.mat.forclime, lambda = lambda_prime, standardize = FALSE)
  re.clime.2 <- clime(x.mat.forclime.2, lambda = tuning.params, standardize = FALSE)
  re.cv.2 <- clime.cv.longit(re.clime.2, loss = "tracel2", fold=10, m=m)
  re.clime.opt.2 <- clime(x.mat.forclime.2, re.cv.2$lambdaopt)
  Phi.2 <- (re.clime.opt.2$Omegalist)[[1]]
  
  neykov_beta.2 <- vector(length = pn)
  neykov_variance.2 <- vector(length = pn)
  lower.2 <- upper.2 <- vector(length = pn)
  
  # Step 4: Solve estimating equation for each single coefficient
  for(i in 1:pn){
    S_hat.theta.gamma_hat <- function(x){ beta_hat_i.2 <- beta_hat.2
    beta_hat_i.2[i] <- x
    x.for.solve.2 <- matrix(0, nrow = n*m, ncol = pn)
    for(j in 1:n){
      x.mat.i.2 <- x.mat.2[(((j*m)-(m-1)):((j*m))),]
      # A_i.neghalf <- diag((V[(((j*m)-(m-1)):((j*m)))])^(-0.5))
      # A_i.neghalf <- diag(m)
      #print(dim(x.for.solve[(((j*m)-(m-1)):((j*m))),]))
      x.for.solve.2[(((j*m)-(m-1)):((j*m))),] <- (A_i.neghalf.2 %*% solve(R_hat.2) %*% A_i.neghalf.2 %*% x.mat.i.2)
    }
    return(t(Phi.2[,i]) %*% t(x.for.solve.2) %*%((x.mat.2%*%beta_hat_i.2)-y.vect.2))} 
    neykov_solve.2 <- uniroot(S_hat.theta.gamma_hat, c(-1000, 1000))
    neykov_beta.2[i] <- neykov_solve.2$root
    # May need to modify for m = 1
    neykov_variance.2[i] <- (1/m) * (quadform(x.piece.2, Phi.2[,i])) 
    
    lower.2[i] <- neykov_beta.2[i] - (qnorm(.975)*sqrt(neykov_variance.2[i]/(n*m)))
    upper.2[i] <- neykov_beta.2[i] + (qnorm(.975)*sqrt(neykov_variance.2[i]/(n*m)))
  }
  
  ######################################
  
  neykov_beta <- mean(c(neykov_beta, neykov_beta.2))
  neykov_variance <- mean(c(neykov_variance, neykov_variance.2))
  lower <- mean(c(lower, lower.2))
  upper <- mean(c(upper, upper.2))
  coverage <- (beta.true >= lower & beta.true <= upper)
  ci.length <- upper - lower
  neykov.var <- neykov_variance
  
  a <- cbind(neykov_beta, neykov.var, lower, upper, coverage)
  
  write.table(a, paste("/n/home08/nhuey/neykov/v7/N100n",n,"p",pn,"s",s,"b",b,"r",rho,"m",m,"/neykov_",j,".txt", sep=""))
