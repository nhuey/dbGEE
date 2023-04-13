############ Needed functions ############ 
check.oracle <- function(alpha.hat, alpha.0){
  nonzeros <- which(alpha.0 != 0)
  zeroes <- which(alpha.0 == 0)
  
  check.nonzero <- (alpha.hat[nonzeros] != 0)
  check.zeroes <- (alpha.hat[zeros] == 0)
  
  return(all(c(check.nonzero, check.zeroes)))
}

check.all.vecs <- function(alpha.hat, gamma.hat, beta.hat){
  a <- check.oracle(alpha.hat, alpha.0)
  g <- check.oracle(gamma.hat, gamma.0)
  b <- check.oracle(beta.hat, beta.0)
  
  return(all(a,g,b))
}

check.non <- function(alpha.hat, gamma.hat, beta.hat){
  nonzeros.a <- which(alpha.0 != 0)
  nonzeros.g <- which(gamma.0 != 0)
  nonzeros.b <- which(beta.0 != 0)
  
  return(mean(c(alpha.hat[nonzeros.a]!=0,gamma.hat[nonzeros.g]!=0,beta.hat[nonzeros.b]!=0)))
}

norm.error <- function(alpha, gamma, beta, norm.type){
  mean(norm((matrix(alpha, nrow = 1) - matrix(alpha.0, nrow = 1)), type = norm.type),
       norm((matrix(gamma, nrow = 1) - matrix(gamma.0, nrow = 1)), type = norm.type),
       norm((matrix(beta, nrow = 1) - matrix(beta.0, nrow = 1)), type = norm.type))
}

reprow <- function(x, n) {
  s <- NROW(x)
  if (length(n) == 1) {
    return(matrix(x[rep(1:s, each = n), ], nrow = s * n, ncol = NCOL(x)))
  }
  matrix(x[rep(1:s, n), ], nrow = sum(n), ncol = NCOL(x))
}

set_up_cv_longit <- function(N, n_folds, id.vect){
  list(
    fold_id = set_up_cv_longit.helper1(N, n_folds, id.vect),
    outlist = as.list(seq(n_folds))
  )
}

set_up_cv_longit.helper1 <- function(N, n_folds, id.vect){
fold_id <- rep(0, length(id.vect))
addition <- sample(1:length(unique(id.vect)), replace = FALSE)
sample.count <- 1
for (fold in 1:(n_folds-1)){
  fold.size <- 0
    while(fold.size <= floor(N/n_folds)){
     fold_id[which(id.vect == addition[sample.count])] <- fold
     fold.size <- fold.size + sum(id.vect == addition[sample.count])
     sample.count <- sample.count + 1
    }
}
fold_id[fold_id==0] <- n_folds
return(fold_id)
}

cv_gds_longit <- function (X, y, family, no_lambda = 10, lambda = NULL, 
          n_folds = 5, weights = rep(1, length(y)), R_hat) 
{
  #browser()
  id.vect <- (as.numeric((X[, 1])))
  X <- X[, -1]
  stopifnot(family %in% c("gaussian", "poisson", "binomial"))
  if (is.null(lambda)) {
    if (no_lambda == 1) {
      lambda <- glmnet::cv.glmnet(X, y, family = family, 
                                  weights = weights)$lambda.min
    }
    else if (no_lambda > 1) {
      lambda <- glmnet::glmnet(X, y, family = family, 
                               nlambda = no_lambda, weights = weights)$lambda
      no_lambda <- length(lambda)
    }
  }
  stopifnot(all(lambda >= 0))
  n <- nrow(X)
  cv_list <- set_up_cv_longit(n, n_folds, id.vect)
  loss <- matrix(nrow = no_lambda, ncol = n_folds)
  for (i in seq(n_folds)) {
    test <- (cv_list$fold_id == i)
    #browser()
    loss[, i] <- as.numeric(lapply(lambda, function(l) {
      fit <- gmus_gee(W = X[!test, , drop = FALSE], y = y[!test], 
                  lambda = l, delta = 0, family = family, weights = weights[!test], id.vect[!test], R_hat = R_hat)
      linpred <- X[test, , drop = FALSE] %*% fit$beta + 
        fit$intercept
      if (family == "binomial") {
        test_pred <- logit(linpred)
      } else if (family == "gaussian") {
        test_pred <- linpred
      } else if (family == "poisson") {
        test_pred <- exp(linpred)
      }
      deviance(as.vector(y[test]), as.vector(test_pred), family, weights[test])
    }))/sum(test)
  }
  cv <- data.frame(lambda = lambda, mean_loss = rowMeans(loss), 
                   sd_loss = apply(loss, 1, stats::sd), upper_1se = rowMeans(loss) + 
                     apply(loss, 1, stats::sd)/sqrt(n_folds), lower_1se = rowMeans(loss) - 
                     apply(loss, 1, stats::sd)/sqrt(n_folds))
  ind1 <- min(which(cv$mean_loss == min(cv$mean_loss)))
  lambda_min <- cv$lambda[ind1]
  loss_min <- cv$mean_loss[ind1]
  ind2 <- min(which(min(cv$mean_loss) > cv$lower_1se))
  lambda_1se <- cv$lambda[ind2]
  loss_1se <- cv$mean_loss[ind2]
  fit <- list(cv = cv, lambda_min = lambda_min, loss_min = loss_min, 
              lambda_1se = lambda_1se, loss_1se = loss_1se, family = family)
  class(fit) <- "cv_gds"
  return(fit)
}

expit <- function(x){
  num <- exp(x)
  denom <- 1+exp(x)
  return(num/denom)
}
##########################################
generate.working.correlation <- function(mydata, structure, betas){
  n <- length(as.vector(table(mydata$id)))
  nt <- as.vector(table(mydata$id))
  cor.struc <- mycor_gee2(N = n, nt = as.vector(table(mydata$id)), y = mydata$y, X = mydata[,-c(1,2)], 
             family = gaussian(link = "identity"), beta_new = betas, corstr = structure, maxclsz = max(as.vector(table(mydata$id))))
  cor.list <-list()
  for(i in 1:n){
    cor.list[[i]] <- (cor.struc$Ehat)[1:nt[i],1:nt[i],i]
  }
  return(cor.list)
}

linprogPD.nathans.version <- function(x0, A, b, epsilon, pdtol=1e-3, pdmaxiter=50) {
  ## Solves
  ## min ||x||_1   subject to  ||Ax-b||_\infty <= epsilon
  ## Adapted from Matlab code for Dantzig Selector by J. Romberg
  N <- length(x0)
  x0 <- matrix(x0, nrow=N)
  b <- matrix(b, nrow=N)
  
  alpha=0.01;
  beta=0.5;
  mu = 10;
  
  gradf0 <- matrix(rep(c(0,1), each=N), nrow=2*N)
  
  if (max(abs(A%*%x0 - b)) > epsilon) {
    stop("Infeasible starting point!")
  }
  
  x <- x0
  u <- 0.95*abs(x0) + 0.1*max(abs(x0))
  
  Atr <- A%*%x - b
  
  fu1 <- x - u
  fu2 <- -x - u
  fe1 <- Atr - epsilon
  fe2 <- -Atr - epsilon
  lamu1 <- -1/fu1
  lamu2 <- -1/fu2
  lame1 <- -1/fe1
  lame2 <- -1/fe2
  
  AtAv <- t(A)%*%(lame1 - lame2)
  
  sdg <- - sum(c(fu1, fu2, fe1, fe2)* c(lamu1, lamu2, lame1, lame2))
  
  tau <- mu*(4*N)/sdg
  
  rdual <- gradf0 + c(lamu1 - lamu2 + AtAv, -lamu1 - lamu2)
  rcent <- -c(lamu1*fu1, lamu2*fu2, lame1*fe1, lame2*fe2) - 1/tau
  resnorm <- sqrt(sum(rdual^2,rcent^2))
  
  pditer <- 0
  done <- (sdg < pdtol) | (pditer >= pdmaxiter)
  
  while(!done) {
    w2 <- -1 - (1/fu1 + 1/fu2)/tau
    
    sig11 <- -lamu1/fu1 - lamu2/fu2
    sig12 <- lamu1/fu1 - lamu2/fu2
    siga <- -(lame1/fe1 + lame2/fe2)
    sigx <- sig11 - sig12^2/sig11
    
    w1 <- -( t(A)%*%(1/fe2 - 1/fe1) + 1/fu2 - 1/fu1)/tau
    w1p <- w1 - (sig12/sig11)*w2
    
    Hp <- t(A)%*%diag(as.vector(siga))%*%A + diag(as.vector(sigx))
    
    if (rcond(Hp) < 1e-14) {
      warning("Ill conditioned matrix.  Previous iterate matrix returned! (May increase perturb/lambda.)")
      xp <- x
      return(xp)
    }
    
    dx <- solve(Hp, w1p)
 
    AtAdx <- A%*%dx
    
    du <- w2/sig11 - (sig12/sig11)*dx
    
    dlamu1 <- -(lamu1/fu1)*(dx-du) - lamu1 - 1/(fu1*tau)
    dlamu2 <- -(lamu2/fu2)*(-dx-du) - lamu2 - 1/(fu2*tau)
    
    dlame1 <- -(lame1/fe1)*(AtAdx) - lame1 - 1/(fe1*tau)
    dlame2 <- (lame2/fe2)*(AtAdx) - lame2 - 1/(fe2*tau)
    
    AtAdv <- t(A)%*%(dlame1 - dlame2)
    
    
    iu1 <- dlamu1 < 0
    iu2 <- dlamu2 < 0
    ie1 <- dlame1 < 0
    ie2 <- dlame2 < 0
    ifu1 <- (dx-du) > 0
    ifu2 <- (-dx-du) > 0
    ife1 <- AtAdx > 0
    ife2 <- AtAdx < 0
    
    smax <- min( -lamu1[iu1]/dlamu1[iu1], -lamu2[iu2]/dlamu2[iu2], -lame1[ie1]/dlame1[ie1], -lame2[ie2]/dlame2[ie2], -fu1[ifu1]/(dx[ifu1] - du[ifu1]), -fu2[ifu2]/(-dx[ifu2] -du[ifu2]), -fe1[ife1]/AtAdx[ife1], -fe2[ife2]/( - AtAdx[ife2])   )
    smax <- min(1,smax)
    
    s <- 0.99*smax
    
    suffdec <- FALSE
    backiter <- 0
    
    while(!suffdec) {
      xp <- x + s*dx
      up <- u + s*du
      Atrp <- Atr + s*AtAdx
      AtAvp <- AtAv+s*AtAdv
      
      fu1p <- fu1 + s*(dx - du)
      fu2p <- fu2 + s*(-dx-du)
      fe1p <- fe1 + s*AtAdx
      fe2p <- fe2 + s*(-AtAdx)
      
      lamu1p <- lamu1 + s*dlamu1
      lamu2p <- lamu2 + s*dlamu2
      lame1p <- lame1 + s*dlame1
      lame2p <- lame2 + s*dlame2
      
      rdp <- gradf0 + c(lamu1p - lamu2p + AtAvp, -lamu1p - lamu2p)
      rcp <- -c(lamu1p*fu1p, lamu2p*fu2p, lame1p*fe1p, lame2p*fe2p) - 1/tau
      suffdec <- sqrt( sum(rdp^2, rcp^2)) < (1 - alpha*s)*resnorm
      
      s <- beta*s
      backiter <- backiter+1
      
      if (backiter > 32) {
        warning("Backtracking stuck.  Previous iterate matrix returned!")
        xp <- x
        return(xp)
      }
    }
    x <- xp
    u <- up
    Atr <- Atrp
    AtAv <- AtAvp
    fu1 <- fu1p
    fu2 <- fu2p
    fe1 <- fe1p
    fe2 <- fe2p
    lamu1 <- lamu1p
    lamu2 <- lamu2p
    lame1 <- lame1p
    lame2 <- lame2p
    
    sdg <- -sum( c(fu1, fu2, fe1, fe2)*c(lamu1, lamu2, lame1, lame2))
    tau <- mu*(4*N)/sdg
    rdual <- rdp
    rcent <- -c(fu1, fu2, fe1, fe2)*c(lamu1, lamu2, lame1, lame2) - 1/tau
    resnorm <- sqrt(sum(rdual^2, rcent^2))
    
    pditer <- pditer+1
    done <- (sdg < pdtol) | (pditer >= pdmaxiter)
  }
  
  return(xp)
}


mycor_gee2 <-
  function(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R = NULL,scale.fix =FALSE,scale.value=1) {
    #browser()
    X <- as.matrix(X)
    eta=X%*%beta_new[-1] +beta_new[1]
    mu=family$linkinv(eta)
    sd=sqrt(family$variance(mu))
    res<-(as.vector(y)-mu)/sd 
    
    if(scale.fix==0) {
      fi<-sum(res^2)/(sum(nt))
    } else
      if(scale.fix==1) {
        fi<-scale.value
      }
    
    aindex=cumsum(nt)
    index=c(0,aindex[-length(aindex)])
    
    if (corstr=="independence") 
    {alfa_hat<-0} else 
      if (corstr=="exchangeable") 
      {
        sum1<-0
        sum3<-0
        for ( i in  1:N)          {
          for ( j in  1:nt[i])      {    
            for ( jj in 1:nt[i])      {    
              if  ( j!=jj)              {
                #cat("i",i,"j",j,"jj",jj,"\n")
                sum2<-res[j+index[i]]*res[jj+index[i]]
                #cat("i",i,"j",j,"jj",jj,"\n")
                sum1<-sum1+sum2
                #cat("i",i,"j",j,"jj",jj,"sum2",sum2,"sum1",sum1,"\n")
              }
            }
          }
          sum4<-nt[i]*(nt[i]-1)
          sum3<-sum3+sum4
        } #i
        alfa_hat<-sum1/(sum3*fi)
      } else
        if (corstr=="AR-1") 
        { 
          sum5<-0
          sum6<-0
          for ( i in  1:N)           {
            for ( j in  1:nt[i])       {  
              for ( jj in 1:nt[i])       {  
                if( j>jj && abs(j-jj)==1)  {
                  #cat("i",i,"j",j,"jj",jj,"\n")
                  sum7<-res[j+index[i]]*res[jj+index[i]]
                  sum5<-sum5+sum7           
                  #cat("i",i,"j",j,"jj",jj,"sum7",sum7,"sum5", sum5, "\n")
                }
              }
            }
            sum8<-(nt[i]-1)
            sum6<-sum6+sum8
          } #i
          alfa_hat<-sum5/(sum6*fi)
        } else
          if (corstr=="stat_M_dep") 
          {  
            alfa_hat=matrix(0,Mv,1)
            for(m in 1:Mv) {
              sum12<-0
              sum14<-0
              for ( i in  1:N)           {
                for ( j in  1:nt[i])       {  
                  for ( jj in 1:nt[i])       {  
                    if( j>jj && abs(j-jj)==m)  {
                      #cat("m",m,"i",i,","j",j,"jj",jj,"\n") 
                      sum11<-res[j+index[i]]*res[jj+index[i]]
                      sum12<-sum12+sum11 
                      #cat("m",m,"i",i,"j",j,"jj",jj,"sum11",sum11,"sum12", sum12, "\n")          
                    } #if
                  }
                }
                sum13<-nt[i]-1
                sum14<-sum14+sum13
              } #i
              alfa_hat[m]<-sum12/(sum14*fi) 
            } #m
          }  else
            if (corstr=="non_stat_M_dep") 
            {  
              alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
              for( m in 1:Mv)            {
                for ( j in  1:nt[1])       {  
                  for ( jj in 1:nt[1])       {  
                    if( j>jj && abs(j-jj)==m)  { 
                      sum16<-0                 
                      for ( i  in 1:N)     {
                        #cat("m",m,"j",j,"jj",jj,"i",i"\n") 
                        sum15<-res[j+index[i]]*res[jj+index[i]]
                        sum16<-sum15+sum16          
                        #cat("j",j,"jj",jj,"i",i,"sum15",sum15,"sum16",sum16,"\n")
                      } #i
                      #cat("j",j,"jj",jj,"sum16",sum16,"\n")
                      alfa_hat[j,jj]<-sum16/(N*fi) 
                    }
                  }
                }
              }
            } else
              if (corstr=="unstructured") 
              {  
                alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
                for ( j in 1:nt[1])  {  
                  for ( jj in 1:nt[1]) {
                    sum20<-0                
                    if (j > jj)          {
                      for ( i  in 1:N )    {
                        #cat("i",i,"j",j,"jj",jj,"\n") 
                        sum21<-res[j+index[i]]*res[jj+index[i]]
                        sum20<-sum21+sum20           
                      } #i
                      #cat("j",j,"jj",jj,"sum20",sum20,"\n")
                      alfa_hat[j,jj]<-sum20/(N*fi) 
                    }
                  }
                }
              } else
                if (corstr=="fixed")
                {alfa_hat=NULL
                }
    
    Ehat<-array(0,c(maxclsz,maxclsz,N))
    
    for(i in 1:N){
      cor1<-matrix(0,nt[i],nt[i])
      if (corstr=="independence")                                        
      {cor1<-diag(nt[i])} else
        if (corstr=="exchangeable")                                        
        { for (t1 in 1:nt[i]) {
          for (t2 in 1:nt[i]) {
            if (t1!=t2) 
            {cor1[t1,t2]<-alfa_hat} else 
            {cor1[t1,t2]<-1}
          }
        }
        } else
          if (corstr=="AR-1")                                      
          { for (t1 in 1:nt[i]) {
            for (t2 in 1:nt[i]) {
              cor1[t1,t2]<-alfa_hat^abs(t1-t2)   
            }
          }
          }  else
            if (corstr=="stat_M_dep")                                     
            { for (t1 in 1:nt[i]) {
              for (t2 in 1:nt[i]) {
                if (abs(t1-t2)==0)
                {cor1[t1,t2]<-1} else
                  for(m in 1:Mv) {
                    if (abs(t1-t2)==m)
                    {cor1[t1,t2]<-alfa_hat[m]} 
                  }
              }
            }
            } else
              if (corstr=="non_stat_M_dep")                                     
              { 
                cor1=alfa_hat+t(alfa_hat)
                diag(cor1)=1
              } else
                if (corstr=="unstructured")                                        
                {cor1=alfa_hat+t(alfa_hat)
                diag(cor1)=1
                } else
                  if (corstr=="fixed")
                  {cor1=R
                  }
      
      Ehat[1:nt[i],1:nt[i],i]<-cor1 
    }
    
    return(list("Ehat"=Ehat,"fi"=fi))
    
  }

clime.nathans.version <- function(x, lambda=NULL,
                  nlambda=ifelse(is.null(lambda),100,length(lambda)),
                  lambda.max=0.8,
                  lambda.min=ifelse(nrow(x)>ncol(x), 1e-4, 1e-2),
                  sigma=FALSE,
                  perturb=TRUE,
                  standardize=TRUE,
                  logspaced=TRUE,
                  linsolver=c("primaldual", "simplex"),
                  pdtol=1e-3, pdmaxiter=50
)
{
  lpfun <- match.arg(linsolver, c("primaldual", "simplex"))
  
  if (sigma) {
    if (is.matrix(x)) {
      Sigma <- x
    } else {
      Sigma <-as.matrix(x)
    }
    p <- ncol(Sigma)
    x <- NULL
  } else {
    n <- nrow(x)
    p <- ncol(x)
    
    if (is.null(lambda)) {
      if (logspaced) {
        lambda <- 10^(seq(log10(lambda.min), log10(lambda.max), length.out=nlambda))
      } else {
        lambda <- seq(lambda.min, lambda.max, length.out=nlambda)
      }
    }
    
    
    if (standardize)  x <- scale(x)
    Sigma <- cov(x)*(1-1/n)
  }
  
  ## Set to perturbed Sigma to have conditional number p
  eigvals <- eigen(Sigma, only.values=T)$values
  if (is.logical(perturb)) {
    if (perturb) { 
      perturb <- max(max(eigvals) - p*min(eigvals), 0)/(p-1)
    } else {
      perturb <- 0
    }
  }
  
  Sigma <- Sigma+diag(p)*perturb
  emat <- diag(p)
  
  Omegalist <- vector("list", nlambda)
  if (lpfun == "simplex") {
    for (jl in 1:nlambda) {
      Omega <- matrix(0, nrow=p, ncol=p)
      lam <- lambda[jl]
      for (j in 1:p) {
        beta <- linprogS(Sigma, emat[,j], lam)
        Omega[,j] <- beta
      }
      Omegalist[[jl]] <- Omega*(abs(Omega)<= abs(t(Omega)))+ t(Omega)*(abs(Omega)> abs(t(Omega)))
    }
  }
  
  if (lpfun == "primaldual") {
    Omega0 <- solve(Sigma)
    
    for (jl in 1:nlambda) {
      Omega <- matrix(0, nrow=p, ncol=p)
      lam <- lambda[jl]
      for (j in 1:p) {
        beta <- linprogPD.nathans.version(Omega0[,j], Sigma, emat[,j], lam, pdtol, pdmaxiter)
        Omega[,j] <- beta
      }
      Omegalist[[jl]] <- Omega*(abs(Omega)<= abs(t(Omega)))+ t(Omega)*(abs(Omega)> abs(t(Omega)))
    }
  }
  
  
  
  
  outlist <- list(Omegalist=Omegalist, x = x, lambda = lambda, perturb=perturb, standardize = standardize, lpfun=lpfun)
  class(outlist) <- c("clime")
  return(outlist)
}
