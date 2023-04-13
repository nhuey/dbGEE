gds_gee <- function(X, y, lambda = NULL, family = "gaussian", weights = NULL, R_hat) {
  
  stopifnot(family %in% c("gaussian", "poisson", "binomial"))
  
  if(!is.null(lambda) & length(lambda) != 1) stop("lambda must be a single value")
  
  stopifnot(lambda >= 0)
  m <- (as.numeric(table(X[, 1])))[1]
  X <- X[, -1]
  fit <- gmus_gee(X, y, lambda = lambda, delta = 0, family = family, weights = weights, R_hat = R_hat, id.vect = id.vect)
  
  # In the Dantzig selector case, delta is not of interest
  fit$delta <- NULL
  
  class(fit) <- "gds"
  return(fit)
}

gmus_gee <- function(W, y, lambda = NULL, delta = NULL,
                 family = "gaussian", weights = NULL, id.vect, R_hat) {
  #browser()
  family <- match.arg(family, choices = c("gaussian", "binomial", "poisson"))
  
  if(!is.null(weights) & length(weights) != nrow(W)) stop("weights vector must be one value per case")
  
  if(is.null(lambda)) {
    lambda <- glmnet::cv.glmnet(W, y, family = family)$lambda.min
  } else {
    stopifnot(all(lambda >= 0))
  }
  if(is.null(delta)) {
    delta <- seq(from = 0, to = 0.5, by = 0.02)
  } else {
    stopifnot(all(delta >= 0))
  }
  
  n <- length(unique(W[,1])) # number of clusters
  id.vect.w <- id.vect
  p <- dim(W)[2] + 1
  W <- scale(W)
  scales <- attr(W, "scaled:scale")
  W <- cbind(rep(1,n), W)
  
  if(family == "gaussian") {
    fit <- sapply(delta, function(delta) mus_glm_gee(W, y, lambda, delta, family, weights, id.vect.w = id.vect.w, R_hat = R_hat))
  } else if(family %in% c("binomial", "poisson")) {
    fit <- sapply(delta, function(delta) mus_glm_gee(W, y, lambda, delta, family, weights, id.vect.w = id.vect.w, R_hat = R_hat))
  }
  
  
  fit <- list(intercept = fit[1, ],
              beta = matrix(fit[2:p, ] / scales, nrow = p - 1),
              family = family,
              delta = delta,
              lambda = lambda,
              num_non_zero = colSums(abs(fit[2:p, , drop = FALSE]) > 1e-10)
  )
  
  class(fit) <- "gmus"
  return(fit)
}


mus_glm_gee <- function(W, y, lambda, delta, family = c("binomial", "poisson", "gaussian"), weights = NULL, id.vect.w, R_hat){
  
  family <- match.arg(family)
  if(family == "gaussian"){
    mu <- identity
  }else if(family == "binomial") {
    mu <- logit
    dmu <- dlogit
  } else if(family == "poisson") {
    mu <- pois
    dmu <- dpois
  }
  
  n <- length(unique(id.vect.w)) 
  p <- dim(W)[2]
  
  bOld <- stats::rnorm(p)/p
  bNew <- stats::rnorm(p)/p
  IRLSeps <- 1e-7
  maxit <- 500
  count <- 1
  Diff1 <- 1
  Diff2 <- 1
  
  while(Diff1 > IRLSeps & Diff2 > IRLSeps & count < maxit){
    bOlder <- bOld
    bOld <- bNew
    if(family == "gaussian"){
      V <- (y - mu(W%*%bOld))^2 # Should a common variance be estimated? Yes, I think so.
      # Function to take in id.vect and return a 'time' vector
      d.V <- data.frame(time = unlist(sapply(as.vector(table(id.vect.w)), seq)), V = V)
      A_i.neghalf <- diag(((d.V %>% group_by(time) %>% summarise(mean = mean(V)))$mean)^(-0.5))
      A_i.poshalf <- diag(((d.V %>% group_by(time) %>% summarise(mean = mean(V)))$mean)^(0.5))
      # A_i.neghalf <- A_i.poshalf <- diag(m)
      #print(A_i.neghalf)
      #print(A_i.poshalf)
    }else{V <- dmu(W%*%bOld)
    }
    z <- vector()
    z.til <- vector()
    W.til <- W
    for(i in 1:n){
      ids <- unique(id.vect.w)
      #print(i)
      #browser()
      Y_i <- y[id.vect.w == ids[i]]
      W_i <- W[id.vect.w == ids[i],]
     # A_i.neghalf <- diag((V[(((i*m)-(m-1)):((i*m)))])^(-0.5))
     # A_i.poshalf <- diag((V[(((i*m)-(m-1)):((i*m)))])^(0.5))
      mu_i <- mu(W%*%bOld)[id.vect.w == ids[i]]
      #print(dim(A_i.neghalf))
      #browser()
      z_i <- W_i%*%bOld + A_i.poshalf[seq(table(id.vect.w)[i]), seq(table(id.vect.w)[i])]%*%solve(R_hat[[ids[i]]])%*%A_i.neghalf[seq(table(id.vect.w)[i]), seq(table(id.vect.w)[i])]%*%(Y_i - mu_i) 
      #z_i <- W_i%*%bOld + solve(R_hat)%*%(Y_i - mu_i) 
      z <- c(z,z_i)
      z.til <- c(z.til, A_i.poshalf[seq(table(id.vect.w)[i]), seq(table(id.vect.w)[i])]%*%z_i )
      W_i.tilde <- A_i.poshalf[seq(table(id.vect.w)[i]), seq(table(id.vect.w)[i])]%*% W_i
      W.til[id.vect.w == ids[i],] <- W_i.tilde
    }
    #Wtilde <- c(sqrt(V)) * W
    #ztilde <- c(sqrt(V)) * c(z)
    if(family == "gaussian"){
      bNew <- musalgorithm_gee(W, z, lambda, delta * sqrt(sum((V)^2)) / sqrt(n), weights = NULL)
    }else{
    bNew <- musalgorithm_gee(W.til, z.til, lambda, delta * sqrt(sum((V)^2)) / sqrt(n), weights = NULL)}
    
    count <- count+1
    Diff1 <- sum(abs(bNew - bOld))
    Diff2 <- sum(abs(bNew - bOlder))
    #print(Diff1)
    #print(bNew)
  }
  if(count >= maxit) print("First Phase did not converge")
  
  #count <- 1
  #Diff1 <- 1
  #Diff2 <- 1
  
  #while(Diff1 > IRLSeps & Diff2 > IRLSeps & count < maxit){
   # bOlder <- bOld
   # bOld <- bNew
   # if(family == "gaussian"){
   #   V <- (y - mu(W%*%bOld))^2 # Should a common variance be estimated? Yes, I think so.
   #   d.V <- data.frame(time = rep(1:m,n),V = V)
   #   A_i.neghalf <- diag(((d.V %>% group_by(time) %>% summarise(mean = mean(V)))$mean)^(-0.5))
   #   A_i.poshalf <- diag(((d.V %>% group_by(time) %>% summarise(mean = mean(V)))$mean)^(0.5))
   #   # A_i.neghalf <- A_i.poshalf <- diag(m)
   #   #print(A_i.neghalf)
   #   #print(A_i.poshalf)
   # }else{V <- dmu(W%*%bOld)
   # }
   # z <- vector()
   # z.til <- vector()
   # W.til <- W
   # for(i in 1:n){
   #   #print(i)
   #   Y_i <- y[(((i*m)-(m-1)):((i*m)))]
   #   W_i <- W[(((i*m)-(m-1)):((i*m))),]
   #   A_i.neghalf <- diag((V[(((i*m)-(m-1)):((i*m)))])^(-0.5))
   #   A_i.poshalf <- diag((V[(((i*m)-(m-1)):((i*m)))])^(0.5))
   #   mu_i <- mu(W%*%bOld)[(((i*m)-(m-1)):((i*m)))]
   #   #print(dim(A_i.neghalf))
   #   #browser()
   #   z_i <- W_i%*%bOld + A_i.neghalf%*%solve(R_hat)%*%A_i.neghalf%*%(Y_i - mu_i) 
   #   #z_i <- W_i%*%bOld + solve(R_hat)%*%(Y_i - mu_i) 
   #   z <- c(z,z_i)
   #   z.til <- c(z.til, A_i.poshalf%*%z_i )
   #   W_i.tilde <- A_i.poshalf%*% W_i
   #   W.til[(((i*m)-(m-1)):((i*m))),] <- W_i.tilde
   # }
   # #Wtilde <- c(sqrt(V)) * W
   # #ztilde <- c(sqrt(V)) * c(z)
   # if(family == "gaussian"){
   #   bNew <- musalgorithm_gee(W, z, lambda, delta * sqrt(sum((V)^2)) / sqrt(n), weights = NULL)
   # }else{
   #   bNew <- musalgorithm_gee(W.til, z.til, lambda, delta * sqrt(sum((V)^2)) / sqrt(n), weights = NULL)}
   # 
   # count <- count+1
   # Diff1 <- sum(abs(bNew - bOld))
   # Diff2 <- sum(abs(bNew - bOlder))
   # #print(Diff1)
   # #print(bNew)
  #}
  #if(count >= maxit) print("Second Phase did not converge")
  return(bNew)
}


musalgorithm_gee <- function(W, y, lambda, delta, weights = NULL){
  # We assume the first column of W is constants, i.e., intercept
  n <- dim(W)[1]
  p <- dim(W)[2]
  obj <- c(rep(1,p),rep(0,p))
  mat <- matrix(0,nrow=4*p, ncol=2*p)
  
  # Weight matrix
  if(is.null(weights)){
    D <- diag(n)
  } else {
    D <- diag(weights)
  }
  
  
  # Inequality constraint, -u_j - beta_j <= 0
  mat[1:p,1:p] <- -diag(p)
  mat[1:p,(p+1):(2*p)] <- -diag(p)
  
  # Inequality constraint, -u_j + beta_j <= 0
  mat[(p+1):(2*p),1:p] <- -diag(p)
  mat[(p+1):(2*p),(p+1):(2*p)] <- diag(p)
  
  # First "score function" constraint
  mat[(2*p+1),1:p] <- matrix(0, nrow=1, ncol=p)
  mat[(2*p+2):(3*p),1:p] <- matrix(-delta, nrow=(p-1), ncol=p)
  mat[(2*p+1):(3*p),(p+1):(2*p)] <- 1/n*(t(W) %*% D %*% W)
  
  # Second "score function" constraint
  mat[(3*p+1),1:p] <- matrix(0, nrow=1, ncol=p)
  mat[(3*p+2):(4*p),1:p] <- matrix(-delta, nrow=(p-1), ncol=p)
  mat[(3*p+1):(4*p),(p+1):(2*p)] <- -1/n*(t(W) %*% D %*% W)
  
  rhs <- rep(0,(4*p))
  rhs[(2*p+1)] <- 1/n*(t(W[,1]) %*% D %*% y)
  rhs[(2*p+2):(3*p)] <- lambda + 1/n*(t(W[,-1]) %*% D %*% y)
  rhs[(3*p+1)] <- -1/n*(t(W[,1]) %*% D %*% y)
  rhs[(3*p+2):(4*p)] <- lambda - 1/n*(t(W[,-1]) %*% D %*% y)
  dir <- rep("<=",4*p)
  bounds <- list(lower=list(ind=1:(2*p), val=rep(-Inf,2*p)),
                 upper=list(ind=1:(2*p), val=rep(Inf,2*p)))
  
  
  
  
  bhat <- Rglpk::Rglpk_solve_LP(obj = obj, mat = mat, dir = dir,
                                rhs = rhs, bounds = bounds)$solution
  
  
  return(bhat[(p+1):(2*p)])
  
  
}
