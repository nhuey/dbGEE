clime.cv.longit <- function (clime.obj, loss = c("likelihood", "tracel2"), fold = 5, id.vect) 
{
  x <- clime.obj$x
  if (is.null(x)) 
    stop("No x in clime object.  Use x instead of sigma for computing clime!")
  n <- nrow(x)
  p <- ncol(x)
  part.list <- cv.part.longit(n, fold, id.vect)
  lambda <- clime.obj$lambda
  lpfun <- clime.obj$lpfun
  nlambda <- length(lambda)
  lossname <- match.arg(loss, c("likelihood", "tracel2"))
  lossfun <- match.fun(lossname)
  loss.re <- matrix(0, nrow = fold, ncol = nlambda)
  for (j in 1:fold) {
    x.train <- x[part.list$trainMat[[j]], ]
    clime.cv <- clime.nathans.version(x.train, lambda, standardize = FALSE, 
                      perturb = clime.obj$perturb, linsolver = lpfun)
    x.test <- x[part.list$testMat[[j]], ]
    ntest <- nrow(x.test)
    for (jl in 1:nlambda) {
      loss.re[j, jl] <- loss.re[j, jl] + lossfun((cov(x.test) * 
                                                    (1 - 1/ntest)), clime.cv$Omegalist[[jl]])
    }
  }
  loss.mean <- apply(loss.re, 2, mean)
  loss.sd <- apply(loss.re, 2, sd)
  lambdaopt <- lambda[which.min(loss.mean)]
  outlist <- list(lambdaopt = lambdaopt, loss = lossname, lambda = lambda, 
                  loss.mean = loss.mean, loss.sd = loss.sd, lpfun = lpfun)
  class(outlist) <- c("cv.clime")
  return(outlist)
}

cv.part.longit <- function(n, k, id.vect) {
  ntest <- floor(length(unique(id.vect))/k) # number of individuals in test set
  ntrain <- length(unique(id.vect))-ntest # number of individuals in training set
  
  ind <- sample(length(unique(id.vect)))
  
  trainMat <- matrix(NA, nrow=ntrain, ncol=k)
  testMat <- matrix(NA, nrow=ntest, ncol=k)
  trainlist.2 <- list()
  testlist.2 <- list()
  nn <- 1:length(unique(id.vect))
  
  for (j in 1:k) {
    sel <- ((j-1)*ntest+1):(j*ntest)
    testMat[,j] <- ind[sel]
    sel2 <-nn[ !(nn %in% sel) ]
    trainMat[,j] <- ind[sel2]
    trainMat.2 <- vector(length = sum((id.vect %in% ind[sel2])))
    testMat.2 <- vector(length = sum((id.vect %in% ind[sel])))
    front <- 1
    for(i in 1:ntrain){
      ns <- length(which(trainMat[i,j] == id.vect)) -1
      trainMat.2[front:(front + ns)] <- which(trainMat[i,j] == id.vect)
      front <- front + ns + 1
    }
    front <- 1
    for(i in 1:ntest){
      ns <- length(which(testMat[i,j] == id.vect)) -1
      testMat.2[front:(front + ns)] <- which(testMat[i,j] == id.vect)
      front <- front + ns + 1
    }
    trainlist.2[[j]] <- trainMat.2
    testlist.2[[j]] <- testMat.2
  }
  return(list(trainMat=trainlist.2, testMat=testlist.2))
}
