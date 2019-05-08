sgPLS <- function(X,Y,ncomp,mode="regression",max.iter=500,tol=1e-06,keepX,keepY=NULL,ind.block.x,ind.block.y=NULL,alpha.x,alpha.y=NULL,upper.lambda=10^5,scale=TRUE){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  q <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) 
    Y.names = "Y"
  else {
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names)) 
      Y.names = paste("Y", 1:q, sep = "")
  }
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) {
    ind.names = dimnames(Y)[[1]]
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Y) = ind.names
  }
  
  mat.c <-matrix(nrow = p, ncol = ncomp)
  mat.d <- matrix(nrow = q, ncol = ncomp)
  mat.e <- matrix(nrow = q, ncol = ncomp)
  mat.t <- matrix(nrow = n, ncol = ncomp)
  mat.u <- matrix(nrow = n, ncol = ncomp)
  
  X.s <- scale(X,scale=scale)
  Y.s <- scale(Y,scale=scale)
  
  sparsity.x <- length(ind.block.x)+1-keepX
  if(is.null(ind.block.y)) {sparsity.y <- rep(0,ncomp)} else {
    if (is.null(keepY)) keepY <- rep(length(ind.block.y)+1,ncomp)
    sparsity.y <- length(ind.block.y)+1-keepY}
  
  
  res.load <- step1.sparse.group.spls.sparsity(X=X.s,Y.s,ind.block.x=ind.block.x,ind.block.y=ind.block.y,sparsity.x=sparsity.x[1],sparsity.y=sparsity.y[1],epsilon=tol,iter.max=max.iter,alpha.x=alpha.x[1],alpha.y=alpha.y[1],upper.lambda=upper.lambda)
  res.deflat <- step2.spls(X=X.s,Y=Y.s,res.load$u.tild.new,res.load$v.tild.new,mode=mode)
  
  mat.c[,1] <- res.deflat$c
  iter <- res.load$iter
  if (mode=="regression") mat.d[,1] <- res.deflat$d else mat.e[,1] <- res.deflat$e
  load.u <- res.load$u.tild.new
  load.v <- res.load$v.tild.new
  mat.t[, 1] <- X.s%*%load.u
  mat.u[, 1] <- Y.s%*%load.v
  if(ncomp>1) {
    
    for (h in 2:ncomp) {
      res.load <- step1.sparse.group.spls.sparsity(X=res.deflat$X.h,Y=res.deflat$Y.h,ind.block.x=ind.block.x,ind.block.y=ind.block.y,sparsity.x=sparsity.x[h],sparsity.y=sparsity.y[h],epsilon=tol,iter.max=max.iter,alpha.x=alpha.x[h],alpha.y=alpha.y[h],upper.lambda=upper.lambda)
      load.u <- cbind(load.u,res.load$u.tild.new)
      load.v <- cbind(load.v,res.load$v.tild.new)
      mat.t[, h] <- res.deflat$X.h%*%res.load$u.tild.new
      mat.u[, h] <- res.deflat$Y.h%*%res.load$v.tild.new 
      res.deflat <- step2.spls(X=res.deflat$X.h,Y=res.deflat$Y.h,res.load$u.tild.new,res.load$v.tild.new,mode=mode)
      mat.c[,h] <- res.deflat$c
      if (mode=="regression") mat.d[,h] <- res.deflat$d else mat.e[,h] <- res.deflat$e
      iter <- c(iter,res.load$iter)}
  }else{
    load.u <- matrix(load.u,ncol=1)
    load.v <- matrix(load.v,ncol=1)
  }
  
  rownames(load.u) <- X.names
  rownames(load.v) <- Y.names
  
  #mat.t <- X.s%*%load.u
  #mat.u <- Y.s%*%load.v
  # mat.c regressor for X
  # mat.d regressor for Y "regression mode"
  # mat.e regressor for Y "canonical mode"
  cl = match.call()
  if (is.null(keepY)){
    if (is.null(ind.block.y)) keepY <- rep(ncol(Y),ncomp) else keepY <- rep(length(ind.block.y)+1,ncomp)
  } 
  result <- list(call = cl,X=X.s,Y=Y.s,ncomp=ncomp,mode=mode,keepX=keepX,keepY=keepY,mat.c=mat.c,mat.d=mat.d,mat.e=mat.e,loadings = list(X = load.u, Y = load.v),variates = list(X = mat.t, Y = mat.u),
                 names = list(X = X.names,Y = Y.names, indiv = ind.names),tol=tol,max.iter=max.iter,iter=iter,ind.block.x=ind.block.x,ind.block.y=ind.block.y,alpha.x=alpha.x,alpha.y=alpha.y,upper.lambda=upper.lambda)
  class(result) = c("sgPLS","sPLS", "spls","pls")
  return(invisible(result))
}


# ========================================================================== #


step1.sparse.group.spls.sparsity <- function(X,Y,ind.block.x,ind.block.y,sparsity.x,sparsity.y,epsilon,iter.max,alpha.x,alpha.y,upper.lambda=upper.lambda){
  n <- dim(X)[1]
  Z <- t(X)%*%Y
  svd.Z <- svd(Z,nu=1,nv=1)
  
  u.tild.old <- svd.Z$u
  v.tild.old <- svd.Z$v
  u.tild.previous <- v.tild.previous <- 0
  iter <- 0
  
  ### Step c
  #|(norm(v.tild.old-v.tild.previous)>epsilon)
  while (((normv(u.tild.old-u.tild.previous)>epsilon) ) & (iter <iter.max))  {
    vecZV <- Z%*%matrix(v.tild.old,ncol=1)
    tab.ind <- c(0,ind.block.x,length(vecZV))
    lamb.x <- NULL
    lamb.max <- upper.lambda
    for (i in 1:(length(ind.block.x)+1)){
      ji <- tab.ind[i+1]-tab.ind[i]  
      vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
      lamb.x <- c(lamb.x,uniroot(lambda.quadra,lower=0,upper=lamb.max,vec=vecx,alpha=alpha.x)$root)
    }   
    if(sparsity.x==0){lambda.x <- sort(lamb.x)[1]-1} else {
      lambda.x <- sort(lamb.x)[sparsity.x]}
    
    ####block to zero
    index.block.zero.x <- which(lamb.x<=lambda.x)
    
    
    u.tild.new <- soft.thresholding.sparse.group(Z%*%matrix(v.tild.old,ncol=1),ind=ind.block.x,lambda=lambda.x,alpha=alpha.x,ind.block.zero=index.block.zero.x)
    
    u.tild.new <- u.tild.new/sqrt(sum(u.tild.new**2))
    
    if(sparsity.y==0) {lambda.y <- 0} else { 
      vecZV <- t(Z)%*%matrix(u.tild.new,ncol=1)
      tab.ind <- c(0,ind.block.y,length(vecZV))
      lamb.y <- NULL
      lamb.max <- 100000
      res <- NULL
      for (i in 1:(length(ind.block.y)+1)){
        ji <- tab.ind[i+1]-tab.ind[i]  
        vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
        lamb.y <- c(lamb.y,uniroot(lambda.quadra,lower=0,upper=lamb.max,vec=vecx,alpha=alpha.y)$root)
      }
      lambda.y <- sort(lamb.y)[sparsity.y]
      index.block.zero.y <- which(lamb.y<=lambda.y)
    }
    
    if(sparsity.y==0) {v.tild.new <- t(Z)%*%matrix(u.tild.new,ncol=1)} else {
      v.tild.new <- soft.thresholding.sparse.group(t(Z)%*%matrix(u.tild.new,ncol=1),ind=ind.block.y,lambda=lambda.y,alpha=alpha.y,ind.block.zero=index.block.zero.y)
    }
    
    v.tild.new <- v.tild.new/sqrt(sum(v.tild.new**2))
    
    u.tild.previous <- u.tild.old
    v.tild.previous <- v.tild.old
    
    u.tild.old <- u.tild.new
    v.tild.old <- v.tild.new
    
    iter <- iter +1
  }  
  res <- list(iter=iter, u.tild.new=u.tild.new,v.tild.new=v.tild.new) 
  
}

# ========================================================================== #

lambda.quadra <- function(x,vec,alpha){
  return(sum(soft.thresholding(vec,x*alpha/2)**2)-length(vec)*((1-alpha)*x)**2)
}

# ========================================================================== #

soft.thresholding <- function(x,lambda){
  tol <- .Machine$double.eps ^ 0.5 
  y <- abs(x)-lambda 
  test  <- y < tol
  return(sign(x)*y*(1-test))  
}  

# ========================================================================== #

soft.thresholding.sparse.group <- function(x,ind,lambda,alpha,ind.block.zero){
  tab.ind <- c(0,ind,length(x))
  res <- NULL
  for (i in 1:(length(ind)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]  
    vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
    if(i%in%ind.block.zero) {vecx <- rep(0,ji)} else{
      temp <- soft.thresholding(vecx,lambda*alpha/2)
      vecx <- temp*(1-lambda*(1-alpha)/2*sqrt(length(vecx))/sqrt(sum(temp**2))) 
    }
    res <- c(res,vecx)  
  }
  return(res)    
}

# ========================================================================== #

tol <- .Machine$double.eps ^ 0.5 
y <- abs(vecx)-lamb.max*alpha.x/2 
test  <- y < tol
sign(x)*y*(1-test)  





#==============================================================#
#==============================================================#
#==============================================================#

sg.lambda.quadratic <- function(vec, lambda, alpha, group.size) {
  
  g.x <- sign(vec) * pmax(0, (abs(vec) - (lambda) * alpha/2))
  
  return(sum(g.x**2) - ((lambda_x * (1 - alpha.x)) ** 2) * group.size)
  
}



#==================== Define function to find lambda threshold ====================#

sg.lambda.quadratic <- function(vect, lambda, alpha, blocks, group.sizes) {
  
  temp <- NULL
  
  for (group in (1:length(group.sizes))) {
    
    vec <- vect[(blocks[[group]])]
    g.soft <- sign(vec) * pmax(0, (abs(vec) - (lambda) * alpha/2))
    temp <- c(temp, sum(g.soft**2) - ((lambda_x * (1 - alpha.x)) ** 2) * group.size[group])
    
  }
  
  
  return(temp)
}
#==============================================================#







