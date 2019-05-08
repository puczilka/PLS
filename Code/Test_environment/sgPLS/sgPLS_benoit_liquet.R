##### gPLS

# function (X, Y, ncomp, mode = "regression", max.iter = 500, tol = 1e-06, 
# keepX, keepY = NULL, ind.block.x, ind.block.y = NULL, scale = TRUE) 
# {

scale = FALSE
tol = 1e-06
max.iter = 1
mode = "regression"
ncomp = 1

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  q <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) 
    {Y.names = "Y"} else {
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
  mat.c <- matrix(nrow = p, ncol = ncomp)
  mat.d <- matrix(nrow = q, ncol = ncomp)
  mat.e <- matrix(nrow = q, ncol = ncomp)
  mat.t <- matrix(nrow = n, ncol = ncomp)
  mat.u <- matrix(nrow = n, ncol = ncomp)
  X.s <- scale(X, scale = scale)
  Y.s <- scale(Y, scale = scale)
  sparsity.x <- length(ind.block.x) + 1 - keepX
  if (is.null(ind.block.y)) {
    sparsity.y <- rep(0, ncomp)
  } else {
    if (is.null(keepY)) 
      keepY <- rep(length(ind.block.y) + 1, ncomp)
    sparsity.y <- length(ind.block.y) + 1 - keepY
  }
  
  ###
  
  res.load <- step1.group.spls.sparsity(X = X.s, Y.s, ind.block.x = ind.block.x, 
                                        ind.block.y = ind.block.y, sparsity.x = sparsity.x[1], 
                                        sparsity.y = sparsity.y[1], epsilon = tol, iter.max = max.iter)
  
  ###
  
  res.deflat <- step2.spls(X = X.s, Y = Y.s, res.load$u.tild.new, 
                           res.load$v.tild.new, mode = mode)
  
  ###
  
  mat.c[, 1] <- res.deflat$c
  iter <- res.load$iter
  if (mode == "regression") 
    {mat.d[, 1] <- res.deflat$d} else mat.e[, 1] <- res.deflat$e
  load.u <- res.load$u.tild.new
  load.v <- res.load$v.tild.new
  mat.t[, 1] <- X.s %*% load.u
  mat.u[, 1] <- Y.s %*% load.v
  if (ncomp > 1) {
    for (h in 2:ncomp) {
      res.load <- step1.group.spls.sparsity(X = res.deflat$X.h, 
                                            Y = res.deflat$Y.h, ind.block.x = ind.block.x, 
                                            ind.block.y = ind.block.y, sparsity.x = sparsity.x[h], 
                                            sparsity.y = sparsity.y[h], epsilon = tol, iter.max = max.iter)
      load.u <- cbind(load.u, res.load$u.tild.new)
      load.v <- cbind(load.v, res.load$v.tild.new)
      mat.t[, h] <- res.deflat$X.h %*% res.load$u.tild.new
      mat.u[, h] <- res.deflat$Y.h %*% res.load$v.tild.new
      res.deflat <- step2.spls(X = res.deflat$X.h, Y = res.deflat$Y.h, 
                               res.load$u.tild.new, res.load$v.tild.new, mode = mode)
      mat.c[, h] <- res.deflat$c
      if (mode == "regression") 
        mat.d[, h] <- res.deflat$d
      else mat.e[, h] <- res.deflat$e
      iter <- c(iter, res.load$iter)
    }
  } else {
    load.u <- matrix(load.u, ncol = 1)
    load.v <- matrix(load.v, ncol = 1)
  }
  
  ###
  
  cl = match.call()
  if (is.null(keepY)) {
    if (is.null(ind.block.y)) 
      keepY <- rep(ncol(Y), ncomp)
    else keepY <- rep(length(ind.block.y) + 1, ncomp)
  }
  rownames(load.u) <- X.names
  rownames(load.v) <- Y.names
  dim = paste("comp", 1:ncomp)
  colnames(load.u) = colnames(load.v) = dim
  result <- list(call = cl, X = X.s, Y = Y.s, ncomp = ncomp, 
                 mode = mode, keepX = keepX, keepY = keepY, mat.c = mat.c, 
                 mat.d = mat.d, mat.e = mat.e, loadings = list(X = load.u, 
                                                               Y = load.v), variates = list(X = mat.t, Y = mat.u), 
                 names = list(X = X.names, Y = Y.names, indiv = ind.names), 
                 tol = tol, max.iter = max.iter, iter = iter, ind.block.x = ind.block.x, 
                 ind.block.y = ind.block.y)
  class(result) = c("gPLS", "sPLS", "spls", "pls")
  # return(invisible(result))
  # }
  
  
  # # =======================================================
  # 
  # step1.group.spls.sparsity <- function(X,Y,ind.block.x,ind.block.y,sparsity.x,sparsity.y,epsilon,iter.max){
  #   n <- dim(X)[1]
  #   Z <- t(X)%*%Y
  #   svd.Z <- svd(Z,nu=1,nv=1)
  #   
  #   u.tild.old <- svd.Z$u
  #   v.tild.old <- svd.Z$v
  #   u.tild.previous <- v.tild.previous <- 0
  #   iter <- 0
  #   
  #   ### Step c
  #   #|(norm(v.tild.old-v.tild.previous)>epsilon)
  #   while (((normv(u.tild.old-u.tild.previous)>epsilon) ) & (iter <iter.max))  {
  #     
  #     vecZV <- Z%*%matrix(v.tild.old,ncol=1)
  #     tab.ind <- c(0,ind.block.x,length(vecZV))      # Add 0 and p into ind.block.x
  #     res <- NULL
  #     # Iterate over each group
  #     for (i in 1:(length(ind.block.x)+1)){
  # 
  #       # Define the size of each group, p_k
  #       ji <- tab.ind[i+1]-tab.ind[i]
  # 
  #       # Return the entries corresponding to each group of the projection of the loadings
  #       vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
  # 
  #       res <- c(res,2*normv(vecx)/sqrt(ji))     # Create vector of 'multipliers' for group lasso computation
  #     }
  # 
  #     # Calculate sparsity based on keepX
  #     if(sparsity.x==0) lambda.x <- 0 else{
  #       lambda.x <- sort(res)[sparsity.x]}
  #     
  #     u.tild.new <- soft.thresholding.group(Z%*%matrix(v.tild.old,ncol=1),ind=ind.block.x,lambda=lambda.x)
  #     u.tild.new <- u.tild.new/sqrt(sum(u.tild.new**2))
  #     
  #     if(sparsity.y==0) {lambda.y <- 0} else { 
  #       vecZV <- t(Z)%*%matrix(u.tild.new,ncol=1)
  #       tab.ind <- c(0,ind.block.y,length(vecZV))
  #       res <- NULL
  #       for (i in 1:(length(ind.block.y)+1)){
  #         ji <- tab.ind[i+1]-tab.ind[i]  
  #         vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
  #         res <- c(res,2*normv(vecx)/sqrt(ji))
  #       }
  #       
  #       lambda.y <- sort(res)[sparsity.y]}
  #     
  #     
  #     
  #     if(sparsity.y==0) {v.tild.new <- t(Z)%*%matrix(u.tild.new,ncol=1)} else {
  #       v.tild.new <- soft.thresholding.group(t(Z)%*%matrix(u.tild.new,ncol=1),ind=ind.block.y,lambda=lambda.y)}
  #     
  #     v.tild.new <- v.tild.new/sqrt(sum(v.tild.new**2))
  #     
  #     u.tild.previous <- u.tild.old
  #     v.tild.previous <- v.tild.old
  #     
  #     u.tild.old <- u.tild.new
  #     v.tild.old <- v.tild.new
  #     
  #     iter <- iter +1
  #   }  
  #   res <- list(iter=iter, u.tild.new=u.tild.new,v.tild.new=v.tild.new) 
  #   
  # }
  # 
  # 
  # # =======================================================
  # 
  # 
  # 
  # soft.thresholding.group <- function(x,ind,lambda){
  #   tab.ind <- c(0,ind,length(x))
  #   tol <- .Machine$double.eps ^ 0.5 
  #   res <- NULL
  #   for (i in 1:(length(ind)+1)){
  #     ji <- tab.ind[i+1]-tab.ind[i]  
  #     vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
  #     y <- 1-(lambda/2)*sqrt(ji)/normv(vecx)
  #     if(y < tol) y <- 0
  #     res <- c(res,vecx*y)  
  #   }
  #   return(res)    
  # }
  