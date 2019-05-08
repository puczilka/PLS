library(mixOmics)
library(sgPLS)
#---------------------------- SIMULATED DATA ---------------------------------#

## Simulation of datasets X and Y with group variables
n <- 100
sigma.gamma <- 1
sigma.e <- 1.5
p <- 400
q <- 500
theta.x1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5,15), 
              rep(0, 5), rep(-1.5, 15), rep(0, 325))
theta.x2 <- c(rep(0, 320), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
              rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))

theta.y1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5, 15),
              rep(0, 5), rep(-1.5, 15), rep(0, 425))
theta.y2 <- c(rep(0, 420), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
              rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))                            


Sigmax <- matrix(0, nrow = p, ncol = p)
diag(Sigmax) <- sigma.e ^ 2
Sigmay <- matrix(0,nrow = q, ncol = q)
diag(Sigmay) <- sigma.e ^ 2

set.seed(125)

gam1 <- rnorm(n)
gam2 <- rnorm(n)

X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                               nrow = 2, byrow = TRUE) + rmvnorm(n, mean = rep(0, p), sigma =
                                                                                                   Sigmax, method = "svd")
Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2), 
                                                               nrow = 2, byrow = TRUE) + rmvnorm(n, mean = rep(0, q), sigma =
                                                                                                   Sigmay, method = "svd")

ind.block.x <- seq(20, 380, 20)
ind.block.y <- seq(20, 480, 20)
keepX = c(4, 4)
keepY = c(4, 4)


X <- as.matrix(scale(X))
Y <- as.matrix(scale(Y))


#-------------------- MANUAL gPLS (penalising groups of vaiables) --------------------#
# Select number of components
n_comp <- 1
tol <- 1e-06

# gPLS <- function(X,Y,ncomp,mode="regression",max.iter=500,tol=1e-06,keepX ,keepY=NULL,ind.block.x,ind.block.y=NULL,scale=TRUE){


X <- as.matrix(X)
Y <- as.matrix(Y)
q <- ncol(Y)
p <- ncol(X)
n <- nrow(X)
# X.names = dimnames(X)[[2]]
# if (is.null(X.names)) 
#   X.names = paste("X", 1:p, sep = "")
# if (dim(Y)[2] == 1) 
#   Y.names = "Y"
# else {
#   Y.names = dimnames(Y)[[2]]
#   if (is.null(Y.names)) 
#     Y.names = paste("Y", 1:q, sep = "")
# }
# ind.names = dimnames(X)[[1]]
# if (is.null(ind.names)) {
#   ind.names = dimnames(Y)[[1]]
#   rownames(X) = ind.names
# }
# if (is.null(ind.names)) {
#   ind.names = 1:n
#   rownames(X) = rownames(Y) = ind.names
# }
# 
# mat.c <-matrix(nrow = p, ncol = n_comp)
# mat.d <- matrix(nrow = q, ncol = n_comp)
# mat.e <- matrix(nrow = q, ncol = n_comp)
# mat.t <- matrix(nrow = n, ncol = n_comp)
# mat.u <- matrix(nrow = n, ncol = n_comp)
# 
# X.s <- scale(X,scale=scale)
# Y.s <- scale(Y,scale=scale)

sparsity.x <- length(ind.block.x)+1-keepX

if(is.null(ind.block.y)) {sparsity.y <- rep(0,ncomp)} else {
  if (is.null(keepY)) keepY <- rep(length(ind.block.y)+1,ncomp)
  sparsity.y <- length(ind.block.y)+1-keepY}


res.load <- step1.group.spls.sparsity(X=X,Y,ind.block.x=ind.block.x,ind.block.y=ind.block.y,sparsity.x=sparsity.x[1],sparsity.y=sparsity.y[1],epsilon=tol,iter.max=1)
# res.deflat <- step2.spls(X=X,Y=Y,res.load$u.tild.new,res.load$v.tild.new,mode=mode)

# mat.c[,1] <- res.deflat$c
# iter <- res.load$iter
# if (mode=="regression") mat.d[,1] <- res.deflat$d else mat.e[,1] <- res.deflat$e
# load.u <- res.load$u.tild.new
# load.v <- res.load$v.tild.new
# mat.t[, 1] <- X.s%*%load.u
# mat.u[, 1] <- Y.s%*%load.v
# if(ncomp>1) {
#   
#   for (h in 2:ncomp) {
#     res.load <- step1.group.spls.sparsity(X=res.deflat$X.h,Y=res.deflat$Y.h,ind.block.x=ind.block.x,ind.block.y=ind.block.y,sparsity.x=sparsity.x[h],sparsity.y=sparsity.y[h],epsilon=tol,iter.max=max.iter)
#     load.u <- cbind(load.u,res.load$u.tild.new)
#     load.v <- cbind(load.v,res.load$v.tild.new)
#     mat.t[, h] <- res.deflat$X.h%*%res.load$u.tild.new
#     mat.u[, h] <- res.deflat$Y.h%*%res.load$v.tild.new 
#     res.deflat <- step2.spls(X=res.deflat$X.h,Y=res.deflat$Y.h,res.load$u.tild.new,res.load$v.tild.new,mode=mode)
#     mat.c[,h] <- res.deflat$c
#     if (mode=="regression") mat.d[,h] <- res.deflat$d else mat.e[,h] <- res.deflat$e
#     iter <- c(iter,res.load$iter)}
# }else{
#   load.u <- matrix(load.u,ncol=1)
#   load.v <- matrix(load.v,ncol=1)
# }

# mat.t <- X.s%*%load.u
# mat.u <- Y.s%*%load.v
# mat.c regressor for X
# mat.d regressor for Y "regression mode"
# mat.e regressor for Y "canonical mode"
# cl = match.call()
# if (is.null(keepY)){
#   if (is.null(ind.block.y)) keepY <- rep(ncol(Y),ncomp) else keepY <- rep(length(ind.block.y)+1,ncomp)
# } 
# rownames(load.u) <- X.names
# rownames(load.v) <- Y.names
# 
# dim = paste("comp", 1:ncomp)
# colnames(load.u) = colnames(load.v) = dim
# # colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = dim 
# 
# result <- list(call = cl,X=X.s,Y=Y.s,ncomp=ncomp,mode=mode,keepX=keepX,keepY=keepY,mat.c=mat.c,mat.d=mat.d,mat.e=mat.e,loadings = list(X = load.u, Y = load.v),variates = list(X = mat.t, Y = mat.u),
#                names = list(X = X.names,Y = Y.names, indiv = ind.names),tol=tol,max.iter=max.iter,iter=iter,ind.block.x=ind.block.x,ind.block.y=ind.block.y)
# class(result) = c("gPLS","sPLS","spls", "pls")
#   return(invisible(result))
# }





# ========================================================================== #




# step1.group.spls.sparsity <- function(X,Y,ind.block.x,ind.block.y,sparsity.x,sparsity.y,epsilon,iter.max){
Z <- t(X)%*%Y
svd.Z <- svd(Z,nu=1,nv=1)

u.tild.old <- svd.Z$u
v.tild.old <- svd.Z$v
u.tild.previous <- v.tild.previous <- 0
iter <- 0

# ========================================================
epsilon <- 1e-06
# ========================================================

### Step c
#|(norm(v.tild.old-v.tild.previous)>epsilon)
# while (((normv(u.tild.old-u.tild.previous)>epsilon) ) & (iter < max.iter))  {

vecZV <- Z%*%matrix(v.tild.old,ncol=1)
tab.ind <- c(0,ind.block.x,length(vecZV))
res <- NULL
for (i in 1:(length(ind.block.x)+1)){
  ji <- tab.ind[i+1]-tab.ind[i]  
  vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
  res <- c(res,2*normv(vecx)/sqrt(ji))
}
if(sparsity.x[1]==0) lambda.x <- 0 else{
  lambda.x <- sort(res)[sparsity.x[1]]}


# u.tild.new <- soft.thresholding.group(Z%*%matrix(v.tild.old,ncol=1),ind=ind.block.x,lambda=lambda.x)
# ========================================================

tab.ind <- c(0,ind.block.x,length(vecZV))
p.tol <- .Machine$double.eps ^ 0.5 
res <- NULL
for (i in 1:(length(ind.block.x)+1)){
  ji <- tab.ind[i+1]-tab.ind[i]  
  vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
  y <- 1-(lambda.x/2)*sqrt(ji)/normv(vecx)
  if(y < p.tol) y <- 0
  res <- c(res,vecx*y)  
}

u.tild.new <- res
# ========================================================

u.tild.new <- u.tild.new/sqrt(sum(u.tild.new**2))



if(sparsity.y[1]==0) {lambda.y <- 0} else { 
  vecZV <- t(Z)%*%matrix(u.tild.new,ncol=1)
  tab.ind <- c(0,ind.block.y,length(vecZV))
  res <- NULL
  for (i in 1:(length(ind.block.y)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]  
    vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
    res <- c(res,2*normv(vecx)/sqrt(ji))
  }
  
  lambda.y <- sort(res)[sparsity.y[1]]}



if(sparsity.y[1]==0) {v.tild.new <- t(Z)%*%matrix(u.tild.new,ncol=1)} else {
  
  # v.tild.new <- soft.thresholding.group(t(Z)%*%matrix(u.tild.new,ncol=1),ind=ind.block.y,lambda=lambda.y)}
  # ========================================================
  tab.ind <- c(0,ind.block.y,length(vecZV))
  p.tol <- .Machine$double.eps ^ 0.5 
  res <- NULL
  for (i in 1:(length(ind.block.y)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]  
    vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
    y <- 1-(lambda.y/2)*sqrt(ji)/normv(vecx)
    if(y < p.tol) y <- 0
    res <- c(res,vecx*y)  
  }
  
  v.tild.new <- res
  # ========================================================
  
  v.tild.new <- v.tild.new/sqrt(sum(v.tild.new**2))
  
  v.new2 <- soft.thresholding.group(t(Z)%*%matrix(u.tild.new,ncol=1),ind=ind.block.y,lambda=lambda.y)
  v.new2 <- v.new2/sqrt(sum(v.new2**2))
}
  


  u.tild.previous <- u.tild.old
  v.tild.previous <- v.tild.old
  
  u.tild.old <- u.tild.new
  v.tild.old <- v.tild.new
  
  iter <- iter +1
# }  
res <- list(iter=iter, u.tild.new=u.tild.new,v.tild.new=v.tild.new) 

# }



soft.thresholding.group <- function(x,ind,lambda){
  tab.ind <- c(0,ind,length(x))
  tol <- .Machine$double.eps ^ 0.5 
  res <- NULL
  for (i in 1:(length(ind)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]  
    vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
    y <- 1-(lambda/2)*sqrt(ji)/normv(vecx)
    if(y < tol) y <- 0
    res <- c(res,vecx*y)  
  }
  return(res)    
}




# step2.spls <- function(X,Y,u.tild.new,v.tild.new,mode){
#   ### Step d
#   xi.h <- X%*% matrix(u.tild.new,ncol=1)/((normv(u.tild.new))**2)
#   w.h  <- Y%*% matrix(v.tild.new,ncol=1)/((normv(v.tild.new))**2)
#   
#   ### Step e
#   c.h <- t(X)%*%matrix(xi.h,ncol=1)/((normv(xi.h))**2)
#   
#   d.rh <- t(Y)%*%matrix(xi.h,ncol=1)/(sum(xi.h*xi.h))
#   
#   d.h <- t(Y)%*%matrix(w.h,ncol=1)/(sum(w.h*w.h))
#   
#   ###Step f and g
#   X.h <- X - xi.h%*%t(c.h)
#   if (mode=="regression") Y.h <- Y - xi.h%*%t(d.rh) else Y.h <- Y - w.h%*%t(d.h)
#   
#   res <- list(X.h=X.h,Y.h=Y.h,c=c.h,d=d.rh,e=d.h)
#   return(res)
# }






model <- gPLS(X,Y,ncomp = n_comp, mode="regression", max.iter=1, tol=1e-06, keepX=keepX, keepY=keepY, ind.block.x=ind.block.x, ind.block.y=ind.block.y, scale=TRUE)

