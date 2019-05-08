# Define ind.block.x/ind.block.y (i.e. vector of indices denoting the end of each group inclusive
# e.g. ind.block.x = c(6, 16) <==> 3 groups s.t. group 1 = 1-6, group 2 = 7-16, group 3 = 17-ncol(X))
ind.block.x <- seq(20, 380, 20)
ind.block.y <- seq(20, 480, 20)

# Select keepX_group/keepY_group variables (i.e. number of groups to keep in each component)
# keepY_group = rep(l, n_components)
keepX = c(4, 4)
keepY = c(4, 4)



# step1.sparse.group.spls.sparsity <- function(X,Y,ind.block.x,ind.block.y,sparsity.x,sparsity.y,epsilon,iter.max,alpha.x,alpha.y,upper.lambda=upper.lambda){

alpha.x = c(0.95, 0.95)
alpha.y = c(0.95, 0.95)
upper.lambda=1e-05

n <- dim(X)[1]
Z <- t(X.s)%*%Y.s
svd.Z <- svd(Z,nu=1,nv=1)

u.tild.old <- svd.Z$u
v.tild.old <- svd.Z$v
u.tild.previous <- v.tild.previous <- 0
iter <- 0

### Step c
#|(norm(v.tild.old-v.tild.previous)>epsilon)
# while (((normv(u.tild.old-u.tild.previous)>epsilon) ) & (iter <iter.max))  {
vecZV <- Z%*%matrix(v.tild.old,ncol=1)
tab.ind <- c(0,ind.block.x,length(vecZV))
lamb.x <- NULL
lamb.max <- upper.lambda
for (i in 1:(length(ind.block.x)+1)){
  ji <- tab.ind[i+1]-tab.ind[i]  
  vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
  lamb.x <- c(lamb.x,uniroot(lambda.quadra,lower=0,upper=lamb.max,vec=vecx,alpha=alpha.x)$root)
}  
# Obtained vector of lambdas s.t. each entry corresponds to the threshold for zeroing out groups (condition 10 in literature)

# We do not know what lambda to choose to penalise groups of variables. This is dependent on the sparse group penalty
# To find lambda, first note that the square of condition 10 gives a quadratic in terms of lambda (hence lambda.quadra).
# We solve the quadratic numerically to obtain a choice of lambda which will penalise the groups according to the sparsity
# parameter

if(sparsity.x==0){lambda.x <- sort(lamb.x)[1]-1} else {
  lambda.x <- sort(lamb.x)[sparsity.x]}
# Can subtract one to account for numerical difference in solution

####block to zero
index.block.zero.x <- which(lamb.x<=lambda.x)
# Vector of Boolean values to determine whether group should be zeroed out.

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
# }  
# res <- list(iter=iter, u.tild.new=u.tild.new,v.tild.new=v.tild.new) 

# }
