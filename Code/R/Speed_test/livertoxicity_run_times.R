# Load data
library(sgPLS)

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

X <- as.matrix(X)
Y <- as.matrix(Y)

X.s <- scale(X)
Y.s <- scale(Y)

# NOTE: Load manual functions into environment before running

#======================= Model parameters ==========================#
n_components <- 2

ind.block.x <- seq(500, 3000, 500)
ind.block.y <- c(3,7)

keepX = c(4, 4)
keepY = c(2, 2)

keepX2 = c(1500, 1500)
keepY2 = c(7, 7)

alpha.x = c(0.95, 0.95)
alpha.y = c(0.95, 0.95)


##### Iterations
N <- 20



#===================== PLS model ============================#
MO_times <- SG_times <- RS_times <- NULL

for (iter in (1:N)){
  # RS
  t0 <- Sys.time()
  RS_pls <- PLSreg(X.s,Y.s,n_components = n_components)
  t1 <- Sys.time()
  RS_times <- c(RS_times, (t1 - t0))
  
  # MixOmics
  t2 <- Sys.time()
  MO_pls <- pls(X, Y, ncomp = n_components, mode = "regression", scale = TRUE)
  t3 <- Sys.time()
  MO_times <- c(MO_times, (t3 - t2))
  
  # sgPLS
  t4 <- Sys.time()
  SG_pls <- sPLS(X, Y, ncomp = n_components, mode = "regression", scale = TRUE)
  t5 <- Sys.time()
  SG_times <- c(SG_times, (t5 - t4))
}

RS.pls <- mean(RS_times)
MO.pls <- mean(MO_times)
SG.pls <- mean(SG_times)

print(RS.pls)
print(MO.pls)
print(SG.pls)



#===================== sPLS model ============================#
MO_times <- SG_times <- RS_times <- NULL

for (iter in (1:N)){
  # RS
  t0 <- Sys.time()
  RS_spls <- sPLSreg(X.s,Y.s,n_components = n_components, keepX = keepX2,
                     keepY = keepY2)
  t1 <- Sys.time()
  RS_times <- c(RS_times, (t1 - t0))
  
  # MixOmics
  t2 <- Sys.time()
  MO_spls <- spls(X, Y, ncomp = n_components, mode = "regression", keepX = keepX2,
                  keepY = keepY2, scale = TRUE)
  t3 <- Sys.time()
  MO_times <- c(MO_times, (t3 - t2))
  
  # sgPLS
  t4 <- Sys.time()
  SG_spls <- sPLS(X, Y, ncomp = n_components, mode = "regression", keepX = keepX2, 
                  keepY = keepY2, scale = TRUE)
  t5 <- Sys.time()
  SG_times <- c(SG_times, (t5 - t4))
}

RS.spls <- mean(RS_times)
MO.spls <- mean(MO_times)
SG.spls <- mean(SG_times)

print(RS.spls)
print(MO.spls)
print(SG.spls)



#===================== gPLS model ============================#
MO_times <- SG_times <- RS_times <- NULL

for (iter in (1:N)){
  # RS
  t0 <- Sys.time()
  RS_gpls <- gPLSreg(X.s,Y.s,n_components = n_components, keepX_group = keepX,
                     keepY_group = keepY, ind.block.x = ind.block.x, ind.block.y = ind.block.y)
  t1 <- Sys.time()
  RS_times <- c(RS_times, (t1 - t0))
  
  # MixOmics
  # t2 <- Sys.time()
  # MO_gpls <- sPLS(X, Y, ncomp = n_components, mode = "regression", keepX = keepX, 
  #                    keepY = keepY, ind.block.x = ind.block.x , ind.block.y = ind.block.y, scale = TRUE)
  # t3 <- Sys.time()
  # MO_times <- c(MO_times, (t3 - t2))
  
  # sgPLS
  t4 <- Sys.time()
  SG_gpls <- gPLS(X, Y, ncomp = n_components, mode = "regression", keepX = keepX, 
                  keepY = keepY, ind.block.x = ind.block.x , ind.block.y = ind.block.y, scale = TRUE)
  t5 <- Sys.time()
  SG_times <- c(SG_times, (t5 - t4))
}

RS.gpls <- mean(RS_times)
# MO <- mean(MO_times)
SG.gpls <- mean(SG_times)

print(RS.gpls)
print(SG.gpls)



#===================== sgPLS model ============================#
MO_times <- SG_times <- RS_times <- NULL

for (iter in (1:N)){
  # RS
  t0 <- Sys.time()
  RS_sgpls <- sgPLSreg(X.s,Y.s,n_components = n_components, keepX_group = keepX,
                       keepY_group = keepY, ind.block.x = ind.block.x, ind.block.y = ind.block.y, alpha.x = alpha.x, alpha.y = alpha.y)
  t1 <- Sys.time()
  RS_times <- c(RS_times, (t1 - t0))
  
  # MixOmics
  # t2 <- Sys.time()
  # MO_gpls <- sPLS(X, Y, ncomp = n_components, mode = "regression", keepX = keepX, 
  #                    keepY = keepY, ind.block.x = ind.block.x , ind.block.y = ind.block.y, scale = TRUE)
  # t3 <- Sys.time()
  # MO_times <- c(MO_times, (t3 - t2))
  
  # sgPLS
  t4 <- Sys.time()
  SG_sgpls <- sgPLS(X, Y, ncomp = n_components, mode = "regression", keepX = keepX, 
                    keepY = keepY, ind.block.x = ind.block.x , ind.block.y = ind.block.y, alpha.x = alpha.x, alpha.y = alpha.y, scale = TRUE)
  t5 <- Sys.time()
  SG_times <- c(SG_times, (t5 - t4))
}

RS.sgpls <- mean(RS_times)
# MO <- mean(MO_times)
SG.sgpls <- mean(SG_times)

print(RS.sgpls)
print(SG.sgpls)