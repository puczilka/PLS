# Load mixOMICs and linnerud dataset
library('mixOmics')

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological

# test <- spls(X, Y, scale = TRUE, ncomp = n_comp, mode = "regression")

# #-------- SUBSETTING --------#
# # Take subset of Y so X and Y have different dimensions (easier for analysing algorithm)
# Y <- Y[,1:2]
# #----------------------------#

# Scale data
X_h <- as.matrix(scale(X))
Y_h <- as.matrix(scale(Y))


#-------------------- MANUAL sPLS (sparsity parameter) --------------------#

# Select number of components
n_comp <- 2

# Select sparsity parameters (penalties)
lambda_u <- 0     # lambda 1
lambda_v <- 0     # lambda 2

# Create items to store results
Eta_mat <- NULL     # X scores
Omega_mat <- NULL     # Y scores
U_mat <- NULL     # X loadings
V_mat <- NULL     # Y loadings
C <- NULL     # Regression coefficient for X
D <- NULL     # Regression coefficient for Y in 'Regression mode'
E <- NULL     # Regression coefficient for Y in 'PLS-mode A'


###### LOOP OVER COMPONENTS ######
# Initiate first component and loop over defined number of components
for (h in (1:n_comp)) {
  
  # Compute matrix M (p x q matrix)
  M <- t(X_h) %*% Y_h
  
  ###### CALIBRATION OF COMPONENT ######
  # Find SVD of M and find loadings from first pair of singular vectors
  M_decomp <- svd(M, nu = 1, nv = 1)
  u_old <- M_decomp$u
  v_old <- M_decomp$v
  
  # Initiate iterables
  counter <- 0
  tol <- 1e-06
  max_iter <- 100
  u_diff <- length(X_h)
  v_diff <- length(X_h)
  
  # Loop until convergence of u and v or max iteration
  while ((sum(abs(u_diff) > tol) != 0 || sum(abs(v_diff) > tol) != 0) && counter < max_iter + 1) {
    
    # Optimise u and v iteratively using soft-tresholding penalisation function (and normalise)
    u_new = sign(M %*% v_old) * max(0, (abs(M %*% v_old) - lambda_v))
    u_new = u_new / sqrt(as.vector(t(u_new) %*% u_new))
    v_new = sign(t(M) %*% u_old) * max(0, (abs(t(M) %*% u_old) - lambda_u))
    v_new = v_new / sqrt(as.vector(t(v_new) %*% v_new))
    
    # Update iterables
    u_diff = u_new - u_old
    v_diff = v_new - v_old
    u_old = u_new
    v_old = v_new
    counter = counter + 1
    
  }
  
  if (counter == max_iter + 1) {
    print(paste0('Warning! Max iteration reached. No convergence for component ', h))
  }
  
  ###### DEFLATION STEP ######
  # Calculate scores/latent variables for X and Y
  eta = as.vector(X_h %*% u_new) / as.vector(t(u_new) %*% u_new)
  omega = as.vector(Y_h %*% v_new) / as.vector(t(v_new) %*% v_new)
  
  # Calculate regression coefficients
  c = as.vector(t(X_h) %*% eta) / as.vector(t(eta) %*% eta)
  d = as.vector(t(Y_h) %*% eta) / as.vector(t(eta) %*% eta)
  e = as.vector(t(Y_h) %*% omega) / as.vector(t(omega) %*% omega)
  
  
  # Deflate X and Y matrices using latent variables and regression coefficients
  X_h <- X_h - (eta %*% t(c))
  Y_h <- Y_h - (eta %*% t(d))
  
  # Store variables
  Eta_mat <- cbind(Eta_mat, eta, deparse.level = 0)
  Omega_mat <- cbind(Omega_mat, omega, deparse.level = 0)
  U_mat <- cbind(U_mat, u_new, deparse.level = 0)
  V_mat <- cbind(V_mat, v_new, deparse.level = 0)
  C <- cbind(C, c, deparse.level = 0)
  D <- cbind(D, d, deparse.level = 0)
  E <- cbind(E, e, deparse.level = 0)
  
}

# Assign row/column names labels
colnames(U_mat) <- colnames(V_mat) <- colnames(Eta_mat) <- colnames(Omega_mat)  <- colnames(C) <- colnames(D) <- colnames(E) <- paste0('comp ', seq(n_comp))
rownames(C) <- colnames(X)
rownames(D) <- rownames(E) <- colnames(Y)
