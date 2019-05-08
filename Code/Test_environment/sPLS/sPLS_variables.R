# Load mixOMICs and linnerud dataset
library('mixOmics')

data(linnerud)
X_data <- linnerud$exercise
Y_data <- linnerud$physiological

# Scale data
X <- scale(X_data)
Y <- scale(Y_data)

# test <- spls(X, Y, scale = TRUE, ncomp = n_comp, mode = "regression")


#-------------------- MANUAL sPLS (penalising number of variables) --------------------#

# Convert data into matrices
X_h <- as.matrix(X)
Y_h <- as.matrix(Y)

# Select number of components
n_comp <- 2

# Define keepX/keepY (i.e. number of variables to keep in each component)
keepX = rep(ncol(X), n_comp)
keepY = rep(ncol(Y), n_comp)
# keepX = rep(ncol(X), n_comp)
# keepY = c(1,2)

# Number of variables to penalise (based on keepX/keepY)
x_sparsity <- dim(X_h)[2] - keepX
y_sparsity <- dim(Y_h)[2] - keepY

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
# h<-1
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
  max_iter <- 500
  u_diff <- length(X_h)
  v_diff <- length(X_h)
  
  # Loop until convergence of u and v or max iteration
  while ((sum(abs(u_diff) > tol) != 0 || sum(abs(v_diff) > tol) != 0) && counter < max_iter + 1) {
    
    # Convert number of penalised variables in X into X sparsity parameter for each component
    if (x_sparsity[h] == 0) {
      lambda_u <- 0
    } else {
      lambda_u <- sort(abs(M %*% matrix(v_old, ncol = 1)))[x_sparsity[h]]
    }
    
    # Optimise u iteratively using soft-tresholding penalisation function (and normalise)
    u_new = sign(M %*% v_old) * pmax(0, (abs(M %*% v_old) - lambda_u))
    u_new = u_new / sqrt(as.vector(t(u_new) %*% u_new))
    
    
    # Convert number of penalised variables in Y into Y sparsity parameter for each component
    if (y_sparsity[h] == 0) {
      lambda_v <- 0
    } else {
      lambda_v <- sort(abs(t(M) %*% matrix(u_new, ncol = 1)))[y_sparsity[h]]
    }
    
    # Optimise v iteratively using soft-tresholding penalisation function (and normalise)
    v_new = sign(t(M) %*% u_new) * pmax(0, (abs(t(M) %*% u_new) - lambda_v))
    v_new = v_new / sqrt(as.vector(t(v_new) %*% v_new))
    
    # Update iterables
    u_diff = u_new - u_old
    v_diff = v_new - v_old
    u_old = u_new
    v_old = v_new
    counter = counter + 1
    
    
    # NOTE: Conversion of sparsity parameter to number of variables is made by calculating the part of the 
    # soft thresholding function which corresponds to the magnitude of the new iteration
    # (i.e. abs(M %*% v_old)     for u;    abs(t(M) %*% u_new)     for v)
    # and then choosing lambda to be the value with the smallest weight to ensure that the variables are
    # penalised according to their contribution to the loadings vectors.
    
  }
  
  if (counter == max_iter + 1) {
    print(paste0('Warning! Max iteration reached. Component ', h))
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





#############
test <- spls(X, Y, ncomp = 2, scale = FALSE, keepX = keepX, keepY = keepY)
