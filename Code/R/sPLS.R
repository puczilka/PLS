# Load mixOMICs and linnerud dataset
library('mixOmics')

data(linnerud)
X_data <- linnerud$exercise
Y_data <- linnerud$physiological

X <- scale(X_data)
Y <- scale(Y_data)


#--------------------------------------- sPLS (variable penalisation) ---------------------------------------#

sPLSreg <- function(X, Y, n_components, keepX, keepY = NULL, tol = 1e-06, max_iter = 100) {
  
  # X = Input data with predictors as columns and observations/samples as rows. This is coerced into a matrix
  # Y = Output data with outcomes as columns and observations/samples as rows. This is coerced into a matrix
  # n_components = The number of components considered for the PLS regression algorithm
  # keepX = A vector of length n_components which enforeces sparsity on X. The hth entry corresponds to how many variables to keep in the X loadings vector for the hth component.
  # keepY = A vector of length n_components which enforeces sparsity on Y. The hth entry corresponds to how many variables to keep in the Y loadings vector for the hth component. Default is NULL (i.e. no sparsity on Y).
  # tol = The tolerance set for the condition of convergence in the iterative step. Default is 10^-6.
  # max_iter = The maximum number of iterations for the iterative process to run. Default is 100 iterations.
  
  
  #==================== Initial checks ====================#
  # Coerce data into matrices to store original data
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Check data
  if (length(dim(X)) != 2) {
    stop('Check dimensions of X')
  }
  if (length(dim(Y)) > 2) {
    stop('Check dimensions of Y')
  }
  
  if (nrow(X) != nrow(Y)) {
    stop('Number of observations in X and Y do not match')
  }
  
  if (n_components > min(nrow(X), ncol(X))) {
    stop('Exceeded maximum number of components')
  }
  
  if (length(keepX) != n_components) {
    stop('Length of keepX does not match number of components')
  }
  
  if (!is.null(keepY) && length(keepY) != n_components) {
    stop('Length of keepY does not match number of components')
  }
  
  if (sum(keepX > ncol(X)) != 0) {
    stop('keepX exceeds the number of variables in X')
  }
  
  if (sum(keepY > ncol(Y)) != 0) {
    stop('keepY exceeds the number of variables in Y')
  }
  
  #==================== Initiate items ====================#
  # Carry out algorithm on X_h, Y_h matices
  X_h <- as.matrix(X)
  Y_h <- as.matrix(Y)
  
  # Dimensions of data
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # Column/row names
  if (is.null(rownames(X))) {
    if (is.null(rownames(Y))) {
      row_names <- c(1:n)
    } else {
      row_names <- rownames(Y)
    }
  } else {
    row_names <- rownames(X)
  }
  
  if (is.null(colnames(X))) {
    x.col_names <- c(paste0('X', seq(p)))
  } else {
    x.col_names <- colnames(X)
  }
  
  if (is.null(colnames(Y))) {
    y.col_names <- c(paste0('Y', seq(q)))
  } else {
    y.col_names <- colnames(Y)
  }
  
  
  # Create items to store results and assign row/column names
  Eta <- matrix(NA, nrow = n, ncol = n_components, dimnames = list(row_names, paste0('comp ', seq(n_components))))     # X scores 
  Omega <- matrix(NA, nrow = n, ncol = n_components, dimnames = list(row_names, paste0('comp ', seq(n_components))))     # Y scores
  U <- matrix(NA, nrow = p, ncol = n_components, dimnames = list(x.col_names, paste0('comp ', seq(n_components))))     # X loadings
  V <- matrix(NA, nrow = q, ncol = n_components, dimnames = list(y.col_names, paste0('comp ', seq(n_components))))     # Y loadings
  C <- matrix(NA, nrow = p, ncol = n_components, dimnames = list(x.col_names, paste0('comp ', seq(n_components))))     # Regression coefficient on latent variables (for X_h)
  D <- matrix(NA, nrow = q, ncol = n_components, dimnames = list(y.col_names, paste0('comp ', seq(n_components))))     # Regression coefficient on latent variables in 'Regression mode' (for Y_h)
  # E <- matrix(NA, nrow = q, ncol = n_components, dimnames = list(y.col_names, paste0('comp ', seq(n_components))))     # Regression coefficient on latent variables in 'PLS-mode A' (for Y_h)
  
  iter.comp <- matrix(NA, nrow = 1, ncol = n_components, dimnames = list(NULL, paste0('comp ', seq(n_components))))     # Stores iterations for each component
  
  # Number of variables to penalise (based on keepX/keepY)
  x_sparsity <- dim(X_h)[2] - keepX
  
  if (is.null(keepY)) {
    keepY = rep(ncol(Y), n_components)
  }
  y_sparsity <- dim(Y_h)[2] - keepY
  
  
  #==================== Loop over components ====================#
  # Loop over defined number of components
  for (h in (1:n_components)) {
    
    
    #==================== Tune component ====================#
    # WARNING: Truncation errors may occur from long format numbers in calculations. These can be carried over during iterative steps.
    
    # NOTE: Conversion of sparsity parameter to number of variables is made by calculating the magnitude 
    # of each variables contribution to the covariance structure between X and Y.
    # If we chose to penalise k parameters (in X w.l.o.g), we choose lambda to be the kth smallest magnitude.
    # Therefore lambda is chosen s.t. its value is sufficiently large enough to penalise the k entries which 
    # contribute least (in the loadings). Similarly for penalisation on Y
    
    # Compute matrix M (p x q matrix)
    M <- t(X_h) %*% Y_h
    
    # Find SVD of M and extract loadings (first pair of left and right singular vectors)
    M_decomp <- svd(M, nu = 1, nv = 1)
    u_old <- M_decomp$u
    v_old <- M_decomp$v
    
    # Initiate iterables
    counter <- 0
    u_diff <- dim(X_h)[2]
    v_diff <- dim(Y_h)[2]
    
    # Loop until convergence of u and v or max iteration
    while ((sum(abs(u_diff) > tol) != 0 || sum(abs(v_diff) > tol) != 0) && counter < max_iter + 1) {
      
      # Calculate the projection of v on M to produce the X loadings candidate
      M_v <- M %*% v_old
      
      # Convert number of penalised variables in X into X sparsity parameter for each component
      if (x_sparsity[h] == 0) {
        lambda_x <- 0
      } else {
        lambda_x <- sort(abs(M_v))[x_sparsity[h]]
      }
      
      # Optimise u iteratively using soft-tresholding penalisation function (and normalise)
      u_new = sign(M_v) * pmax(0, (abs(M_v) - lambda_x))
      u_new = u_new / sqrt(sum(u_new^2))
      
      
      # Calculate the projection of u on M to produce the Y loadings candidate
      M_u <- t(M) %*% u_new
      
      # Convert number of penalised variables in Y into Y sparsity parameter for each component
      if (y_sparsity[h] == 0) {
        lambda_y <- 0
      } else {
        lambda_y <- sort(abs(M_u))[y_sparsity[h]]
      }
      
      # Optimise v iteratively using soft-tresholding penalisation function (and normalise)
      v_new = sign(M_u) * pmax(0, (abs(M_u) - lambda_y))
      v_new = v_new / sqrt(sum(v_new^2))
      
      
      # Update iterables
      u_diff = u_new - u_old
      v_diff = v_new - v_old
      u_old = u_new
      v_old = v_new
      counter = counter + 1
      
    }
    
    # Check convergence
    if (counter == max_iter + 1) {
      warning(paste0('Warning! Max iteration reached. No convergence for component ', h))
    }
    
    # Add number of iterations to vector
    iter.comp[, h] <- counter
    
    
    #==================== Deflation step ====================#
    # Calculate scores/latent variables for X and Y
    eta = as.vector(X_h %*% u_new) / sum(u_new^2)
    omega = as.vector(Y_h %*% v_new) / sum(v_new^2)
    
    # Calculate regression coefficients
    c = as.vector(t(X_h) %*% eta) / sum(eta^2)
    d = as.vector(t(Y_h) %*% eta) / sum(eta^2)
    # e = as.vector(t(Y_h) %*% omega) / sum(omega^2)
    
    # Deflate X and Y matrices using latent variables and regression coefficients
    X_h <- X_h - (eta %*% t(c))
    Y_h <- Y_h - (eta %*% t(d))
    
    # Store variables
    Eta[, h] <- eta
    Omega[, h] <- omega
    U[, h] <- u_new
    V[, h] <- v_new
    C[, h] <- c
    D[, h] <- d
    # E[, h] <- e
    
  }
  
  
  #==================== Form predictions using results ====================#
  # Create function for prediction
  
  # ??????
  
  #========================================================================#
  
  
  # Return final outputs
  cl = match.call()
  result <- list(call = cl, n_components = n_components, keepX = keepX, keepY = keepY, data = list(X = X, Y = Y),
                 scores = list(X.scores = Eta, Y.scores = Omega), loadings = list(X.loadings = U, Y.loadings = V),
                 defl.coefs = list(C = C, D = D), iterations = iter.comp,
                 names = list(sample = row_names, X.columns = x.col_names, Y.columns = y.col_names), tol = tol)
  
  return(invisible(result))
  
}

# Test function
keepX = c(1,2)
keepY = c(2,1)

test <- sPLSreg(X, Y, n_components = 2, keepX = keepX, keepY = keepY)
mO_test <- spls(X_data, Y_data, ncomp = 2, scale = TRUE, keepX = keepX, keepY = keepY)