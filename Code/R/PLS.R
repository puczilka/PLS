# Load mixOMICs and linnerud dataset
library('mixOmics')

data(linnerud)
X_data <- linnerud$exercise
Y_data <- linnerud$physiological

X <- scale(X_data)
Y <- scale(Y_data)


#--------------------------------------------- PLS ---------------------------------------------#

PLSreg <- function(X, Y, n_components, tol = 1e-06, max_iter = 100) {
  
  # X = Input data with predictors as columns and observations/samples as rows. This is coerced into a matrix.
  # Y = Output data with outcomes as columns and observations/samples as rows. This is coerced into a matrix.
  # n_components = The number of components considered for the PLS regression algorithm.
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
  
  if (length(dim(Y)) == 1) {
    type <- 'PLS1'
  } else {
    type <- 'PLS2'}
  
  if (nrow(X) != nrow(Y)) {
    stop('Number of observations in X and Y do not match')
  }
  
  if (n_components > min(nrow(X), ncol(X))) {
    stop('Exceeded maximum number of components')
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
  W <- matrix(NA, nrow = p, ncol = n_components, dimnames = list(x.col_names, paste0('comp ', seq(n_components))))     # X weights
  C <- matrix(NA, nrow = q, ncol = n_components, dimnames = list(y.col_names, paste0('comp ', seq(n_components))))     # Y weights
  T <- matrix(NA, nrow = n, ncol = n_components, dimnames = list(row_names, paste0('comp ', seq(n_components))))     # X scores
  U <- matrix(NA, nrow = n, ncol = n_components, dimnames = list(row_names, paste0('comp ', seq(n_components))))     # Y scores
  P <- matrix(NA, nrow = p, ncol = n_components, dimnames = list(x.col_names, paste0('comp ', seq(n_components))))     # X loadings
  Q <- matrix(NA, nrow = q, ncol = n_components, dimnames = list(y.col_names, paste0('comp ', seq(n_components))))     # Y loadings
  B <- matrix(NA, nrow = 1, ncol = n_components, dimnames = list(NULL, paste0('comp ', seq(n_components))))     # Regression coefficients
  
  iter.comp <- matrix(NA, nrow = 1, ncol = n_components, dimnames = list(NULL, paste0('comp ', seq(n_components))))     # Stores iterations for each component
  
  
  #==================== Loop over components ====================#
  # Loop over defined number of components
  for (h in (1:n_components)) {
    
    
    #==================== Tune component (PLS1/PLS2) ====================#
    # NOTE: Any commented sections with '(OPTIONAL)' beside them are optional rescailings of quantities which are
    # consistent with the PLS tutorial: Partial Least Squares Regression: A tutorial - Paul Geladi & 
    # Bruce Kowalski (1986)
    
    # PLS1 if Y is univariate, PLS2 if Y is multivariate
    if (type == 'PLS1') {
      
      # Initiate u as the only column of Y
      u <- as.vector(Y_h)
      
      # Update X weights and normalise
      w = as.vector(u %*% X_h) / as.vector(t(u) %*% u)
      w = w / sqrt(sum(w^2))
      
      # Calculate X scores by regressing X on X weights
      t =  as.vector(X_h %*% w)
      # t =  t / as.vector(t(w) %*% w)     # (OPTIONAL)
      
      # Y has one variable and therefore it is assigned a weight of 1. Similar for loadings
      c = 1
      q = 1
      
    } else {
      
      # Initiate u as 1st column of Y
      u <- as.vector(Y_h[,1])
      
      
      # Initiate iterables
      counter <- 0
      t_diff <- dim(X_h)[2]
      old_t <- 0
      
      # Loop until convergence or max iteration (for multivariate Y)
      while (sum(abs(t_diff) > tol) != 0 && counter < max_iter + 1) {
        
        # Update X weights and normalise
        w = as.vector(u %*% X_h) / sum(u^2)
        w = w / sqrt(sum(w^2))
        
        # Calculate X scores by regressing X on X weights
        t =  as.vector(X_h %*% w)
        # t =  t / sum(w^2)     # (OPTIONAL)
        
        # Update Y weights and normalise
        c = as.vector(t %*% Y_h) / sum(t^2)
        c = c / sqrt(sum(c^2))
        
        # Calculate Y scores by regressing Y on Y weights
        u =  as.vector(Y_h %*% c)
        # u =  u / sum(c^2)     # (OPTIONAL)
        
        
        # Update iterables
        t_diff = t - old_t
        old_t = t
        counter = counter + 1
        
      }
      
      # Check convergence
      if (counter == max_iter + 1) {
        warning(paste0('Warning! Max iteration reached. No convergence for component ', h))
      }
      
      # Add number of iterations to vector
      iter.comp[, h] <- counter
      
    }
    
    
    #==================== Deflation step ====================#
    # Calculate X loadings
    p_old = as.vector(t %*% X_h) / sum(t^2)
    p = p_old
    # p = p_old / sqrt(sum(p_old^2))     # (OPTIONAL)
    # t = t * sqrt(sum(p_old^2))     # (OPTIONAL)
    # w = w * sqrt(sum(p_old^2))    # (OPTIONAL)
    
    # Calculate Y loadings
    q_old = as.vector(u %*% Y_h) / sum(u^2)
    q = q_old
    # q = q_old / sqrt(sum(q_old^2))     # (OPTIONAL)
    # u = u * sqrt(sum(q_old^2))     # (OPTIONAL)
    # c = c * sqrt(sum(q_old^2))     # (OPTIONAL)
    
    # Calculate regression coefficient for regression mode
    b = as.vector(t(u) %*% t) / sum(t^2)
    
    # Deflate matrices
    X_h <- X_h - (t %*% t(p))
    Y_h <- Y_h - b * (t %*% t(c))
    
    
    #==================== Store results for each component ====================#
    # Store results
    W[, h] <- w
    C[, h] <- c
    T[, h] <- t
    U[, h] <- u
    P[, h] <- p
    Q[, h] <- q
    B[, h] <- b
    
  }
  
  
  #==================== Form predictions using results ====================#
  # Create function for prediction
  
  # ??????
  
  #========================================================================#
  
  
  # Return final outputs
  cl = match.call()
  result <- list(call = cl, n_components = n_components, data = list(X = X, Y = Y),
                 scores = list(X.scores = T, Y.scores = U), weights = list(X.weights = W, Y.weights = C),
                 loadings = list(X.loadings = P, Y.loadings = Q), regression.coef = B,
                 iterations = iter.comp, algorithm = type,
                 names = list(sample = row_names, X.columns = x.col_names, Y.columns = y.col_names), tol = tol)
  
  return(invisible(result))
  
}

# Test function
test <- PLSreg(X = X, Y = Y, n_components = 2)
# linn_pls <- pls(X, Y, ncomp = 2, mode = "regression")