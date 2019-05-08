####### gPLS
# library(mixOmics)
library(sgPLS)
#---------------------------- DATA ---------------------------------#
# Load data
setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH030 - Translational Data Sciences/Project/Data')

X <- read.csv('Xdata_simulated.csv')
Y <- read.csv('Ydata_simulated.csv')


# X <- scale(X)
# Y <- scale(Y)

X.s <- scale(as.matrix(X))
Y.s <- scale(as.matrix(Y))


#-------------------- MANUAL gPLS (penalising groups of vaiables) --------------------#

gPLSreg <- function(X, Y, n_components, keepX_group, keepY_group = NULL, ind.block.x, ind.block.y = NULL, tol = 1e-06, max_iter = 100) {
  
  # X = Input data with predictors as columns and observations/samples as rows. This is coerced into a matrix
  # Y = Output data with outcomes as columns and observations/samples as rows. This is coerced into a matrix
  # n_components = The number of components considered for the PLS regression algorithm
  # keepX_group = A vector of length n_components which enforeces sparsity on X. The hth entry corresponds to how many groups to keep for the hth component.
  # keepY_group = A vector of length n_components which enforeces sparsity on Y. The hth entry corresponds to how many groups to keep for the hth component. Default is NULL (i.e. no sparsity on Y).
  # ind.block.x = A vector of column indices denoting the split points (inclusive) of variables for each group in X.
  # ind.block.y = A vector of column indices denoting the split points (inclusive) of variables for each group in Y.
  # tol = The tolerance set for the condition of convergence in the iterative step. Default is 10^-6.
  # max_iter = The maximum number of iterations for the iterative process to run. Default is 100 iterations.
  
  # OPTIONAL
  # trunc_adjust = Boolean to enable adjustment for truncation error. TRUE gives a harder threshold within the soft-thresholding function. FALSE gives the original soft-thresholding function. Default is TRUE.
  
  
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
  
  if (length(keepX_group) != n_components) {
    stop('Length of keepX does not match number of components')
  }
  
  if (!is.null(keepY_group) && length(keepY_group) != n_components) {
    stop('Length of keepY does not match number of components')
  }
  
  if (sum(keepX_group > (length(ind.block.x) + 1)) != 0) {
    stop('keepX_group exceeds the number of groups in X')
  }
  
  if (sum(keepY_group > (length(ind.block.y) + 1)) != 0) {
    stop('keepY_group exceeds the number of groups in Y')
  }
  
  
  
  #==================== Initiate items ====================#
  # Carry out algorithm on X_h, Y_h matices
  X_h <- as.matrix(X)
  Y_h <- as.matrix(Y)
  
  # Dimensions of data
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  k = length(ind.block.x) + 1
  l = length(ind.block.y) + 1
  
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
  
  # Number of groups to penalise
  x_sparsity <- rep(k, n_components) - keepX_group
  
  if (is.null(ind.block.y)) {
    y_sparsity <- rep(0, n_components)
  } else {
    y_sparsity <- rep(l, n_components) - keepY_group
  }
  
  
  #==================== Create blocks of X and Y data ====================#
  # Calculate group indices from ind.block.x/ind.block.y (can more intuitively get group information from this vector
  # instead of ind.block.x/ind.block.y)
  ind.x <- c(0,ind.block.x,ncol(X_h))
  ind.y <- c(0,ind.block.y,ncol(Y_h))
  
  # p_k/q_l holds the number of variables in each group for X and Y respectively
  # x.blocks/y.blocks holds indices for the groups of variables in X and Y respectively
  x.blocks <- list()
  p_k <- NULL
  for (index in 1:k) {
    p_k[index] = ind.x[index + 1] - ind.x[index]
    x.blocks[[index]] = c((ind.x[index] + 1):ind.x[index + 1])
  }
  
  y.blocks <- list()
  q_l <- NULL
  for (index in 1:l) {
    q_l[index] = ind.y[index + 1] - ind.y[index]
    y.blocks[[index]] = ((ind.y[index] + 1):ind.y[index + 1])
  }
  
  # # Set penalty tolerance (i.e. heavyside function). 'trunc_adjust' = TRUE gives tolerance
  # if (trunc_adjust == TRUE) {
  #   pen_tol <- .Machine$double.eps ^ 0.5
  # } else {
  #   pen_tol <- 0
  # }
  
  pen_tol <- 0
  
  
  #==================== Loop over components ====================#
  # Initiate first component and loop over defined number of components
  for (h in (1:n_components)) {
    
    # Compute matrix M (p x q matrix)
    M <- t(X_h) %*% Y_h
    
    #==================== Tune component ====================#
    # WARNING: Truncation errors may occur from long format numbers in calculations. These can be carried over during iterative steps.
    
    # NOTE: Unlike sPLS, we penalise the number of groups involved in producing the component. A similar procedure occurs
    # when calculating the penalties although the group penalties are dependent on the size of the group it belongs to
    # so this is calculated for each block structure.
    
    # Find SVD of M and find loadings from first pair of singular vectors
    M_decomp <- svd(M, nu = 1, nv = 1)
    u_old <- M_decomp$u
    v_old <- M_decomp$v
    
    # Initiate iterables
    counter <- 0
    u_diff <- dim(X_h)[2]
    v_diff <- dim(X_h)[2]
    
    # Loop until convergence of u and v or max iteration
    while ((sum(abs(u_diff) > tol) != 0 || sum(abs(v_diff) > tol) != 0) && counter < max_iter + 1) {
      
      # Calculate the projection of v on M to produce the X loadings candidate
      M_v <- M %*% v_old
      
      # Calculate group lasso penalties using the entries correspeonding to each group from the projection vector
      # and add to x.penalty vector
      x.penalties <- NULL
      for (group in 1:k) {
        vec <- M_v[(x.blocks[[group]])]
        x.penalties <- c(x.penalties, 2*sqrt(sum(vec^2))/sqrt(p_k[group]))
      }
      
      # Convert number of penalised groups in X into sparsity parameter based on group lasso penalties
      if (x_sparsity[h] == 0) {
        lambda_x <- 0
      } else {
        lambda_x <- sort(x.penalties)[x_sparsity[h]]
      }
      
      # Optimise u iteratively for each group (and normalise)
      tmp <- NULL
      for (group in 1:k) {
        vec <- M_v[(x.blocks[[group]])]
        pen <- 1 - (lambda_x/x.penalties[group])
        if (pen < pen_tol) {pen <- 0}
        tmp <- c(tmp, pen*vec)
      }
      u_new = tmp / sqrt(sum(tmp^2))
      
      
      # Calculate the projection of u on M to produce the Y loadings candidate
      M_u <- t(M) %*% u_new
      
      # Calculate group lasso penalties using the entries correspeonding to each group from the projection vector
      # and add to y.penalty vector
      y.penalties <- NULL
      for (group in 1:l) {
        vec <- M_u[(y.blocks[[group]])]
        y.penalties <- c(y.penalties, 2*sqrt(sum(vec^2))/sqrt(q_l[group]))
      }
      
      # Convert number of penalised groups in Y into sparsity parameter based on group lasso penalties
      if (y_sparsity[h] == 0) {
        lambda_y <- 0
      } else {
        lambda_y <- sort(y.penalties)[y_sparsity[h]]
      }
      
      # Optimise v iteratively for each group (and normalise)
      tmp <- NULL
      for (group in 1:l) {
        vec <- M_u[(y.blocks[[group]])]
        pen <- 1 - (lambda_y/y.penalties[group])
        if (pen < pen_tol) {pen <- 0}
        tmp <- c(tmp, pen*vec)
      }
      v_new = tmp / sqrt(sum(tmp^2))
      
      
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
  x.block = list(ind = x.blocks, size = p_k)
  y.block = list(ind = y.blocks, size = q_l)
  result <- list(call = cl, n_components = n_components, keepX_group = keepX_group, keepY_group = keepY_group,
                 ind.block.x = ind.block.x, ind.block.y = ind.block.y,
                 data = list(X = X, Y = Y), blocks = list(x.block = x.block, y.block = y.block),
                 scores = list(X.scores = Eta, Y.scores = Omega),
                 loadings = list(X.loadings = U, Y.loadings = V), defl.coefs = list(C = C, D = D), iterations = iter.comp,
                 names = list(sample = row_names, X.columns = x.col_names, Y.columns = y.col_names), tol = tol)
  
  return(invisible(result))
  
}



#=================================================#

n_components <- 2

# Define ind.block.x/ind.block.y (i.e. vector of indices denoting the end of each group inclusive
# e.g. ind.block.x = c(6, 16) <==> 3 groups s.t. group 1 = 1-6, group 2 = 7-16, group 3 = 17-ncol(X))
ind.block.x <- seq(20, 380, 20)
ind.block.y <- seq(20, 480, 20)

# Select keepX_group/keepY_group variables (i.e. number of groups to keep in each component)
# keepY_group = rep(l, n_components)
keepX_group = c(4, 4)
keepY_group = c(4, 4)


#### gPLS model
t0 <- Sys.time()
test <- gPLSreg(X.s,Y.s,n_components = n_components, keepX_group = keepX_group, keepY_group = keepY_group, ind.block.x = ind.block.x, ind.block.y = ind.block.y, trunc_adjust = TRUE)
t1 <- Sys.time()
print(t1 - t0)

t2 <- Sys.time()
model.gPLS <- gPLS(X, Y, ncomp = n_components, mode = "regression", keepX = keepX_group, 
                   keepY = keepY_group, ind.block.x = ind.block.x , ind.block.y = ind.block.y, scale = TRUE)
t3 <- Sys.time()
print(t3 - t2)

test2 <- gPLSreg(X.s,Y.s,n_components = n_components, keepX_group = keepX_group, keepY_group = keepY_group, ind.block.x = ind.block.x, ind.block.y = ind.block.y, trunc_adjust = FALSE)
