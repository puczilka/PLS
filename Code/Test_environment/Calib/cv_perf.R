# Load mixOMICs and linnerud dataset
library('mixOmics')

data(linnerud)
X_data <- linnerud$exercise
Y_data <- linnerud$physiological

# X <- scale(X_data)
# Y <- scale(Y_data)
X <- (X_data)
Y <- (Y_data)

#--------------------------------------- PLS calibration ---------------------------------------#

PLSperf <- function(X, Y, k_folds)

k_folds <- 4

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

if (nrow(X) < 2) {
  stop('Not enough observations')
}

if (k_folds < 2) {
  stop('The number of folds must be at least 2')
}

if (k_folds >= nrow(X)) {
  stop('k_folds must not be larger than the number of observations')
}


#==================== Initiate items ====================#
# Dimensions of data
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
max_components <- min(n, p)

# Initiate matrix to store Q^2 values
col_names <- paste0(seq(max_components), ' component(s)')
row_names <- paste0('fold ', seq(k_folds))
Q2_matrix <- matrix(NA, nrow = k_folds, ncol = max_components, dimnames = list(row_names, col_names))


#==================== Cross validation ====================#
# Create folds
# set.seed(101)
folds <- sample(cut(seq(1:n), breaks = k_folds, labels = 1:k_folds), size = n)

# Loop over each fold
for (fold_id in 1:k_folds) {
  
  # Assign a fold to be validation indices
  val_indices <- which(folds == fold_id, arr.ind = TRUE)
  
  # Define training and validation sets
  X_val <- X[val_indices, , drop = FALSE]
  Y_val <- Y[val_indices, , drop = FALSE]
  X_train <- X[-val_indices, , drop = FALSE]
  Y_train <- Y[-val_indices, , drop = FALSE]
  
  
  #==================== Define model ====================#
  # Define model and get predictions
  model <- pls(X_train, Y_train, ncomp = max_components, scale = TRUE, mode = "regression", )
  model_pred <- predict(model, X_val)
  
  
  #==================== Define function ====================#
  # Define Q2 function
  Q2 = function(object, Y_test, n_components) {
    
    # Vector to store Q^2 values for each component
    vec <- rep(0, n_components)
    
    # Loop over each component
    for (component in 1:n_components) {
      
      # Calculate Q^2 for each component
      Y.hat <- object$predict[ , , component]
      vec[component] = 1 - sum((Y_test - Y.hat)^2)/sum((Y_test - mean(Y_test))^2)
      
    }
    
    return(vec) 
  }
  
  # Store Q2 vectors for each fold in the Q2 matrix
  Q2_matrix[fold_id,] = Q2(model_pred, Y_val, max_components)
  
}


#==================== Results ====================#
# Average results over folds
avg.results <- apply(Q2_matrix, MARGIN = 2, FUN = mean)

# Find optimal values
opt_val <- max(avg.results)
opt_param <- which(avg.results == opt_Q2)
names(opt_param) <- NULL

# Return final outputs
cl = match.call()
result<- list(call = cl, results = Q2_matrix, avg.results = avg.results, optimal.value = opt_val, optimal.param = opt_param)

return(invisible(result))

# }