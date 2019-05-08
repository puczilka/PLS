library(mixOmics)

# Get some test data
dataset<-as.matrix(data(linnerud))
X_data <- as.matrix(as.data.frame(linnerud$exercise))
Y_data <- as.matrix(as.data.frame(linnerud$physiological))

# Manually scale data
X <- scale(X_data)
Y <- scale(Y_data)


# test <- pls(X, Y, ncomp = 3, scale=FALSE, mode='regression')
# test_pred <- predict(test, X)

#------------------------ FUNCTIONS ------------------------#

# Define MSEP function
X_variable_calibration = function(X_train, Y_train, X_test, Y_test, n_components = min(ncol(X), nrow(X))) {
  
  # Define the total number of variables in the data
  total_vars = ncol(X_test)
  
  # Create objects to store values
  col_names <- paste0('comp ', seq(n_components))
  row_names <- paste0('Number of variables = ', seq(total_vars))
  X.MSEP_matrix <- matrix(NA, nrow = total_vars, ncol = n_components, dimnames = list(row_names, col_names))
  keep_X <- NULL
  tmp <- NULL
  
  for (component in (1:n_components)) {
    
    # Create vector to store MSEPs for each component
    MSEP.vec <- rep(0, total_vars)
    
    for (variable in (1:total_vars)) {
      
      # Initiate objects (or clear existing objects) to store variables
      tmp <- c(keep_X, variable)     # Candidate for keepX argument
      
      model <- spls(X_train, Y_train, ncomp = component, keepX = tmp, scale = FALSE, mode = "regression")
      model_pred <- predict(model, X_test)
      MSEP.vec[variable] <- mean((Y_test -  model_pred$predict[ , , component])^2)
      
    }
    
    # Store the number of variables with the lowest MSEP and update tmp vector.
    # Ensures only 1 entry is returned even if min is not unique
    X.MSEP_matrix[, component] <- MSEP.vec
    keep_X[component] <- which(MSEP.vec == min(MSEP.vec), arr.ind = TRUE)[1]
    
    # Return message if minimum MSEP is not unique
    if (length(which(MSEP.vec == min(MSEP.vec), arr.ind = TRUE)) != 1) {
      
      paste0('Minimum is NOT unique : (component, variable) = (',component, ', ', variable, ')') 
    }
    
  }
  result <- list(X.MSEP_matrix, keep_X)
  names(result) <- c('MSEP', 'keepX')
  
  return(result)
  
}


# Define Q2 function
Q2 = function(model, X_test, Y_test, n_components) {
  
  vec <- rep(0, n_components)     # Vector to store Q^2 values for each component
  model_prediction <- predict(model, X_test)     # Create object for prediction based on trained model
  
  for (component in 1:n_components) {
    vec[component] = 1 - sum((Y_test - model_prediction$predict[ , , component])^2)/sum((Y_test - mean(Y_test))^2)
    
  }
  
  return(vec) 
}


#------------------------ CROSS-VALIDATION ------------------------#

##### SPARSITY ON X

# Define number of folds
no_of_folds <- 5

# Define max number of components (s.t. latent variable matrices are invertible)
ncomp <- min(ncol(X), nrow(X))

# Define number of folds
folds <- sample(cut(seq(1, nrow(X)), breaks = no_of_folds, labels = 1:no_of_folds), size=nrow(X))

# Set max number of components (such that the latent variable matrices are invertible)
max_comp <- min(ncol(X), nrow(X))
total_X_vars <- ncol(X)

# Initiate matrix to store Q^2 values
col_names <- paste0('comp ', seq(max_comp))
row_names <- paste0('fold ', seq(no_of_folds))

Q2_matrix <- matrix(NA, nrow = no_of_folds, ncol = max_comp, dimnames = list(row_names, col_names))


##### PERFORM CROSS VALIDATION

# NOTE: 'Validation' set is sometimes referred to as a 'Test' set but these are slightly different
# in theory (Test data is unseen by model and performance on test data indicates generalisability.
# Validation data is a segment of data used to test model performance during parameter tuning in 
# cross-validation)

# Loop over each fold
for (fold in 1:no_of_folds) {
  
  # Assign a fold to be validation indices
  val_indices <- which(folds == fold, arr.ind = TRUE)
  
  # Define training and validation sets
  X_val <- X[val_indices, , drop = FALSE]
  Y_val <- Y[val_indices, , drop = FALSE]
  X_train <- X[-val_indices, , drop = FALSE]
  Y_train <- Y[-val_indices, , drop = FALSE]
  
  # Run calibration on number of variables to extract calibrated keepX vector
  X_calib <- X_variable_calibration(X_train, Y_train, X_val, Y_val)
  
  # Train sPLS regression model on optimal number of variables and predict on validation set
  pls_model <- pls(X_train, Y_train, ncomp = ncomp, keepX = X_calib$keepX, scale = FALSE, mode = 'regression')
  
  # Store Q^2 values in matrix
  Q2_matrix[fold, ] <- Q2(model = pls_model, X_test = X_val, Y_test = Y_val, n_components = ncomp)
  
}