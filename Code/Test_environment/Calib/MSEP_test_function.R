# Define MSEP function
X_variable_calibration = function(X_train, Y_train, X_test, Y_test, n_components = min(ncol(X), nrow(X)), no_of_folds = 5) {
  
  # Define the total number of variables in the data
  total_vars = ncol(X_test)
  
  # # Create objects to store values
  # col_names <- paste0('Number of variables = ', seq(total_vars))
  # row_names <- paste0('fold ', seq(no_of_folds))
  # mat <- matrix(NA, nrow = total_vars, ncol = no_of_folds, dimnames = list(row_names, col_names))
  
  # Loop over components (i.e. Cross validate for each component)
  for (component in (1:n_components)) {
  
  # Create objects to store values
  col_names <- paste0('Number of variables = ', seq(total_vars))
  row_names <- paste0('Fold ', seq(no_of_folds))
  mat <- matrix(NA, nrow = total_vars, ncol = no_of_folds, dimnames = list(row_names, col_names))
  keep_X <- NULL
  tmp <- NULL
  
  # Create vector to store MSEPs for each component
  MSEP.vec <- rep(0, total_vars)
  
    # Loop over each fold
    for (fold in 1:no_of_folds) {
      
      # Assign a fold to be validation indices
      val_indices <- which(folds == fold, arr.ind = TRUE)
      
      # Define training and validation sets
      X_val <- X[val_indices, , drop = FALSE]
      Y_val <- Y[val_indices, , drop = FALSE]
      X_train <- X[-val_indices, , drop = FALSE]
      Y_train <- Y[-val_indices, , drop = FALSE]
      
      
      for (variable in (1:total_vars)) {
        
        # Initiate objects (or clear existing objects) to store variables
        tmp <- c(keep_X, variable)     # Candidate for keepX argument
        
        model <- spls(X_train, Y_train, ncomp = component, keepX = tmp, scale = FALSE, mode = "regression")
        model_pred <- predict(model, X_test)
        
        # Store MSEP for the corresponding number of variables
        mat[variable, fold] <- mean((Y_test -  model_pred$predict[ , , component])^2)
        
      }
      
    }
  
  # Average across all folds
  avg.MSEP <- apply(mat, MARGIN = 2, FUN = mean)
  
  
  
  
  
  
  
  # Ensures only 1 entry is returned even if min is not unique
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


