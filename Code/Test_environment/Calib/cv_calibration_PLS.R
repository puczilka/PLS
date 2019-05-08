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
#------------------------ CROSS-VALIDATION ------------------------#


# Define Q2 function
Q2 = function(model, X_test, Y_test, n_components) {
  
  vec <- rep(0, n_components)     # Vector to store Q^2 values for each component
  model_prediction <- predict(model, X_test)     # Create object for prediction based on trained model
  
  for (component in 1:n_components) {
    vec[component] = 1 - sum((Y_test - model_prediction$predict[ , , component])^2)/sum((Y_test - mean(Y_test))^2)
    
  }
  
  return(vec) 
}

##### CREATE FOLDS OF DATA

# Define number of folds
no_of_folds <- 5

# Define max number of components (s.t. latent variable matrices are invertible)
ncomp <- min(ncol(X), nrow(X))

# Define number of folds
folds <- sample(cut(seq(1, nrow(X)), breaks = no_of_folds, labels = 1:no_of_folds), size=nrow(X))

# Initiate matrix to store Q^2
col_names <- paste0('comp ', seq(ncomp))
row_names <- paste0('fold ', seq(no_of_folds))
Q2_matrix <- matrix(NA, nrow = no_of_folds, ncol = ncomp, dimnames = list(row_names, col_names))

##### PERFORM CROSS VALIDATION

# Loop over each fold
for (fold in 1:no_of_folds) {
  
  # Assign a fold to be validation indices
  val_indices <- which(folds == fold, arr.ind = TRUE)
  
  # Define training and validation sets
  X_val <- X[val_indices, , drop = FALSE]
  Y_val <- Y[val_indices, , drop = FALSE]
  X_train <- X[-val_indices, , drop = FALSE]
  Y_train <- Y[-val_indices, , drop = FALSE]
  
  # Train PLS regression model and predict on validation set
  pls_model <- pls(X_train, Y_train, ncomp = ncomp, scale = FALSE, mode = 'regression')
  
  # Store Q^2 values in matrix
  Q2_matrix[fold, ] <- Q2(model = pls_model, X_test = X_val, Y_test = Y_val, n_components = ncomp)
  
}
