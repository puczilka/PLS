library(mixOmics)

# Get some test data
data(linnerud)
X <- as.matrix(as.data.frame(linnerud$exercise))
Y <- as.matrix(as.data.frame(linnerud$physiological))

# Create folds
kfolds <- 10

set.seed(101)
folds <- sample(cut(seq(1,nrow(X)), breaks = kfolds, labels = 1:kfolds), size = nrow(X))

# Set max number of components (such that the latent variable matrices are invertible)
max_comp <- min(ncol(X), nrow(X))
ncomp <- c(1:max_comp)
max_X_vars <- ncol(X)
max_Y_vars <- ncol(Y)

# # Initiate objects to store metrics
# Q2.matrix <- NULL
# Q2.total <- c()
MSEP.matrix <- NULL

k <- 1
# Reshuffle folds
# for (k in 1:kfolds) {

# Split data into training/validation sets
val_indices <- which(folds == k, arr.ind=TRUE)

X_val <- X[val_indices, ,drop=FALSE]
Y_val <- Y[val_indices, ,drop=FALSE]

X_train <- X[-val_indices, ,drop=FALSE]
Y_train <- Y[-val_indices, ,drop=FALSE]


# ------------------- TEST MSEP ------------------- #
test <- spls(X_train, Y_train, ncomp = 2, keepX = c(2,3) , mode = "regression")
test_pred <- predict(test, X_val)
test_pred <- test_pred$predict[ , , 2]

MSEP_test <- mean((Y_val - test_pred)^2)

test2 <- spls(X_train, Y_train, ncomp = 2, keepX = c(2,2) , mode = "regression")
test_pred2 <- predict(test2, X_val)
test_pred2 <- test_pred2$predict[ , , 2]

MSEP_test2 <- mean((Y_val - test_pred2)^2)

# ------------------------------------------------- #


# NOTE: 'Validation' set is sometimes referred to as a 'Test' set but these are slightly different
# in theory (Test data is unseen by model and performance on test data indicates generalisability.
# Validation data is a segment of data used to test model performance during parameter tuning in 
# cross-validation)


##### SPARSE ON X ONLY

# Train model using training data and test performance on validation data
keep_X <- NULL
tmp <- NULL
comp <- 1

for (comp in (1:max_comp)) {
  
  # Return the MSEP for each number of variables
  MSEP.vec <- c()

  for (var_no in (1:max_X_vars)) {
    
    # Initiate objects (or clear existing objects) to store variables
    tmp <- c(keep_X, var_no)
    
    model_vars <- spls(X_train, Y_train, ncomp = comp, keepX = tmp, mode = "regression")
    pred <- predict(model_vars, X_val)
    MSEP.vec[var_no] <- mean((Y_val - pred$predict[ , , comp])^2)
    
  }
  
  # Store the number of variables with the lowest MSEP and update tmp vector
  MSEP.matrix <- cbind(MSEP.matrix, MSEP.vec)
  keep_X[comp] <- which(MSEP.vec == min(MSEP.vec), arr.ind = TRUE)
  
  comp <- comp + 1
  
}

#   for (j in 1:max(ncomp)) {
#     MSEP[j]<-mean((y_test-pred$predict[,,j])^2) 
#     Q2[j]<-1-sum((y_test-pred$predict[,,j])^2)/(sum((y_test-mean(y_test))^2))
#   } 
#   
#   Q2.matrix<-rbind(Q2.matrix,Q2)
#   MSEP.matrix<-rbind(MSEP.matrix, MSEP)
#   
# }
# Q2.total<-apply(Q2.matrix, 2, mean)
# MSEP.total<-apply(MSEP.matrix, 2, mean)
# 
# 
# 
# modelPLS <- pls(X,Y, ncomp=1)
# ## Run PLS on all the folds
# for (j in 1:i) {
#   pls_train <- pls(trainData$X,trainData$Y, ncomp=1, scale=TRUE,mode='regression') #(fit on the train data)
#   pls_test <- predict(pls_train,newdata=testData) #Get the predicitons for the validation set (from the model just fit on the train data)
# }
# # For each fold
# for(j in 1:k){
#   # Fit the model with each subset of predictors on the training part of the fold
#   best.fit=regsubsets(Salary~.,data=Hitters[folds!=j,], nvmax=19) 
#   # For each subset
#   for(i in 1:19){
#     # Predict on the hold out part of the fold for that subset
#     pred=predict(best.fit, Hitters[folds==j,],id=i)
#     # Get the mean squared error for the model trained on the fold with the subset
#     cv.errors[j,i]=mean((Hitters$Salary[folds==j]-pred)^2)
#   }
# }
# 
# 
# ######################
# 
# ## Predict values on test set
# predict(test)
# ## compute RMSEP 
# MSEP<-mean((Y-y_pred)^2/nrow(Y))
# ## take the average of MSEP
# 
# ### Go over all number of components and return MSEP values for every comp
# 
# fitted_models <- lapply(ncomp,FUN=function(k) { model<-pls(x_train,y_train, ncomp=k)
# pred <- predict(model,x_test)
# MSEP<-seq(1:k)
# for (j in 1:k) { MSEP[j]<-mean((y_test-pred$predict[,,j])^2)} 
# return(MSEP)})
# 
# 
# MSEP<-seq(1:max(ncomp))
# 
