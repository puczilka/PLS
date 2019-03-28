library(mixOmics)

# Get some test data
dataset<-as.matrix(data(linnerud))
X <- as.matrix(as.data.frame(linnerud$exercise))
Y <- as.matrix(as.data.frame(linnerud$physiological))

## Create folds
kfolds <- 10
folds <- sample(cut(seq(1,nrow(X)),breaks=kfolds,labels=1:kfolds),size=nrow(X))
ncomp<-c(1:min(ncol(X),nrow(X)))
Q2.matrix<-NULL
Q2.total<-c()
MSEP.total<-c()
#Perform 10 fold cross validation
for(k in 1:kfolds) {
  #Segement your data by fold using the which() function 
  testIndeces <- which(folds==k,arr.ind=TRUE)
  x_train <- X[-testIndeces,, drop=FALSE]
  y_train <- Y[-testIndeces,, drop=FALSE]
  x_test <- X[testIndeces,,drop=FALSE]
  y_test <- Y[testIndeces,, drop=FALSE]
  model<-pls(x_train,y_train, ncomp=max(ncomp))
  pred <- predict(model,x_test)
  MSEP<-seq(1:max(ncomp))
  Q2<-seq(1:max(ncomp))
  for (j in 1:max(ncomp)) {
    MSEP[j]<-mean((y_test-pred$predict[,,j])^2) 
    Q2[j]<-1-sum((y_test-pred$predict[,,j])^2)/(sum((y_test-mean(y_test))^2))
  } 
  
  Q2.matrix<-rbind(Q2.matrix,Q2)
  MSEP.matrix<-rbind(MSEP.matrix, MSEP)
  
}
Q2.total<-apply(Q2.matrix, 2, mean)
MSEP.total<-apply(MSEP.matrix, 2, mean)



modelPLS <- pls(X,Y, ncomp=1)
## Run PLS on all the folds
for (j in 1:i) {
  pls_train <- pls(trainData$X,trainData$Y, ncomp=1, scale=TRUE,mode='regression') #(fit on the train data)
  pls_test <- predict(pls_train,newdata=testData) #Get the predicitons for the validation set (from the model just fit on the train data)
}
# For each fold
for(j in 1:k){
  # Fit the model with each subset of predictors on the training part of the fold
  best.fit=regsubsets(Salary~.,data=Hitters[folds!=j,], nvmax=19) 
  # For each subset
  for(i in 1:19){
    # Predict on the hold out part of the fold for that subset
    pred=predict(best.fit, Hitters[folds==j,],id=i)
    # Get the mean squared error for the model trained on the fold with the subset
    cv.errors[j,i]=mean((Hitters$Salary[folds==j]-pred)^2)
  }
}


######################

## Predict values on test set
predict(test)
## compute RMSEP 
MSEP<-mean((Y-y_pred)^2/nrow(Y))
## take the average of MSEP

### Go over all number of components and return MSEP values for every comp

fitted_models <- lapply(ncomp,FUN=function(k) { model<-pls(x_train,y_train, ncomp=k)
pred <- predict(model,x_test)
MSEP<-seq(1:k)
for (j in 1:k) { MSEP[j]<-mean((y_test-pred$predict[,,j])^2)} 
return(MSEP)})

