# Load mixOMICs and linnerud dataset
library('mixOmics')

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological

#-------- SUBSETTING --------#
# Take subset of Y so X and Y have different dimensions (easier for analysing algorithm)
Y <- Y[,1:2]
#----------------------------#


#### POSSIBILITY OF PLS1 AND PLS2 IN ONE FUNCTION?


#-------------------- MANUAL PLS2 --------------------#

# Select number of components
n_comp <- 2

# Scale data
X_h <- as.matrix(scale(X))
Y_h <- as.matrix(scale(Y))

# Create items to store results
W_mat <- NULL     # X weights
T_mat <- NULL     # X scores
U_mat <- NULL     # Y scores
P_mat <- NULL     # X loadings
Q_mat <- NULL     # Y loadings/weights
b <- NULL     # Regression coefficients

###### LOOP OVER COMPONENTS ######
# Initiate first component and loop over defined number of components
h <- 1
while (h <= n_comp) {
  
  ###### CALIBRATION OF COMPONENT ######
  # Initiate u as 1st column of Y
  u <- as.vector(Y_h[,1])
  
  # Initiate iterables
  counter <- 0
  tol <- 1e-06
  max_iter <- 100
  t_diff <- length(X_h)
  old_t <- 0
  
  # Loop until convergence or max iteration
  while (sum(abs(t_diff) > tol) != 0 && counter < max_iter + 1) {
    
    # (Normalised) X weights calculated using Y scores (weights give single representation 
    # for each parameter)
    w = as.vector(u %*% X_h) / as.vector(t(u) %*% u)
    w = w / sqrt(as.vector(t(w) %*% w))
    
    # X scores given by regression of X on weights
    t =  as.vector(X_h %*% w)
    # t =  t / as.vector(t(w) %*% w)
    
    
    # (Normalised) Y weights/loadings calculated using X scores
    q = as.vector(t %*% Y_h) / as.vector(t(t) %*% t)
    q = q / sqrt(as.vector(t(q) %*% q))
    
    # Update Y score by regressing against Y weights/loadings
    u =  as.vector(Y_h %*% q) / as.vector(t(q) %*% q)
    # u =  u / as.vector(t(q) %*% q)
    
    
    # NOTES:
    
    # NORMALISING SCORE
    # The calculations of the scores (for X and Y) are based on the projection of the
    # respective data on the UNIT vector of weights (hence why u,t is divided by the norm of the
    # associated weight). Some algorithms do not use this rescaling as we only care about the
    # RELATIVE magnitudes of the individual scores. The directions are therefore preserved.
    # i.e.
    # t =  as.vector(X_h %*% w)    and      u =  as.vector(Y_h %*% q)
    
    # DISTINCTION BETWEEN LOADINGS/WEIGHTS
    # Some confusion occurs when defining loadings and weights. Loadings are very different from
    # weights for X but not much distinction between loadings and weights is made for Y.
    # This is because loadings in X have the propoerty of orthogonality whilst weights are not
    # necessarily orthogonal in general. Since we want to explain Y in terms of X, a suitable
    # regression (of the scores of X on the scores of Y) would require orthogonality in the X
    # component to remove shared variance among predictors.
    
    
    # Update iterables
    t_diff = t - old_t
    old_t = t
    counter = counter + 1
    
  }
  
  # Calculate X loadings and rescale scores and weights
  p_old = as.vector(t %*% X_h) / as.vector(t(t) %*% t)
  p = p_old     # Remove if including rescaling
  # p = p_old / sqrt(as.vector(t(p_old) %*% p_old))
  # t = t * sqrt(as.vector(t(p_old) %*% p_old))
  # w = w * sqrt(as.vector(t(p_old) %*% p_old))

  
  # NOTE: This rescaling gives a vector of loadings with a norm of 1. Rescaling our associated weights
  # and scores preserves the ratio since this scaling factor cancels out iin the deflation step.
  # Some algorithms therefore do not use this rescaling
  
  
  ###### DEFLATION STEP ######
  # Calculate regression coefficient
  b_h = as.vector(t(u) %*% t) / as.vector(t(t) %*% t)
  
  # NOTE:
  # 'Regression mode' assumes model is of form 'Y ~ X'. Therefore, a prediction of Y based on X
  # can be estimated by regressing the component of X on the coefficients, b_h
  
  # Deflate X and Y matrices
  X_h <- X_h - (t %*% t(p))
  Y_h <- Y_h - b_h * (t %*% t(q))
  
  # Store scores and loadings
  W_mat <- rbind(W_mat, w, deparse.level = 0)
  T_mat <- cbind(T_mat, t, deparse.level = 0)
  U_mat <- cbind(U_mat, u, deparse.level = 0)
  P_mat <- rbind(P_mat, p, deparse.level = 0)
  Q_mat <- rbind(Q_mat, q, deparse.level = 0)
  b <- cbind(b,b_h,  deparse.level = 0)
  
  # Update component number
  h <- h + 1
  
}

# Transpose matrices for consistency with literature
W_mat <- t(W_mat)
P_mat <- t(P_mat)
Q_mat <- t(Q_mat)

# Assign row/column names labels
colnames(T_mat) <- colnames(U_mat) <- colnames(P_mat) <- colnames(Q_mat) <- colnames(W_mat)<- paste0('comp ', seq(n_comp))
rownames(P_mat) <- rownames(W_mat) <- colnames(X)
rownames(Q_mat) <- colnames(Y)


# # ----------------- Visualise regression of scores -------------------- #
# plot(t, u, main = paste('X_scores (t) vs. Y_scores (u): Correlation = ', 
#                         cor(t, u, method = "pearson")), xlab = 'X_scores', ylab = 'Y_scores')
# # abline(lm(u ~ t), col = "red")
# # Should see correlation between X scores and Y scores increasing until convergence
# # -------------------------------------------------------------------- #

# linn_pls <- pls(X, Y, ncomp = 2, mode = "regression")

