#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 17:53:40 2019

@author: raphaelsinclair
"""


import numpy as np
import pandas as pd





# Load dataframes
path = '/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH029 - Machine Learning/Project'
df = pd.read_csv(path+'/Combat_filtered_exprs.txt',  sep = '\t').T
df_clin_info = pd.read_csv(path+'/Clinicalinformation.txt', index_col = 0, sep = '\t')

##### DATA PREPROCESSING
# Rename inconvenient columns
df_clin_info.columns = ['Dataset', 'Disease Status', 'Gender', 'Race', 'Age', 'Stage', 'Histology', 'Overall survival (month)', 'Death', 'Smoking', 'TNM stage (T)', 'TNM stage (N)', 'TNM stage (M)', 'Recurrence', 'Others']

# Recode 'Healthy' samples with missing subtype
df_clin_info.loc[((df_clin_info['Histology'].isnull()) & (df_clin_info['Disease Status'] == 'Normal')) ,'Histology'] = 'Healthy' # Replace 'false' missing values

# Recode subtypes with low number of samples as missing
df_clin_info['Histology 2'] = df_clin_info['Histology']
low_sample_subtypes = ['Large cell Neuroendocrine carcinoma', 'Adenosquamous carcinoma', 'Other', 'NSCLC-favor adenocarcinoma', 'NSClarge cell carcinoma-mixed']     # list of subtypes to remove
for subtype in low_sample_subtypes:
    df_clin_info.loc[(df_clin_info['Histology 2'] == subtype), 'Histology 2'] = np.NaN

# Recode disease status labels (0 = Normal, 1 = NSCLC)
df_clin_info['Disease status labels'] = np.where(df_clin_info['Disease Status'] == 'NSCLC', 1, 0)

# Create merged dataframe with covariates and labels. Then create data arrays
clin_info_columns = ['Disease Status', 'Disease status labels', 'Histology', 'Histology 2']     # Create list of column headers from df_clin_info that we want to merge
df_merged = df.merge(df_clin_info[clin_info_columns], how = 'inner', left_index = True, right_index = True)

# Create arrays for covariates and disease status labels
Y = df_merged['Disease status labels'].values
X = df_merged.drop(['Disease Status', 'Disease status labels', 'Histology', 'Histology 2'], axis = 1).values



#-------------------- MANUAL PLS2 --------------------#

# Select number of components
no_of_components = 2


# Scale data (ignoring missing values)
X_h = (X - np.nanmean(X))/np.nanstd(X)
Y_h = (Y - np.nanmean(Y))/np.nanstd(Y)


# Find number of variables and observations for X and Y
if len(X.shape) == 2:
    No_of_X_vars = X.shape[1]
elif len(X.shape) == 1:
    No_of_X_vars = 1
else:
    print('Check dimensions of X')

if len(Y.shape) == 2:
    No_of_Y_vars = Y.shape[1]
elif len(Y.shape) == 1:
    No_of_Y_vars = 1
else:
    print('Check dimensions of Y')
    
if X.shape[0] == Y.shape[0]:
    No_of_obs = X.shape[0]
else:
    print('Number of observations in X and Y do not match')
    
    
# Initiate matrices to store results
W = np.zeros((No_of_X_vars, no_of_components))     # X weights
C = np.zeros((No_of_Y_vars, no_of_components))     # Y weights
T = np.zeros((No_of_obs, no_of_components))        # X scores
U = np.zeros((No_of_obs, no_of_components))        # Y scores
P = np.zeros((No_of_X_vars, no_of_components))     # X loadings
Q = np.zeros((No_of_Y_vars, no_of_components))     # Y loadings
B = np.zeros((no_of_components))                   # 'Regression mode' coefficients


###### LOOP OVER COMPONENTS ######
# Initiate first component and loop over defined number of components
h = 1

while h <= no_of_components:
  
  ###### CALIBRATION OF COMPONENT ######
  # Initiate u as 1st column of Y
  u = Y_h[:,0]
  
  # Initiate iterables
  counter = 0
  tol = 1e-06
  max_iter = 100
  t_diff = X_h.shape[0]
  old_t = 0
  
  # Loop until convergence or max iteration
  while (len(t_diff[np.absolute(t_diff) > tol]) != 0 and counter < max_iter + 1) {
    
    # (Normalised) X weights calculated using Y scores (weights give single representation 
    # for each parameter)
    w = as.vector(u %*% X_h) / as.vector(t(u) %*% u)
    w = w / sqrt(as.vector(t(w) %*% w))
    
    # X scores given by regression of X on weights
    t =  as.vector(X_h %*% w)
    # t =  t / as.vector(t(w) %*% w)
    
    
    # (Normalised) Y weights/loadings calculated using X scores
    c = as.vector(t %*% Y_h) / as.vector(t(t) %*% t)
    c = c / sqrt(as.vector(t(c) %*% c))
    
    # Update Y score by regressing against Y weights/loadings
    u =  as.vector(Y_h %*% c)
    # u =  u / as.vector(t(c) %*% c)
    
    
    # NOTES:
    
    # NORMALISING SCORE
    # The calculations of the scores (for X and Y) are based on the projection of the
    # respective data on the UNIT vector of weights (hence why u,t is divided by the norm of the
    # associated weight). Some algorithms do not use this rescaling as we only care about the
    # RELATIVE magnitudes of the individual scores. The directions are therefore preserved.
    # i.e.
    # t =  as.vector(X_h %*% w)    and      u =  as.vector(Y_h %*% c)
    
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
  
  # Calculate loadings
  p_old = as.vector(t %*% X_h) / as.vector(t(t) %*% t)
  p = p_old
  
  q_old = as.vector(u %*% Y_h) / as.vector(t(u) %*% u)
  q = q_old
  
  
  # # Rescale scores and weights
  # p = p_old / sqrt(as.vector(t(p_old) %*% p_old))
  # t = t * sqrt(as.vector(t(p_old) %*% p_old))
  # w = w * sqrt(as.vector(t(p_old) %*% p_old))

  
  # NOTE: This rescaling gives a vector of loadings with a norm of 1. Rescaling our associated weights
  # and scores preserves the ratio since this scaling factor cancels out iin the deflation step.
  # Some algorithms therefore do not use this rescaling
  
  
  ###### DEFLATION STEP ######
  # Calculate regression coefficient
  b = as.vector(t(u) %*% t) / as.vector(t(t) %*% t)
  
  # NOTE:
  # 'Regression mode' assumes model is of form 'Y ~ X'. Therefore, a prediction of Y based on X
  # can be estimated by regressing the component of X on the coefficients, b_h
  
  # Deflate matrices
  X_h <- X_h - (t %*% t(p))
  Y_h <- Y_h - b * (t %*% t(c))     # Regression
  # Y_h <- Y_h - b_h * (u %*% t(q))     # Canonical
  
  
  # Store scores and loadings
  W_mat <- rbind(W_mat, w, deparse.level = 0)
  C_mat <- rbind(C_mat, c, deparse.level = 0)
  T_mat <- cbind(T_mat, t, deparse.level = 0)
  U_mat <- cbind(U_mat, u, deparse.level = 0)
  P_mat <- rbind(P_mat, p, deparse.level = 0)
  Q_mat <- rbind(Q_mat, q, deparse.level = 0)
  B <- cbind(B, b,  deparse.level = 0)
  
  # Update component number
  h <- h + 1
  
}

# Transpose matrices for consistency with literature
W_mat <- t(W_mat)
C_mat <- t(C_mat)
P_mat <- t(P_mat)
Q_mat <- t(Q_mat)

# Assign row/column names labels
colnames(T_mat) <- colnames(U_mat) <- colnames(P_mat) <- colnames(Q_mat) <- colnames(W_mat) <- colnames(C_mat) <- paste0('comp ', seq(n_comp))
rownames(P_mat) <- rownames(W_mat) <- colnames(X)
rownames(Q_mat) <- rownames(C_mat) <- colnames(Y)


# # ----------------- Visualise regression of scores -------------------- #
# plot(t, u, main = paste('X_scores (t) vs. Y_scores (u): Correlation = ', 
#                         cor(t, u, method = "pearson")), xlab = 'X_scores', ylab = 'Y_scores')
# # abline(lm(u ~ t), col = "red")
# # Should see correlation between X scores and Y scores increasing until convergence
# # -------------------------------------------------------------------- #

# linn_pls <- pls(X, Y, ncomp = 2, mode = "regression")

