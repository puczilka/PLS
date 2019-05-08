#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 17:53:40 2019

@author: raphaelsinclair
"""


import numpy as np
import pandas as pd



# Load dataframes (NOTE: linnerud dataset is incorrectly labelled. Error is preserved for consistency)
path = '/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH030 - Translational Data Sciences/Project'

X_data = pd.read_csv(path+'/Data/linnerud_physiological.csv', sep = ' ')
Y_data = pd.read_csv(path+'/Data/linnerud_exercise.csv', sep = ' ')



# Select number of components
no_of_components = 2


# Scale data ignoring missing values (Use unbiased estimator of standard deviation i.e. divisor = (n-1))
X = X_data.values
Y = Y_data.values

X_h = (X - np.nanmean(X, axis = 0))/np.nanstd(X, axis = 0, ddof = 1)
Y_h = (Y - np.nanmean(Y, axis = 0))/np.nanstd(Y, axis = 0, ddof = 1)


# =============================================================================
# def scale(data):
#     return (data.values - np.nanmean(data.values, axis = 0), ddof = 1)/np.nanstd(data.values, axis = 0, ddof = 1)

# X_h = scale(X_data)
# Y_h = scale(Y_data)
# =============================================================================





#-------------------- MANUAL PLS2 --------------------#


# Find number of variables and observations for X and Y
if len(X_h.shape) == 2:
    No_of_X_vars = X_h.shape[1]
elif len(X_h.shape) == 1:
    No_of_X_vars = 1
else:
    print('Check dimensions of X')

if len(Y_h.shape) == 2:
    No_of_Y_vars = Y_h.shape[1]
elif len(Y_h.shape) == 1:
    No_of_Y_vars = 1
else:
    print('Check dimensions of Y')
    
if X_h.shape[0] == Y_h.shape[0]:
    No_of_obs = X.shape[0]
else:
    print('Number of observations in X and Y do not match')
    
    
# Initiate matrices to store results
W = np.zeros((No_of_X_vars, no_of_components))     # X weights
C = np.zeros((No_of_Y_vars, no_of_components))     # Y weights
T = np.zeros((No_of_obs, no_of_components))        # X scores
U = np.zeros((No_of_obs, no_of_components))        # Y scores
P = np.zeros((No_of_X_vars, no_of_components))     # X loadings
Q = np.zeros((No_of_Y_vars, no_of_components))     # Y loadings 'Canonical mode'
B = np.zeros((no_of_components))                   # 'Regression mode' coefficients


###### LOOP OVER COMPONENTS ######
# Loop over defined number of components
for h in range(0,no_of_components):


    ###### CALIBRATION OF COMPONENT ######
    # Initiate u as 1st column of Y
    u = Y_h[:,0]
    
    # Initiate iterables
    counter = 0
    tol = 1e-06
    max_iter = 100
    t_diff = np.ones(X_h.shape[0])
    old_t = 0
    
    
    # Loop until convergence or max iteration
    while (len(t_diff[np.absolute(t_diff) > tol]) != 0 and counter < max_iter + 1):
    
    
        # (Normalised) X weights calculated using Y scores (weights give single representation
        # for each parameter)
        w = X_h.T.dot(u) / u.dot(u)
        w = w / (w.dot(w) ** 0.5)
        
        # X scores given by regression of X on weights (Normalisation of scores is optional)
        t =  X_h.dot(w)
        # t =  t / w.dot(w)
        
        # (Normalised) Y weights calculated using X scores
        c = Y_h.T.dot(t) / t.dot(t)
        c = c / (c.dot(c) ** 0.5)
        
        # Update Y score by regressing against Y weights (Normalisation of scores is optional)
        u =  Y_h.dot(c)
        # u =  u / c.dot(c)
        
        
        # =============================================================================
        #  NOTES:
        # 
        #  NORMALISING SCORE
        #  The calculations of the scores (for X and Y) are based on the projection of the
        #  respective data on the UNIT vector of weights (hence why u,t is divided by the norm of the
        #  associated weight). Some algorithms do not use this rescaling as we only care about the
        #  RELATIVE magnitudes of the individual scores. The directions are therefore preserved.
        #  i.e.
        #  t =  as.vector(X_h %*% w)    and      u =  as.vector(Y_h %*% c)
        # 
        #  DISTINCTION BETWEEN LOADINGS/WEIGHTS
        #  Some confusion occurs when defining loadings and weights. Loadings are very different from
        #  weights for X but not much distinction between loadings and weights is made for Y.
        #  This is because loadings in X have the propoerty of orthogonality whilst weights are not
        #  necessarily orthogonal in general. Since we want to explain Y in terms of X, a suitable
        #  regression (of the scores of X on the scores of Y) would require orthogonality in the X
        #  component to remove shared variance among predictors.
        # 
        # =============================================================================
        
        
        # Update iterables
        t_diff = t - old_t
        old_t = t
        counter = counter + 1
        
          
        # Calculate loadings
        p_old = X_h.T.dot(t) / t.dot(t)
        p = p_old
              
        q_old = Y_h.T.dot(u) / u.dot(u)
        q = q_old
          
          
        # Rescale scores and weights (optional, given that scores are normalised)
        # p = p_old / (p_old.dot(p_old) ** 0.5)
        # t = t * (p_old.dot(p_old) ** 0.5)
        # w = w * (p_old.dot(p_old) ** 0.5)
        
          
        # =============================================================================
        #  NOTE: This rescaling gives a vector of loadings with a norm of 1. Rescaling our associated weights
        #  and scores preserves the ratio since this scaling factor cancels out iin the deflation step.
        #  Some algorithms therefore do not use this rescaling
        # =============================================================================
          
      
    ###### DEFLATION STEP ######
    # Calculate regression coefficient
    b = u.dot(t) / t.dot(t)
      
      
    # Deflate matrices
    X_h = X_h - np.outer(t, p)
    Y_h = Y_h - b * np.outer(t, c)     # Regression
    # Y_h = Y_h - np.outer(u, q)     # Canonical
    
    
    # =============================================================================
    #  NOTE:
    #  'Regression mode' assumes model is of form 'Y ~ X'. Therefore, a prediction of Y based on X
    #  can be estimated by regressing the component of X on the coefficients, b
    # =============================================================================
    
    
    # Store scores and loadings
    W[:,h] = w
    C[:,h] = c
    T[:,h] = t
    U[:,h] = u
    P[:,h] = p
    Q[:,h] = q
    B[h] = b


# Assign row/column labels
x_columns = X_data.columns
y_columns = Y_data.columns
x_rows = X_data.index
y_rows = Y_data.index

column_names = []
for index in range(1, no_of_components + 1):
    column_names.append('comp '+ str(index))


W = pd.DataFrame(W, index = x_columns, columns = column_names)
C = pd.DataFrame(C, index = y_columns, columns = column_names)
T = pd.DataFrame(T, index = x_rows, columns = column_names)
U = pd.DataFrame(U, index = y_rows, columns = column_names)
P = pd.DataFrame(P, index = x_columns, columns = column_names)
Q = pd.DataFrame(Q, index = y_columns, columns = column_names)
B = pd.DataFrame(B.reshape(1, no_of_components), columns = column_names)
