#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 18:00:54 2019

@author: raphaelsinclair
"""






# Import libraries
import numpy as np
import pandas as pd



# Load dataframes (NOTE: linnerud dataset is incorrectly labelled. Error is preserved for consistency)
path = '/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH030 - Translational Data Sciences/Project'

X_data = pd.read_csv(path+'/Data/linnerud_physiological.csv', sep = ' ')
Y_data = pd.read_csv(path+'/Data/linnerud_exercise.csv', sep = ' ')


# Scale data ignoring missing values (Use unbiased estimator of standard deviation i.e. divisor = (n-1))

X_data = (X_data - np.nanmean(X_data, axis = 0))/np.nanstd(X_data, axis = 0, ddof = 1)
Y_data = (Y_data - np.nanmean(Y_data, axis = 0))/np.nanstd(Y_data, axis = 0, ddof = 1)



#--------------------------------------------- PLS ---------------------------------------------#

# PLSreg = def(X, Y, n_components, tol = 1e-06, max_iter = 100):
  
# X = Input data with predictors as columns and observations/samples as rows. This is coerced into a matrix.
# Y = Output data with outcomes as columns and observations/samples as rows. This is coerced into a matrix.
# n_components = The number of components considered for the PLS regression algorithm.
# tol = The tolerance set for the condition of convergence in the iterative step. Default is 10^-6.
# max_iter = The maximum number of iterations for the iterative process to run. Default is 100 iterations.


#==================== Initial checks ====================#
# Coerce data into matrices to store original data
X_h = X.values
Y_h = Y.values

# Check data
if len(X_h.shape) != 2:
    return('Check dimensions of X')
      
if len(Y_h.shape) > 2:
    return('Check dimensions of Y')
    
if X_h.shape[0] != Y_h.shape[0]: 
    return('Number of observations in X and Y do not match')
  
if n_components > min(X_h.shape[0], X_h.shape[1]):
    return('Exceeded maximum number of components')


#==================== Initiate items ====================#

# Dimensions of data
n = X_h.shape[0]
p = X_h.shape[1]
if len(Y_h.shape) == 1:
    type = 'PLS1'
    q = 1
elif len(Y_h.shape) == 2:
    type = 'PLS2'
    q = Y_h.shape[1]

# Assign row/column labels
x_columns = X.columns
y_columns = Y.columns
id_rows = X.index

comp_no = []
for index in range(1, no_of_components + 1):
    comp_no.append('comp '+ str(index))


# Initiate matrices to store results
W = np.zeros((p, n_components))     # X weights
C = np.zeros((q, n_components))     # Y weights
T = np.zeros((n, n_components))        # X scores
U = np.zeros((n, n_components))        # Y scores
P = np.zeros((p, n_components))     # X loadings
Q = np.zeros((q, n_components))     # Y loadings
B = np.zeros((n_components))                   # Regression coefficients

iter_comp = np.zeros((n_components))     # Stores iterations for each component


#==================== Loop over components ====================#
# Loop over defined number of components
for h in range(0, n_components):


    #==================== Tune component (PLS1/PLS2) ====================#
    # NOTE: Any commented sections with '(OPTIONAL)' beside them are optional rescailings of quantities which are
    # consistent with the PLS tutorial: Partial Least Squares Regression: A tutorial - Paul Geladi & 
    # Bruce Kowalski (1986)
    
    # PLS1 if Y is univariate, PLS2 if Y is multivariate
    if type == 'PLS1':
      
      # Initiate u as the only column of Y
      u = Y_h
      
      # Update X weights and normalise
      w = X_h.T.dot(u) / np.sum(u ** 2)
      w = w / (np.sum(w ** 2) ** 0.5)
      
      # Calculate X scores by regressing X on X weights
      t =  X_h.dot(w)
      # t =  t / np.sum(w ** 2)     # (OPTIONAL)
      
      # Y has one variable and therefore it is assigned a weight of 1. Similar for loadings
      c = 1
      q = 1
      
    else:
      
        # Initiate u as 1st column of Y
        u = Y_h[:,0]
        
        # Initiate iterables
        counter = 0
        t_diff = X_h.shape[1]
        old_t = 0

        # Loop until convergence or max iteration (for multivariate Y)
        while ((len(t_diff[np.absolute(t_diff) > tol]) != 0 or len(u_diff[np.absolute(u_diff) > tol]) != 0) and counter < max_iter + 1):
            
            # Update X weights and normalise
            w = X_h.T.dot(u) / np.sum(u ** 2)
            w = w / (np.sum(w ** 2) ** 0.5)
          
            # Calculate X scores by regressing X on X weights
            t =  X_h.dot(w)
            # t =  t / np.sum(w ** 2)     # (OPTIONAL)
            
            # Update Y weights and normalise
            c = Y_h.T.dot(t) / np.sum(t ** 2)
            c = c / (np.sum(c ** 2) ** 0.5)
            
            # Calculate Y scores by regressing Y on Y weights
            u =  Y_h.dot(c)
            # u =  u / np.sum(c ** 2)     # (OPTIONAL)
      
      
            # Update iterables
            t_diff = t - old_t
            old_t = t
            counter = counter + 1


    # Check convergence
    if counter == (max_iter + 1):
        print('Warning! Max iteration reached. No convergence for component ' + str(h))
 
    # Add number of iterations to vector
    iter_comp[:, h] = counter
    
    
    #==================== Deflation step ====================#
    # Calculate X loadings
    p_old = X_h.T.dot(t) / np.sum(t ** 2)
    p = p_old
    # p = p_old / (np.sum(p_old ** 2) ** 0.5)     # (OPTIONAL)
    # t = t * (np.sum(p_old ** 2) ** 0.5)     # (OPTIONAL)
    # w = w * (np.sum(p_old ** 2) ** 0.5)    # (OPTIONAL)
    
    # Calculate Y loadings
    q_old = Y_h.T.dot(u) / np.sum(u ** 2)
    q = q_old
    # q = q_old / (np.sum(q_old ** 2) ** 0.5)     # (OPTIONAL)
    # u = u * (np.sum(q_old ** 2) ** 0.5)     # (OPTIONAL)
    # c = c * (np.sum(q_old ** 2) ** 0.5)     # (OPTIONAL)
    
    # Calculate regression coefficient for regression mode
    b = u.dot(t) / np.sum(t ** 2)
    
    # Deflate matrices
    X_h = X_h - np.outer(t, p)
    Y_h = Y_h - b * np.outer(t, c)
    
    
    #==================== Store results for each component ====================#
    # Store results
    W[:, h] = w
    C[:, h] = c
    T[:, h] = t
    U[:, h] = u
    P[:, h] = p
    Q[:, h] = q
    B[h] = b
    
    # Add row/column names
    W = pd.DataFrame(W, index = x_columns, columns = comp_no)
    C = pd.DataFrame(C, index = y_columns, columns = comp_no)
    T = pd.DataFrame(T, index = id_rows, columns = comp_no)
    U = pd.DataFrame(U, index = id_rows, columns = comp_no)
    P = pd.DataFrame(P, index = x_columns, columns = comp_no)
    Q = pd.DataFrame(Q, index = y_columns, columns = comp_no)
    B = pd.DataFrame(B.reshape(1, no_of_components), columns = comp_no)
  
  
    #==================== Form predictions using results ====================#
    # Create function for prediction
  
    # ??????
  
    #========================================================================#
  
  
    # Return final outputs
    n_components
    tol
    max_iter
    data = list(X = X, Y = Y)
    scores = list(X.scores = T, Y.scores = U)
    weights = list(X.weights = W, Y.weights = C)
    loadings = list(X.loadings = P, Y.loadings = Q)
    iterations = iter_comp
    regression_coef = B
    algorithm = type
    names = list(sample = row_names, X.columns = x.col_names, Y.columns = y.col_names
    
    
    
    return(???)
  