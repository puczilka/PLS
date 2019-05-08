#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 00:31:15 2019

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



#--------------------------------------- sPLS ---------------------------------------#

# sPLSreg = def(X, Y, n_components, keepX, keepY = None, tol = 1e-06, max_iter = 100):


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

if len(keepX) != n_components:
    return('Length of keepX does not match number of components')
  
if keepY != None and len(keepY) != n_components:
    return('Length of keepY does not match number of components')
    
if len(keepX[keepX > X_h.shape[1]]) != 0:
    return('keepX exceeds the number of variables in X')
  
if len(keepY[keepY > Y_h.shape[1]]) != 0:
    return('keepY exceeds the number of variables in Y')


#==================== Initiate items ====================#


# Dimensions of data
n = X_h.shape[0]
p = X_h.shape[1]
if len(Y_h.shape) == 1:
    q = 1
elif len(Y_h.shape) == 2:
    q = Y_h.shape[1]

# Assign row/column labels
x_columns = X.columns
y_columns = Y.columns
id_rows = X.index

comp_no = []
for index in range(1, n_components + 1):
    comp_no.append('comp '+ str(index))
  

# Initiate matrices to store results
Eta = np.zeros((n, n_components))     # X scores
Omega = np.zeros((n, n_components))     # Y scores
U = np.zeros((p, n_components))     # X loadings
V = np.zeros((q, n_components))     # Y loadings
C = np.zeros((p, n_components))     # Regression coefficient on latent variables (for X_h)
D = np.zeros((q, n_components))     # Regression coefficient on latent variables in 'Regression mode' (for Y_h)

iter_comp = np.zeros((n_components))     # Stores iterations for each component


# Number of variables to penalise (based on keepX/keepY)
x_sparsity = X_h.shape[1] - keepX

if keepY == None:
    keepY = q * np.ones((n_components))

y_sparsity = Y_h.shape[1] - keepY


#==================== Loop over components ====================#
# Loop over defined number of components
for h in range(0, n_components):
    
    #==================== Tune component ====================#
    # WARNING: Truncation errors may occur from long format numbers in calculations. These can be carried over during iterative steps.
    
    # NOTE: Conversion of sparsity parameter to number of variables is made by calculating the magnitude 
    # of each variables contribution to the covariance structure between X and Y.
    # If we chose to penalise k parameters (in X w.l.o.g), we choose lambda to be the kth smallest magnitude.
    # Therefore lambda is chosen s.t. its value is sufficiently large enough to penalise the k entries which 
    # contribute least (in the loadings). Similarly for penalisation on Y
    
    # Compute matrix M (p x q matrix)
    M = X_h.T.dot(Y_h)
    
    # Find SVD of M and extract loadings (first pair of left and right singular vectors)
    U_tilda, Diag, Vh_tilda = np.linalg.svd(M, full_matrices = False)
    
    u_old = U_tilda[:, [0]]
    v_old = Vh_tilda.T[:, [0]]
    
    # Initiate iterables
    counter = 0
    u_diff = X_h.shape[1]
    v_diff = Y_h.shape[1]

    # Loop until convergence or max iteration (for multivariate Y)
    while ((len(u_diff[np.absolute(u_diff) > tol]) != 0 or len(v_diff[np.absolute(v_diff) > tol]) != 0) and counter < max_iter + 1):
        
        # Calculate the projection of v on M to produce the X loadings candidate
        M_v = M.dot(v_old)
      
        # Convert number of penalised variables in X into X sparsity parameter for each component
        if x_sparsity[h] == 0:
            lambda_x = 0
        else:
            index = x_sparsity[h]
            lambda_x = sorted(np.absolute(M_v))[index]

        # Optimise u iteratively using soft-tresholding penalisation function (and normalise)
        u_new = np.multiply(np.sign(M_v), np.clip(np.absolute(M_v) - lambda_x, a_min = 0, a_max = None))
        u_new = u_new / (np.sum(u_new ** 2) ** 0.5)
        
        
        # Calculate the projection of u on M to produce the Y loadings candidate
        M_u = M.T.dot(u_new)
        
        # Convert number of penalised variables in Y into Y sparsity parameter for each component
        if y_sparsity[h] == 0:
            lambda_y = 0
        else:
            index = y_sparsity[h]
            lambda_y = sorted(np.absolute(M_u))[index]

        # Optimise u iteratively using soft-tresholding penalisation function (and normalise)
        v_new = np.multiply(np.sign(M_u), np.clip(np.absolute(M_u) - lambda_y, a_min = 0, a_max = None))
        v_new = v_new / (np.sum(v_new ** 2) ** 0.5)
        
        
        # Update iterables
        u_diff = u_new - u_old
        v_diff = v_new - v_old
        u_old = u_new
        v_old = v_new
        counter = counter + 1
        
        
    # Check convergence
    if counter == (max_iter + 1):
        print('Warning! Max iteration reached. No convergence for component ' + str(h))
 
    # Add number of iterations to vector
    iter_comp[:, h] = counter
        
    
    #==================== Deflation step ====================#
    # Calculate scores/latent variables for X and Y
    eta = X_h.dot(u_new) / np.sum(u_new ** 2)
    omega = Y_h.dot(v_new) / np.sum(v_new ** 2)
    
    # Calculate regression coefficients
    c = X_h.T.dot(eta) / np.sum(eta ** 2)
    d = Y_h.T.dot(eta) / np.sum(eta ** 2)
    # e = Y_h.T.dot(omega) / np.sum(omega ** 2)     # Canonical
    
    # Deflate X and Y matrices using latent variables and regression coefficients
    X_h = X_h - np.outer(eta, c)
    Y_h = Y_h - np.outer(eta, d)
    
    
    #==================== Store results for each component ====================#
    # Store variables
    Eta[:, h] = eta
    Omega[:, h] = omega
    U[:, h] = u_new
    V[:, h] = v_new
    C[:, h] = c
    D[:, h] = d
    # E[:, h] = e
    
    
    # Add row/column names
    Eta = pd.DataFrame(Eta, index = id_rows, columns = comp_no)
    Omega = pd.DataFrame(Omega, index = id_rows, columns = comp_no)
    U = pd.DataFrame(U, index = x_columns, columns = comp_no)
    V = pd.DataFrame(V, index = y_columns, columns = comp_no)
    C = pd.DataFrame(C, index = x_columns, columns = comp_no)
    D = pd.DataFrame(D, index = y_columns, columns = comp_no)
  
  
    #==================== Form predictions using results ====================#
    # Create function for prediction
    
    # ??????
    
    #========================================================================#
  
    
    # Return final outputs
    n_components
    keepX
    keepY
    tol
    max_iter
    data = list(X = X, Y = Y)
    scores = list(X.scores = Eta, Y.scores = Omega)
    loadings = list(X.loadings = U, Y.loadings = V)
    defl.coefs = list(C = C, D = D)
    iterations = iter_comp
    names = list(sample = row_names, X.columns = x_columns, Y.columns = y_columns
    
    
    
    return(???)
  
  
        
        