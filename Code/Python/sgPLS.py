#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 21:23:14 2019

@author: raphaelsinclair
"""



# Import libraries
import numpy as np
import pandas as pd
from scipy.optimize import brentq



# Load dataframes (NOTE: linnerud dataset is incorrectly labelled. Error is preserved for consistency)
path = '/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH030 - Translational Data Sciences/Project'

X_data = pd.read_csv(path+'/Data/Xdata_simulated.csv', sep = ',')
Y_data = pd.read_csv(path+'/Data/Ydata_simulated.csv', sep = ',')


# Scale data ignoring missing values (Use unbiased estimator of standard deviation i.e. divisor = (n-1))

X_data = (X_data - np.nanmean(X_data, axis = 0))/np.nanstd(X_data, axis = 0, ddof = 1)
Y_data = (Y_data - np.nanmean(Y_data, axis = 0))/np.nanstd(Y_data, axis = 0, ddof = 1)



#=================================================#

X = X_data
Y = Y_data


n_components = 2

# Define ind.block.x/ind.block.y (i.e. vector of indices denoting the end of each group inclusive
# e.g. ind.block.x = c(6, 16) <==> 3 groups s.t. group 1 = 1-6, group 2 = 7-16, group 3 = 17-ncol(X))
ind_block_x = np.array(range(20, 400, 20))
ind_block_y = np.array(range(20, 500, 20))

# Select keepX_group/keepY_group variables (i.e. number of groups to keep in each component)
# keepY_group = rep(l, n_components)
keepX_group = np.array([4, 4])
keepY_group = np.array([4, 4])


alpha_x = np.array([0.95, 0.95])
alpha_y = np.array([0.95, 0.95])

max_lambda = 1e+05
#=================================================#




#--------------------------------------- sgPLS ---------------------------------------#

# def sgPLSreg(X, Y, n_components, keepX_group, keepY_group = None, ind_block_x, ind_block_y = None, alpha.x, alpha.y = None, tol = 1e-06, max_iter = 100, max_lambda = 1e+05):


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

if len(keepX_group) != n_components:
    return('Length of keepX does not match number of components')
  
if keepY_group.any != None and len(keepY_group) != n_components:
    return('Length of keepY does not match number of components')

if (len(ind_block_x) + 1) >= X_h.shape[1]:
    return('Length of ind_block_x exceeds the number of variables of X')

# Condition for length(indblocky) < q and indblocky = None if q=1


if len(keepX_group[keepX_group > (len(ind_block_x) + 1)]) != 0:
    return('keepX exceeds the number of groups in X')
  
if len(keepY_group[keepY_group > (len(ind_block_y) + 1)]) != 0:
    return('keepY exceeds the number of groups in Y')


# length and values of alpha 





#==================== Initiate items ====================#
# Dimensions of data
n = X_h.shape[0]
p = X_h.shape[1]
if len(Y_h.shape) == 1:
    q = 1
elif len(Y_h.shape) == 2:
    q = Y_h.shape[1]
k = len(ind_block_x) + 1
l = len(ind_block_y) + 1

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


# Number of variables to penalise (based on keepX_group/keepY_group)
x_sparsity = k - keepX_group

if keepY_group.any == None:
    keepY_group = k * np.ones((n_components))

y_sparsity = l - keepY_group


#==================== Create blocks of X and Y data ====================#
# Add {0, total variables} onto ind_block vectors (can more intuitively get group information from this vector)
ind_x = np.insert(ind_block_x, obj = [0, k - 1], values = [0, p])
ind_y = np.insert(ind_block_y, obj = [0, l - 1], values = [0, q])
  
# p_k/q_l holds the number of variables in each group for X and Y respectively
# x_blocks/y_blocks holds indices for the groups of variables in X and Y respectively
x_blocks = []
p_k = []
for index in range(0, k):
    p_k.append(ind_x[index + 1] - ind_x[index])
    x_blocks.append(range(ind_x[index], ind_x[index + 1])) 

y_blocks = []
q_l = []
for index in range(0, l):
    q_l.append(ind_y[index + 1] - ind_y[index])
    y_blocks.append(range(ind_y[index], ind_y[index + 1])) 
    

#==================== Define key functions ====================#
# Soft thresholding function
def soft_threshold(vec, lamb):
    return np.multiply(np.sign(vec), np.clip(np.absolute(vec) - lamb, a_min = 0, a_max = None))

# lambda quadratic function (want to solve for each group)
def lambda_quadratic(lamb, vec, alpha):
    g_soft = soft_threshold(vec, (lamb * alpha)/2)
    return np.sum(g_soft ** 2) - ((lamb * (1 - alpha)) ** 2) * len(vec)


#==================== Loop over components ====================#
# Loop over defined number of components
for h in range(0, n_components):
    
    #==================== Tune component ====================#
    # To penalise the number of groups according to sparse group penalties, the lambda threshold for group penalisation must
    # be found for each group. After determining the thresholds and the number of groups to be penalised, the sparse group
    # sparsity parameter is chosen (for the given alpha and block structure) s.t. the groups with the smallest contribution
    # to the covariance matrix are penalised.
    
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
      
        # Calculate group lasso penalties using the entries correspeonding to each group from the projection vector
        x_penalties = []
        for group in range(0, k):
            vec = M_v[x_blocks[group]]
            x_lam = brentq(lambda_quadratic, a = 0, b = max_lambda, args = (vec, alpha_x[h]), rtol = (np.finfo(float).eps**0.25), maxiter = 1000)
            x_penalties.append(x_lam)
            
        # Convert number of penalised groups in X into X sparsity parameter for each component
        if x_sparsity[h] == 0:
            lambda_x = sorted(x_penalties)[index] - 1
        else:
            index = x_sparsity[h] - 1     # Adjust for zero-indexing
            lambda_x = sorted(x_penalties)[index]

        # Optimise u iteratively for each group (and normalise)
        tmp = []
        for group in range(0, k):
            if x_penalties[group] <= lambda_x:
                u_k = np.zeros(p_k[group])
            else:
                vec = M_v[x_blocks[group]] 
                g = soft_threshold(vec, lambda_x * alpha_x[h]/2)
                u_k = 0.5 * (g - ((1 - alpha_x[h]) * (p_k[group] ** 0.5 / (np.sum(g ** 2) ** 0.5))) * g)
            
            tmp = np.concatenate((tmp, u_k), axis = None)
        u_new = tmp / (np.sum(tmp ** 2) ** 0.5)
            

        # Calculate the projection of u on M to produce the Y loadings candidate
        M_u = M.T.dot(u_new)
        
        # Calculate group lasso penalties using the entries correspeonding to each group from the projection vector
        y_penalties = []
        for group in range(0, l):
            vec = M_u[y_blocks[group]]
            y_lam = brentq(lambda_quadratic, a = 0, b = max_lambda, args = (vec, alpha_y[h]), rtol = (np.finfo(float).eps**0.25), maxiter = 1000)
            y_penalties.append(y_lam)
             
        # Convert number of penalised groups in Y into Y sparsity parameter for each component
        if y_sparsity[h] == 0:
            lambda_y = sorted(y_penalties)[index] - 1
        else:
            index = y_sparsity[h] - 1     # Adjust for zero-indexing
            lambda_y = sorted(y_penalties)[index]

        # Optimise v iteratively for each group (and normalise)
        tmp = []
        for group in range(0, l):
            if y_penalties[group] <= lambda_y:
                v_l = np.zeros(q_l[group])
            else:
                vec = M_u[y_blocks[group]] 
                g = soft_threshold(vec, lambda_y * alpha_y[h]/2)
                v_l = 0.5 * (g - ((1 - alpha_y[h]) * (q_l[group] ** 0.5 / (np.sum(g ** 2) ** 0.5))) * g)
            
            tmp = np.concatenate((tmp, v_l), axis = None)
        v_new = tmp / (np.sum(tmp ** 2) ** 0.5)
        
        
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
    keepX_group
    keepY_group
    ind_block_x
    ind_block_y
    alpha_x
    alpha_y
    tol
    max_iter
    max_lambda
    data = list(X = X, Y = Y)
    scores = list(X.scores = Eta, Y.scores = Omega)
    loadings = list(X.loadings = U, Y.loadings = V)
    defl.coefs = list(C = C, D = D)
    iterations = iter_comp
    names = list(sample = row_names, X.columns = x_columns, Y.columns = y_columns
    
    
    
    return(???)
  
  
