####### gPLS

#---------------------------- SIMULATED DATA ---------------------------------#

## Simulation of datasets X and Y with group variables
n <- 100
sigma.gamma <- 1
sigma.e <- 1.5
p <- 400
q <- 500
theta.x1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5,15), 
              rep(0, 5), rep(-1.5, 15), rep(0, 325))
theta.x2 <- c(rep(0, 320), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
              rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))

theta.y1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5, 15),
              rep(0, 5), rep(-1.5, 15), rep(0, 425))
theta.y2 <- c(rep(0, 420), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
              rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))                            


Sigmax <- matrix(0, nrow = p, ncol = p)
diag(Sigmax) <- sigma.e ^ 2
Sigmay <- matrix(0,nrow = q, ncol = q)
diag(Sigmay) <- sigma.e ^ 2

set.seed(125)

gam1 <- rnorm(n)
gam2 <- rnorm(n)

X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
                                                               nrow = 2, byrow = TRUE) + rmvnorm(n, mean = rep(0, p), sigma =
                                                                                                   Sigmax, method = "svd")
Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2), 
                                                               nrow = 2, byrow = TRUE) + rmvnorm(n, mean = rep(0, q), sigma =
                                                                                                   Sigmay, method = "svd")



#-------------------- MANUAL gPLS (penalising groups of vaiables) --------------------#

# gPLS ~= sPLS
# Loadings vectors are tuned with a different expression (not the soft thresholding seen in sPLS) and tuning
# is applied to all subgroups of X and/or subgroups of Y for their respective loadings vectors


# Convert data into matrices
X_h <- as.matrix(scale(X))
Y_h <- as.matrix(scale(Y))
# X_h <- as.matrix(X)
# Y_h <- as.matrix(Y)

# Select number of components
n_comp <- 1

# Define ind.block.x/ind.block.y (i.e. vector of indices denoting the end of each group inclusive
# e.g. ind.block.x = c(6, 16) <==> 3 groups s.t. group 1 = 1-6, group 2 = 7-16, group 3 = 17-ncol(X))
ind.block.x <- seq(20, 380, 20)
ind.block.y <- seq(20, 480, 20)

# Define number of groups
k = length(ind.block.x) + 1
l = length(ind.block.y) + 1

# Select keepX_group/keepY_group variables (i.e. number of groups to keep in each component)
# keepX_group = rep(k, n_comp)
# keepY_group = rep(l, n_comp)
keepX_group = c(4, 4)
keepY_group = c(4, 4)

# Number of groups to penalise
x_sparsity <- rep(k, n_comp) - keepX_group

if (is.null(ind.block.y)) {
  y_sparsity <- rep(0, n_comp)
} else {
  y_sparsity <- rep(l, n_comp) - keepY_group
}


# Create items to store results
Eta_mat <- NULL     # X scores
Omega_mat <- NULL     # Y scores
U_mat <- NULL     # X loadings
V_mat <- NULL     # Y loadings
C <- NULL     # Regression coefficient for X
D <- NULL     # Regression coefficient for Y in 'Regression mode'
E <- NULL     # Regression coefficient for Y in 'PLS-mode A'



# Calculate group indices from ind.block.x/ind.block.y (can more intuitively get group information from this vector
# instead of ind.block.x/ind.block.y)
ind.x <- c(0,ind.block.x,ncol(X_h))
ind.y <- c(0,ind.block.y,ncol(Y_h))


# Return vectors containing the number of variables in each group (the set of all p_k and q_l) and the indices for the
# blocks of variables for X and Y
x.blocks <- list()
p_k <- NULL
for (index in 1:k) {
  p_k[index] = ind.x[index + 1] - ind.x[index]
  x.blocks[[index]] = c((ind.x[index] + 1):ind.x[index + 1])
}

y.blocks <- list()
q_l <- NULL
for (index in 1:l) {
  q_l[index] = ind.y[index + 1] - ind.y[index]
  y.blocks[[index]] = ((ind.y[index] + 1):ind.y[index + 1])
}


###### LOOP OVER COMPONENTS ######
# Initiate first component and loop over defined number of components
# for (h in (1:n_comp)) {
h<-1

# Compute matrix M (p x q matrix)
M <- t(X_h) %*% Y_h

###### CALIBRATION OF COMPONENT ######
# Find SVD of M and find loadings from first pair of singular vectors
M_decomp <- svd(M, nu = 1, nv = 1)
u_old <- M_decomp$u
v_old <- M_decomp$v

# Initiate iterables
counter <- 0
max_iter <- 1
u_diff <- dim(X_h)[1]
v_diff <- dim(X_h)[1]
tol <- 1e-06
pen_tol <- .Machine$double.eps ^ 0.5

# Loop until convergence of u and v or max iteration
# while ((sum(abs(u_diff) > tol) != 0 || sum(abs(v_diff) > tol) != 0) && counter < max_iter) {

# Calculate the projection of v on M to produce the X loadings candidate
M_v <- M %*% v_old

# Calculate group lasso penalties using the entries correspeonding to each group from the projection vector
# and add to x.penalty vector
x.penalties <- NULL
for (group in 1:k) {
  vec <- M_v[(x.blocks[[group]])]
  x.penalties <- c(x.penalties, 2*normv(vec)/sqrt(p_k[group]))
}

# Convert number of penalised groups in X into sparsity parameter based on group lasso penalties
if (x_sparsity[h] == 0) {
  lambda_u <- 0
} else {
  lambda_u <- sort(x.penalties)[x_sparsity[h]]
}

# Optimise u iteratively for each group (and normalise)
# NOTE: 1st entry of pmax is used for non-unique maxima
tmp <- NULL
for (group in 1:k) {
  vec <- M_v[(x.blocks[[group]])]
  pen <- 1 - (lambda_u/x.penalties[group])
  if (pen < pen_tol) {pen <- 0}
  tmp <- c(tmp, pen*vec)
}
u_new = tmp / normv(tmp)



# Calculate the projection of u on M to produce the Y loadings candidate
M_u <- t(M) %*% u_new

# Calculate group lasso penalties using the entries correspeonding to each group from the projection vector
# and add to y.penalty vector
y.penalties <- NULL
for (group in 1:l) {
  vec <- M_u[(y.blocks[[group]])]
  y.penalties <- c(y.penalties, 2*normv(vec)/sqrt(q_l[group]))
}

# Convert number of penalised groups in Y into sparsity parameter based on group lasso penalties
if (y_sparsity[h] == 0) {
  lambda_v <- 0
} else {
  lambda_v <- sort(y.penalties)[y_sparsity[h]]
}

# Optimise v iteratively for each group (and normalise)
# NOTE: 1st entry of pmax is used for non-unique maxima
tmp <- NULL
for (group in 1:l) {
  vec <- M_u[(y.blocks[[group]])]
  pen <- 1 - (lambda_v/y.penalties[group])
  if (pen < pen_tol) {pen <- 0}
  tmp <- c(tmp, pen*vec)
}
v_new = tmp / normv(tmp)


# Update iterables
u_diff = u_new - u_old
v_diff = v_new - v_old
u_old = u_new
v_old = v_new
counter = counter + 1


# NOTE: Conversion of sparsity parameter to number of variables is made by calculating the part of the 
# the magnitude of the projection of the penalty on each entry of the loadings
# (i.e. abs(M %*% v_old)     for u;    abs(t(M) %*% u_new)     for v)
# and then choosing lambda to be the value with the smallest penalty to ensure that the variables are
# penalised according to their contribution to the loadings vectors.

# }

if (counter == max_iter + 1) {
  print(paste0('Warning! Max iteration reached. Component ', h))
}

###### DEFLATION STEP ######
# Calculate scores/latent variables for X and Y
eta = as.vector(X_h %*% u_new) / as.vector(t(u_new) %*% u_new)
omega = as.vector(Y_h %*% v_new) / as.vector(t(v_new) %*% v_new)


# Calculate regression coefficients
c = as.vector(t(X_h) %*% eta) / as.vector(t(eta) %*% eta)
d = as.vector(t(Y_h) %*% eta) / as.vector(t(eta) %*% eta)
e = as.vector(t(Y_h) %*% omega) / as.vector(t(omega) %*% omega)


# Deflate X and Y matrices using latent variables and regression coefficients
X_h <- X_h - (eta %*% t(c))
Y_h <- Y_h - (eta %*% t(d))

# Store variables
Eta_mat <- cbind(Eta_mat, eta, deparse.level = 0)
Omega_mat <- cbind(Omega_mat, omega, deparse.level = 0)
U_mat <- cbind(U_mat, u_new, deparse.level = 0)
V_mat <- cbind(V_mat, v_new, deparse.level = 0)
C <- cbind(C, c, deparse.level = 0)
D <- cbind(D, d, deparse.level = 0)
E <- cbind(E, e, deparse.level = 0)

# }

# Assign row/column names labels
colnames(U_mat) <- colnames(V_mat) <- colnames(Eta_mat) <- colnames(Omega_mat)  <- colnames(C) <- colnames(D) <- colnames(E) <- paste0('comp ', seq(n_comp))
rownames(C) <- colnames(X)
rownames(D) <- rownames(E) <- colnames(Y)


#=================================================#

#### gPLS model
model.gPLS <- gPLS(X, Y, ncomp = n_comp, mode = "regression", keepX = keepX_group, 
                   keepY = keepY_group, ind.block.x = ind.block.x , ind.block.y = ind.block.y, max.iter = max_iter, scale = TRUE)
