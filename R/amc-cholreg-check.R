Y = matrix(rnorm(100*4), ncol=4)

A = matrix(rnorm(16), ncol=4)
Y = Y %*% A + rnorm(4)
Ybar = apply(Y, 2, mean)

Sigma = cov(Y)

# this gives recursion but results are stored in reverse order :-(
S = Sigma[4:1, 4:1]


U = chol(solve(S))
U.inv = solve(U)

t(U) %*% U  # is Sigma^-1
solve(S)
(U.inv) %*% t(U.inv) # Sigma
S

betas = diag(rep(1,4))-diag((1/diag(U))) %*% U

alpha = Ybar[4:1] - betas %*% Ybar[4:1]
# Conditional mean 
# E[Y_j | Y <j] = E[Y_j] + S_{j<j} S{<j<j}^{-1}(Y_{<j} - E[Y_{<j}])
#               = E[Y_j] + Beta (Y_{<j} - E[Y_{<j}])
#     alpha     = E[Y_j] - Beta E[Y_{<j}] + Beta Y_{<j}
cbind(alpha, betas)
# verified that is correct for slopes and intercepts
coef(lm(Y[,4] ~ Y[,3:1]))
coef(lm(Y[,3] ~ Y[,2:1]))
coef(lm(Y[,2] ~ Y[,1]))


# Sigma = L^T L
# Sigma^-1 = (L-1 L^-T) = L U   


S = Sigma  # do not permute!  
U.inv = solve(chol(S))
betas = diag(rep(1,4)) - diag(1/diag(U.inv)) %*% t(U.inv)
alpha = Ybar - betas %*% Ybar
# E[Y_j] = Ybar_j + betas(Y_<j - Ybar_<j)

cbind(alpha, betas)
# verified that is correct for slopes and intercepts

coef(lm(Y[,2] ~ Y[,1]))
coef(lm(Y[,3] ~ Y[,1:2]))
coef(lm(Y[,4] ~ Y[,1:3]))


#note this gives the right solution too and in the correct order!!!

# C/FORTRAN Check
#SS: 200 iterations
#0  0.955000  191.000000 159.000000 14.000000 46.000000 
#1  0.805000  0.000000 161.000000 11.000000 16.000000 
#2  0.080000  0.000000 0.000000 16.000000 6.000000 
#3  0.265000  0.000000 0.000000 0.000000 53.000000 

if (interactive()) {
mean.cov = scan()
0  0.950000  190.000000 137.000000 55.000000 73.000000 
1  0.715000  0.000000 143.000000 39.000000 25.000000 
2  0.315000  0.000000 0.000000 63.000000 34.000000 
3  0.410000  0.000000 0.000000 0.000000 82.000000 

mean.cov = matrix(mean.cov, ncol=6, byrow = T)
ybar = mean.cov[,2]
SS = mean.cov[1:4, 3:6]
m = 200; lambda = 1
Cov = (SS + t(SS))
diag(Cov) = diag(SS)
Cov = (Cov - m*ybar %*% t(ybar))/m + diag(lambda,4)
U = chol(Cov)
U.inv = solve(U)
diag(1, 4) - diag(1/diag(U.inv))%*%t(U.inv)
}