# MLEs and hypothesis -----------------------------------------------------
p = 250
n = 100

# Suppose to have the following linkages: F = {(1,2), (3,6), (9,30), (200,249)}
D = matrix(0, p, p)
D[1, 2] = 1
D[3, 6] = 1
D[9, 30] = 1
D[200, 249] = 1

# So X must have these linkages.
# Since X came from an adiacency matrix A, let's generate it 
A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = T)] = 0
A[1, 2] = 1
A[3, 6] = 1
A[9, 30] = 1
A[200, 249] = 1
# This it the REAL (UNOBSERVED) adiacency matrix A



# UNDER H0
A.0 = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A.0[upper.tri(A, diag = T)] = 0
A.0[1, 2] = 0
A.0[3, 6] = 0
A.0[9, 30] = 0
A.0[200, 249] = 0

# Random A without constraint (all theta, not only theta 0)
A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = T)] = 0

library(clrdag)
MLEdag()

# likelihood and sigma estimation -----------------------------------------
log_likeihood = function(X, A, sigma.square){
  #The formula is: tot = -(sum(middle/const1 + const2))
  # where middle = sum((i-th row of X - inner)^2)
  #  where inner= sum(j-th row of A \ j-th col * i-th row of X \ j-th col)
  tot = 0
  n = dim(X)[1]
  p = dim(X)[2]
  for(j in 1:p){
    
    # Finding the indexes K = {k s.t. k is in [1,p]\j}
    K = which((1:p)!=j)
    # in k we have the indexes of the row without j
      
    inners = apply(X[, K] * A[j, K], 1, sum) # n length vector containing the inner sums
    middle = sum((X[, j] - inners)^2)
    tot = tot - (middle/(2*sigma.square) + (n/2*log(sigma.square)))
  }
}

sigma.estimator = function(X, A){
  n = dim(X)[1]
  p = dim(X)[2]
  sigma = 0
  for(j in 1:p){
    K = which((1:p)!=j)
    inners = apply(X[, K] * A[j, K], 1, sum)
    sigma = sigma + sum((X[,j]-inners)^2)
  }
  return(sigma/(n*p))
}


generate.X = function(p, n, A=NULL, D=NULL, h0=T, linkages=T){
  # If A is not provided I'll generate it accordingly to eventual links and hypotesis sake
  if(is.null(A)){
    A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
    A[upper.tri(A, diag = T)] = 0
    if(!is.null(D))
      if(h0)
        A[D==1] = 1-as.integer(linkages) #A[j,k] = 0 for ALL linkages in F if we want to test the linkages, else 1
    else{
      A[sample(which(D==1), sample(sum(D==1), 1))] = as.integer(linkages) #A[j,k] = 0 for ALL the links in F if we are NOT testing the linkages
      cat('for invalidate h0 A has in positions', which(D==1), ' ', sum(A[which=D==1]), ' ones')
    }
  }
  
  return(matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A)))
}


sigma.estimator = function(X, A){
  n = dim(X)[1]
  p = dim(X)[2]
  sigma = 0
  for(j in 1:p){
    K = which((1:p)!=j)
    inners = apply(X[, K] * A[j, K], 1, sum)
    sigma = sigma + sum((X[,j]-inners)^2)
  }
  return(sigma/(n*p))
}


n = 100
p = 10
sparsity <- 1/p
m = 10
for ( i in 1:n){
  # Generating the hidden Adaicency matrix
  A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
  A[upper.tri(A, diag = T)] = 0
  # X is generated gaussian mantaining the dependencies described by the adiacency matrix A
  X = matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A) )
  sigma = sigma.estimator(X,A)
  cat('sigma: ', sigma, ' vectorized: ', vectorized.sigma.estimator(X, A),'\n')
}



