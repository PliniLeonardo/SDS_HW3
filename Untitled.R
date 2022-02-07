# Part 1 ------------------------------------------------------------------

library(clrdag)

# Part 2 ------------------------------------------------------------------

set.seed(2018)
p <- 25; n <- 2000;sparsity <- 2/p
## generate a random lower triangular adjacnecy matrix
A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1))*runif(p*p, 0.7, 1), p, p)
A[upper.tri(A, diag=TRUE)] <- 0

## num of edges in A
sum(A != 0) # 26
## data matrix
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))

#Loglikelihood
loglikelihood=function(A,sigma_square,X){
  n=ncol(X)
  p=nrow(X)
  tot=0
  for (j in 1:n){
    somma2=0
    for (i in 1:p){
      somma1=0
      for (k in 1:n){
        if (k!=j){
          somma1=somma1+A[j,k]*X[j,k]
        }
      }
      partial=X[i,j]-somma1
    }
    somma2=somma2+partial^2
  }
  tot=tot+(1/(2*sigma_square)*somma2+n*log(sigma_square)/2)
  return(-tot)
}
'''
log_likelihood = function(X, A, sigma.square){
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
'''
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

#generate D matrix 
#number of hyphothesized edges
'''
q=20
D <- matrix(0, p, p)
index1=sample(1:p,q,replace=TRUE)
index2=sample(1:p,q,replace=TRUE)
for (i in 1:q){
  D[index1[i],index2[i]]=1
}
'''

D = matrix(0, p, p)
D[2, 1] = 1
D[6, 3] = 1
D[10, 9] = 1
'''
#MLEdag method to recover A.H0 and A.H1
out=MLEdag(X=X,D=D,tau=0.3, mu=1, rho=1.2)
A_costrained=out$A.H0
A_uncostrained=out$A.H1

#loglikelihood as a function of just sigma squared
l_sigma=function(A,sigma_square,X){
  #perform the loglikelihood in function of sigma squared
  function1=function(sigma_square){
    loglikelihood(A=A,sigma_square=sigma_square,X=X)
  }
  return(function1)
}

#sigma_square,loglikelihood costrained

out_costrained=optim(par=c(1),fn=l_sigma(A_costrained,sigma_square,X),method=c("Brent"),lower = 0, upper = 10000,control = list(fnscale = -1))
best_sigma_squared_costrained=out_costrained$par
maximum_value_l_costrained=out_costrained$value
likelihood_costrained=exp(maximum_value_l_costrained)

#sigma_squared, loglikelihood uncostrained

out_uncostrained=optim(par=c(1),fn=l_sigma(A_uncostrained,sigma_square,X),method=c("Brent"),lower = 0, upper = 10000,control = list(fnscale = -1))
best_sigma_squared_uncostrained=out_uncostrained$par
maximum_value_l_uncostrained=out_uncostrained$value
likelihood_uncostrained=exp(maximum_value_l_uncostrained)


#Built the statistic U_n
U=likelihood_uncostrained/likelihood_costrained

#Comment
# U=1, so for alpha=0.05, we can not reject the null hypothesis
'''

# Part 3 ------------------------------------------------------------------

U_n.linkages = function(X, D){
  A.h0 = MLEdag(X = X.tr, D = D, tau = 0.35, mu = 1, 
                rho = 1.2, trace_obj = FALSE)$A.H0
  A = MLEdag(X = X, tau = 0.35, mu = 1, 
             rho = 1.2, trace_obj = FALSE)$A
  sigma.h0 = sigma.estimator(X.tr, A.h0)
  sigma.unconstrained = sigma.estimator(X.te, A)
  return(exp(loglikelihood( A, sigma.unconstrained,X.tr))/ exp(loglikelihood( A.h0, sigma.h0,X.tr)))
}


#SIZE

#proportion of rejections we get by applying our test on M dataset generated from
#a model compatible with H0=A[j,k]=0 for all (j,k) belong to F

### Parameters
set.seed(2018)
p <- 25
n <- 200
M <- 1000

### H0: F = { (p,2),(2,2) }, and A[F] = 0
D <- matrix(0, p, p)
D[p,2] = 1
D[2,2] = 1

### Adjacency Matrix >> Hub
# All connected to 1, NO EDGE between p-2 >> **COMPATIBLE** with H0
A      <- matrix(0, p, p)     
A[, 1] <- sign( runif( p, min = -1, max = 1 ) )
A[1,1] <- 0

# Simulation --------------------------------------------------------------

p_bar <- FALSE # progress-bar yes/no?
if (p_bar){
  # install.packages("svMisc")  
  suppressMessages( require(svMisc, quietly = T ) ); ?progress
}

str_time <- Sys.time()
cont=0
alpha=0.05
ris <- rep(NA, M)
for(i in 1:M) {
  X   <- matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A) )
  X.tr = X[1:(n%/%2),]
  X.te = X[(n%/%2+1):n,]
  ris[i]=U_n.linkages(X,D)
  
  if (ris[i]>(1/alpha)) cont=cont+1
  if (p_bar) progress(i, M)
}
stp_time <- Sys.time()
stp_time - str_time

cat('alpha is equal to ', cont/M)
# Take a look: (M = 1000, n = 200)
# Are we simulating under H0 as expected here? 
# If so, the pval should be Unif(0,1)

hist(pval, prob = T, border = "white")  # close enough!
abline(h = 1, lty = 2, col = "red")
mean(pval < 0.05)   # proportion under 0.05
mean(pval < 0.01)   # proportion under 0.01


# POWER -------------------------------------------------------------------
#Beta= the probability of failing to reject the null hypothesis when the null 
# hypothesis is false

### H0: F = {(,2) }, and A[F] = 0
D <- matrix(0, p, p)
D[,1] = 1
D[,10] = 1
D[,20] = 1

### Adjacency Matrix >> Hub
# All connected to 1, NO EDGE between p-2 >> ** NOT COMPATIBLE** with H0
A      <- matrix(0, p, p) 
A[, 1] = sign( runif( p, min = -1, max = 1 ) )
A[, 10] = sign( runif( p, min = -1, max = 1 ) )
A[, 20] = sign( runif( p, min = -1, max = 1 ) )

str_time <- Sys.time()
cont=0
M=10
ris <- rep(NA, M)
for(i in 1:M) {
  X   <- matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A) )
  #X.tr = X[1:(n%/%2),]
  #X.te = X[(n%/%2+1):n,]
  ris[i]=U_n.linkages(X,D)
  
  # we fail to reject if Un<(1/0.05)
  if (ris[i]<(1/0.05)) cont=cont+1
  if (p_bar) progress(i, M)
}
stp_time <- Sys.time()
stp_time - str_time

cat('1-Beta is equal to ', 1-cont/M)
