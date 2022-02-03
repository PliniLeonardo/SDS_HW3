# Part 1 ------------------------------------------------------------------

library(clrdag)
set.seed(2018)
p <- 50; n <- 1000; sparsity <- 2/p

## generate a random lower triangular adjacnecy matrix
A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1))*runif(p*p, 0.7, 1), p, p)
A[upper.tri(A, diag=TRUE)] <- 0

## permute the order of adjacency matrix
idx <- sample(1:p)
A <- A[idx,idx]

## num of edges in A
sum(A != 0) # 43

## data matrix
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))

## estimate the graph
t <- proc.time()
out <- MLEdag(X=X, tau=0.3, mu=1, rho=1.2)
proc.time() - t 

## Frobenius distance to the truth adjacency matrix
sum((out$A - A)^2) # 0.0285789

## Hamming distance to the truth graph
sum(abs((out$A != 0) - (A != 0))) # 0

## test edge 1 --> 3
D <- matrix(0, p, p)
D[3,1] <- 1
out <- MLEdag(X=X, D=D, tau=0.3, mu=1, rho=1.2)
out$pval # 0.7496623

## test edge 7 --> 4
D <- matrix(0, p, p)
D[4,7] <- 1
out <- MLEdag(X=X, D=D, tau=0.3, mu=1, rho=1.2)
out$pval # 8.827349e-155

#randomly split
'''
size=floor(0.5*nrow(X))
indexes=sample(seq_len(nrow(X)),size=size)
D_train=X[indexes,]
D_test=X[-indexes,]
out=MLEdag(D_train,tau=0.3,mu=1,rho=1.2)
'''
# Part 2 ------------------------------------------------------------------

set.seed(2018)
p <- 50; n <- 1000;sparsity <- 2/p
## generate a random lower triangular adjacnecy matrix
A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1))*runif(p*p, 0.7, 1), p, p)
A[upper.tri(A, diag=TRUE)] <- 0

## num of edges in A
sum(A != 0) # 43
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

#generate D matrix 
#number of hyphothesized edges
q=43
D <- matrix(0, p, p)
index1=sample(1:p,q,replace=TRUE)
index2=sample(1:p,q,replace=TRUE)
for (i in 1:q){
  D[index1[i],index2[i]]=1
}

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
# U=35.95, so for alpha=0.05, we reject the null hypothesis

# Part 3 ------------------------------------------------------------------

#SIZE



