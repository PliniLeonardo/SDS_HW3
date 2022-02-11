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
'''
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


# Part 3 ------------------------------------------------------------------

# FUNZIONI ALE

compute.middle = function(X,A,j, p){
  # Finding the indexes K = {k s.t. k is in [1,p]\j}
  K = which((1:p)!=j)
  # in k we have the indexes of the row without j
  
  inners = apply(X[, K] * A[j, K], 1, sum) # n length vector containing the inner sums for each i-th sample
  middle = sum((X[, j] - inners)^2)  # as a matter of fact in the middle we are working on the j-th column
  return(middle)
}

sigma.estimator = function(X, A){
  n = dim(X)[1]
  p = dim(X)[2]
  return(sum(sapply(1:p, function(j)compute.middle(X,A,j, p)))/(n*p))
}

log_likelihood = function(X, A, sigma.square){
  n = dim(X)[1]
  p = dim(X)[2]
  return(-sum(sapply(1:p, function(j) (compute.middle(X,A,j,p)/(2*sigma.square)+(n/2*log(sigma.square))))))
}


inner_LRT_function.paths = function(X, X.tr, X.te, D){
  
  # Computing all the A.h0 over the sparsity parameters 'k'
  A.mles = lapply(1:sum(D), function(mu)MLEdag(X, D=D, tau=0.35, mu=mu, rho=1.2, trace_obj = F))
  A.h0s  = lapply(A.mles, function(a)a$A.H0)
  A.h1s  = lapply(A.mles, function(a)a$A.H1)
  
  # For each A.h0 compute the associated sigma
  sigmas.0 = t(sapply(A.h0s, function(a)sigma.estimator(X.tr, a)))
  sigmas.1 = t(sapply(A.h1s, function(a)sigma.estimator(X.te, a)))
  
  # Computing the likelihoods of each pair (A.h0, sigma.h0)
  likelihoods.0 = t(sapply(1:length(sigmas.0), function(idx)log_likelihood(X.tr, A.h0s[[idx]], sigmas.0[idx])))
  likelihoods.1 = t(sapply(1:length(sigmas.1), function(idx)log_likelihood(X.tr, A.h1s[[idx]], sigmas.1[idx])))
  
  # Finding the idx that give us the maximum likelihood
  mle_idx.0 = which(likelihoods.0==max(likelihoods.0))[1] # maximum likelihood over h0
  mle_idx.1 = which(likelihoods.1==max(likelihoods.1))[1] # maximum likelihood over h0
  
  return (likelihoods.1[mle_idx.1]-likelihoods.0[mle_idx.0])
  #return(list(U.n=exp(log_likelihood(X.tr, A, sigma.unconstrained)-likelihoods[mle_idx]), A.h0=A.h0s[[mle_idx]]))
  
}

inner_LRT_function.links = function(X, X.tr, X.te, D){
  tmp = MLEdag(X = X, D = D, tau = 0.35, mu = 1, 
               rho = 1.2, trace_obj = FALSE)
  A.h0 = tmp$A.H0
  A = tmp$A.H1
  
  sigma.h0 = sigma.estimator(X.tr, A.h0)
  sigma.unconstrained = sigma.estimator(X.te, A)
  
  return(log_likelihood(X.tr, A, sigma.unconstrained)-log_likelihood(X.tr, A.h0, sigma.h0))
}

log.LRT = function(X,D, links=T){
  n=dim(X)[1]
  X.tr = X[1:(n%/%2),]
  X.te = X[(n%/%2+1):n,]
  
  if(links){
    U_n = inner_LRT_function.links(X, X.tr, X.te, D)
    U_n.swap = inner_LRT_function.links(X, X.te, X.tr, D)
  }
  else{
    U_n = inner_LRT_function.paths(X, X.tr, X.te, D)
    U_n.swap = inner_LRT_function.paths(X, X.te, X.tr, D)
  }
  return(list(links = links, U_n = U_n, W_n = (log((exp(U_n)+exp(U_n.swap))/2))))
}


# SIZE --------------------------------------------------------------------

#proportion of rejections we get by applying our test on M dataset generated from
#a model compatible with H0=A[j,k]=0 for all (j,k) belong to F

### Parameters
set.seed(2018)
p <- 10
n <- 200
M <- 1000

### H0: F = { (p,2) }, and A[F] = 0
D <- matrix(0, p, p)
D[p,2] = 1

### Adjacency Matrix >> Hub
# All connected to 1, NO EDGE between p-2 >> **COMPATIBLE** with H0
A      <- matrix(0, p, p)     
A[, 1] <- sign( runif( p, min = -1, max = 1 ) )
A[1,1] <- 0

# Simulation --------------------------------------------------------------

M=100
cont=0
alpha=0.05
U_n <- rep(NA, M)
W_n <- rep(NA, M)
for(i in 1:M) {
  X   <- matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A) )
  temp=log.LRT(X,D, links=T)
  U_n[i]=temp$U_n
  W_n[i]=temp$W_n
  
  if (U_n[i]>log(1/alpha)) contU=contU+1
  if (W_n[i]>log(1/alpha)) contW=contW+1
}

cat('alpha is equal to for U ', contU/M,'alpha is equal to for W',contW/M)


# POWER -------------------------------------------------------------------

#Beta= the probability of failing to reject the null hypothesis when the null 
# hypothesis is false

### H0: F = { (p,1) }, and A[F] = 0
D <- matrix(0, p, p)
D[p,1] = 1

### Adjacency Matrix >> Hub
# All connected to 1, NO EDGE between p-2 >> **NOT COMPATIBLE** with H0
A      <- matrix(0, p, p)     
A[2:p, 1] <- sign( runif( p-1, min = -1, max = 1 ) )
A[3,4]=1
A[4,5]=1

'''
D.link = matrix(0, p, p)
D.link[2, 1] = 1
D.link[6, 3] = 1
D.link[10, 9] = 1
# Building A accordingly to h0
A.link.hidden = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A.link.hidden[upper.tri(A.link.hidden, diag = T)] = 0
A.link.hidden[2, 1] = 1#as.integer(!h0)
A.link.hidden[6, 3] = 1
A.link.hidden[10, 9] =1
'''
# SIMULATION
contU=0
contW=0
alpha=0.05
M=10
U_n <- rep(NA, M)
W_n <- rep(NA, M)
for(i in 1:M) {
  X   <- matrix( rnorm(n*p), n, p) %*% t(solve(diag(p) - A) )
  temp=log.LRT(X,D, links=T)
  U_n[i]=temp$U_n
  W_n[i]=temp$W_n
  
  if (U_n[i]>log(1/alpha)) contU=contU+1
  if (W_n[i]>log(1/alpha)) contW=contW+1
}

cat('1-beta is equal to for U: ', contU/M,'   ,1-beta is equal to for W:',contW/M)


# PART 4 ------------------------------------------------------------------
#import concatenated
data=read.csv("data/cell_signaling/concatenated.csv",header=T,sep=",")
head(data,n=5)
X=data.matrix(data, rownames.force = NA)
#data will be our X matrix
p=11

# Let's check the missed connections with the linkage type hypothesis
# So, we are working with a matrix NOT COMPATIBLE with H0:
# we already know that we have some problems in this case 

#FIRST CONNECTION: PIP2->PKC (columns 4  and column 9)
### H0: F = { (4,9) }, and A[F] = 0
D <- matrix(0, p, p)
D[9,4] = 1
alpha=0.05

temp=log.LRT(X,D, links=T)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK, we accept the null hypothesis

#FIRST CONNECTION: PIP2->PKC (columns 4  and column 9)
### H0: F = { (4,9) }, and A[F] = 0
D <- matrix(0, p, p)
D[4,9] = 1
alpha=0.05

temp=log.LRT(X,D, links=T)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK

#Let's continue
#SECOND CONNECTION: Plc-gamma -> PKC (column 3 and 9)
### H0: F = { (3,9) }, and A[F] = 0
D <- matrix(0, p, p)
D[3,9] = 1
alpha=0.05

temp=log.LRT(X,D, links=T)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK

#Let's continue
#THIRD CONNECTION: PIP3-> AKT (column 5 and 7)
### H0: F = { (4,7) }, and A[F] = 0
D <- matrix(0, p, p)
D[4,7] = 1
alpha=0.05

temp=log.LRT(X,D, links=T)
U_n=temp$U_n
W_n=temp$W_n

if(U_n>log(1/alpha)) cat(" We can reject the null hypothesis") else cat(" We can not reject the null hypothesis")
## OKKKK

#OSSERVAZIONI: Wn schizza ad infinito due volte
 

#PATHWAY LINKAGE
#check all the connections from PKC
#PKC->(RAF;MEK,ERK) which means columns (9->1,2,6)
D <- matrix(0, p, p)
D[9,1] = 1
D[1,2] = 1
D[2,6] = 1
alpha=0.05

temp=log.LRT(X,D, links=F)
U_n=temp$U_n
W_n=temp$W_n
cat("U_n:", U_n, "   W_n:", W_n)

##GOOD

# PART 5 ------------------------------------------------------------------

data=read.csv("data/cell_signaling/pma.csv",header=T,sep=",")
head(data,n=5)
X=data.matrix(data, rownames.force = NA)
hist(X[,9],breaks=100)

# non mi sembrano tutti normalmente distribuiti--> utilizzo scale o boxcox
# scale non mi pare porti lontano
# boxcox