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


