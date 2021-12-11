# Auxiliary Functions

# Find all pairs (n is number of indices)
lst_pairs = function(n){
  full.pairs <- matrix(c(rep(seq(1:n), each = n), rep(seq(1:n), n)), ncol = 2)
  pairs <- full.pairs[full.pairs[,1] > full.pairs[,2],]
  return(pairs[,2:1])
}

