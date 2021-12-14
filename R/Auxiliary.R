# Auxiliary Functions

# Find all pairs (n is number of indices)
lst_pairs = function(n){
  full.pairs <- matrix(c(rep(seq(1:n), each = n), rep(seq(1:n), n)), ncol = 2)
  pairs <- full.pairs[full.pairs[,1] > full.pairs[,2],]
  pairs <- to_matrix(pairs)
  
  return(pairs[,2:1])
}

# helper function to convert output into matrix:
to_matrix <- function(vec){
  mat <- as.matrix(vec)
  n <- dim(mat)[2]
  if (n == 1) {
    mat = t(mat)
  } 
  return(mat)
}
