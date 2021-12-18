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

# Prpnet score calculated function which is paralleled.
Prpnet = function(K, Leig, W, lenW){
  n_z1 = dim(W)[1]
  U = Leig[,n_z1:(n_z1-K+1)]
  C = kmeans(U,K, iter.max = 40) # C = MyKmeans(U,K)
  indx <- A_net <- net_V <- C_net <- rep(0, length = length(K))
  for (k in 1:K) {
    indx_k = which(C$cluster==k) # indices of nodes in cluster k
    indx[k:(k+length(indx_k)-1)] = indx_k # reordered indices
    net_V[k] = length(indx_k) # number of nodes in cluster k
    WC = W[indx_k, indx_k] # corr matrix of cluster k
    ##### Improvement
    C_net[k] = sum(WC[lower.tri(WC)])
    # C_net[k] = sum(sum(WC))/2 # sum of weights in cluster k
    # A_net[k] = net_V[k]*(net_V[k]-1)/2 # number of potential connections in cluster k
  }
  A_net = net_V*(net_V-1)/2
  Prpnet = (sum(C_net)/lenW)*(sum(C_net)/sum(A_net))
  return(Prpnet)
}
