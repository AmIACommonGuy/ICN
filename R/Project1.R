NICE = function(corr_m, thres) {
  if (corr_m[1,1] == 1) {
    W = corr_m - diag(nrow(corr_m))
  } else {
    W = corr_m 
  }
  W[which(W<thres)] = 0 # In simu1 and thres=0.2, sparsity= 89.8%
  
  # (add later) use sparse property if sparsity is large. Don't use if not.
  
  z1 = which(sum(W)>0)
  W = W[z1, z1]
  
  ##  Build laplacian matrix
  degs = sum(W)
  D = diag(degs)
  Lsp = D-W
  
  ## Find the best K (number of groups)
  Prp_net = rep(0, nrow(Lsp)-2)
  iter = rep(0, nrow(Lsp)-2)
  lenW = length(which(W>0))/2 # total number of edges
  time = rep(0, nrow(Lsp)-2)
  for (K in 2:(nrow(Lsp)-1)) {
    time0 = as.numeric(proc.time()[1])
    L_eig = eigs(Lsp, K, which = "SM", use.arpack = TRUE) # K smallest eig value. Descent order.
    U = L_eig$vectors
    C = kmeans(U,K) # C = MyKmeans(U,K)
    indx = rep(0, length = K)
    A_net = rep(0, length = K)
    net_V = rep(0, length = K)
    C_net = rep(0, length = K)
    for (k in 1:K) {
      indx_k = which(C$cluster==k) # indices of nodes in cluster k
      indx[k:(k+length(indx_k)-1)] = indx_k # reordered indices
      net_V[k] = length(indx_k) # number of nodes in cluster k
      WC = W[indx_k, indx_k] # corr matrix of cluster k
      C_net[k] = sum(sum(WC[which(WC>0)]))/2
      A_net[k] = net_V[k]*(net_V[k]-1)/2 # number of potential connections in cluster k
    }
    time[K] = as.numeric(proc.time()[1]) - time0
    Prp_net[K]=(sum(C_net)/lenW)*(sum(C_net)/sum(A_net))
    iter[K] = C$iter
  }
  K = which(Prp_net == max(Prp_net))
  K = K[1] # Best K. In case several k's give the same Prp_net value
  K_selected = K
  
  ## Run keams again using the best K
  L_eig = eigs(Lsp, K, which = "SM") 
  U = L_eig$vectors
  C = kmeans(U,K) # C = MyKmeans(U,K)
  indx = rep(0, length = K)
  net_V = rep(0, length = K)
  C_net = rep(0, length = K)
  A_net = rep(0, length = K)
  for (k in 1:K) {
    indx_k = which(C$cluster==k) # indices of nodes in cluster k
    indx[k:(k+length(indx_k)-1)] = indx_k # reordered indices
    net_V[k] = length(indx_k) # number of nodes in cluster k
    WC = W[indx_k, indx_k] # corr matrix of cluster k
    C_net[k] = sum(sum(WC[which(WC>0)]))/2 # sum of weights in cluster k
    A_net[k] = net_V[k]*(net_V[k]-1)/2 # number of potential connections
  }
  
  ## Reorder the groups (purpose of display) in terms of group power
  diagscore = (C_net)^2/(A_net)/lenW
  diagscore[is.na(diagscore)]=0
  diagscore_sortID=order(diagscore,decreasing = TRUE)
  
  ## Reorder the original indices
  inx_imporance = vector()
  for (i in 1:K) {
    inx_imporance = c(inx_imporance, which(C$cluster==diagscore_sortID[i]))
  }
  
  Cindx = seq(1,nrow(corr_m))
  Cindx[z1] = C$cluster
  Cindx[-z1] = -1
  CID = diagscore_sortID
  Clist = z1[inx_imporance]
  Clist = c(Clist, Cindx[-z1])
  return(list(Cindx = Cindx, CID = CID, Clist = Clist, K_selected = K_selected, time = time, iter = iter, Prp_net = Prp_net))
}


### Data Simulation  ####

## Simulate Two Communities (100 nodes)

# Simulate true corr matrix
true1 = matrix(0, nrow = 100, ncol = 100)
true1[1:20,1:20] = 0.5
true1[21:30,21:30] = 0.5
true1[1:30,1:30] = true1[1:30,1:30] + 0.5*diag(30)
true1[31:100,31:100] = true1[31:100,31:100] + diag(70);
heatmap.2(true1, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", breaks=seq(0,1,0.01),
          density.info = "none")

# Simulate samples based on true corr matrix, and then calculate the sample corr matrix.
nvars = 100
nobs = 50
R = chol(true1) # upper triangular
Z_vec = rnorm(nvars * nobs)
Z = matrix(Z_vec, nrow = 100)
samples =  crossprod(R,Z) # can use crossprod
samples = t(samples) # For each node (100 nodes in total), draw 50 samples.
simu1 = cor(samples) # Get corr between nodes
heatmap.2(simu1, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

# Shuffle the sample corr matrix
shuffle_index = sample(1:100)
simu1_shuf = simu1[shuffle_index, shuffle_index]
heatmap.2(simu1_shuf, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

result = NICE(simu1_shuf, 0.2)
system.time(NICE(simu1_shuf, 0.2))

# Run simulation again when encounter error
# When result$Clist has -1, code below will have error. -1 corresponds to zero row in W (Line 79).
simu1_reordered = simu1_shuf[result$Clist, result$Clist]

# Sample corr matrix under true order
heatmap.2(simu1, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

# Sample corr matrix odered by NICE, after shuffling the true order
heatmap.2(simu1_reordered, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

