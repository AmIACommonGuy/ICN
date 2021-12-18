#' NICE_fast_cut
#'
#' @param corr_m correlation matrix
#' @param threshold
#' @param cutter A sequence jump during iteration of K, default 1.
#' @return Cindx: the cluster index of every non-isolated node
#' @return CID: the cluster index of every cluster in a power 
#' descending order. i.e. CID(1) will be the cluster index of the 
#' most concentrated cluster
#' @return Clist: the reordered node index, nodes in the same 
#' cluster are permuted together.
#' @return K_selected: number of communities identified
#' @return t.overall: Time counter for NICE 
#' @return Prp_net: evaluation index of K
#' 
#' @export
#'
#' @examples
#' check sample in vignette in detail
NICE_fast_cut = function(corr_m, thres, cutter = 1) {
  # time count overall 
  t.overall1 = as.numeric(proc.time()[1])
  
  if (corr_m[1,1] == 1) {
    W = corr_m - diag(nrow(corr_m))
  } else {
    W = corr_m 
  }
  W[which(W<thres)] = 0 # In simu1 and thres=0.2, sparsity= 89.8%
  
  # (add later) use sparse property if sparsity is large. Don't use if not.
  z1 = which(sum(W)>0)
  W = W[z1, z1]
  n_z1 = length(z1)
  
  ##  Build laplacian matrix
  degs = sum(W)
  D = diag(degs)
  Lsp = D-W
  Leig = eigen(Lsp, n_z1)
  
  ## Find the best K (number of groups)
  lenW = length(which(W>0))/2 # total number of edges
  row = nrow(Lsp)
  # Values of K used for iteration
  K_vec = seq(from = 2, to = (nrow(Lsp) - 1), by = cutter)
  #### time counting vectors
  t.eigen <- t.kmeans <- t.obj <- rep(0, length=nrow(Lsp) - 2)
  
  vectors = Leig$vectors
  
  ncores = 4
  cl <- makeCluster(ncores)
  clusterEvalQ(cl, {library('rARPACK')
    library('matlab')})
  Prp_net = parSapply(cl, K_vec, Prpnet, Leig = vectors, W = W,  lenW = lenW)
  stopCluster(cl)

  best_idx = which(Prp_net == max(Prp_net))
  K = K_vec[best_idx]
  K = K[1] # Best K. In case several k's give the same Prp_net value
  K_selected = K
  
  ## Run keams again using the best K
  Leig = eigs(Lsp, K, which = "SM") 
  U = Leig$vectors
  C = kmeans(U,K, iter.max = 40) # C = MyKmeans(U,K)
  indx <- A_net <- net_V <- C_net <- rep(0, length = K)
  for (k in 1:K) {
    indx_k = which(C$cluster==k) # indices of nodes in cluster k
    indx[k:(k+length(indx_k)-1)] = indx_k # reordered indices
    net_V[k] = length(indx_k) # number of nodes in cluster k
    WC = W[indx_k, indx_k] # corr matrix of cluster k
    C_net[k] = sum(WC[lower.tri(WC)]) # sum of weights in cluster k
  }
  A_net = net_V*(net_V-1)/2 # number of potential connections
  
  ## Reorder the groups (purpose of display) in terms of group power
  diagscore = (C_net)^2/(A_net)/lenW
  diagscore[is.na(diagscore)] = 0
  diagscore_sortID=order(diagscore,decreasing = TRUE)
  
  ## Reorder the original indices
  inx_imporance = vector()
  for (i in 1:K) {
    inx_imporance = c(inx_imporance, which(C$cluster==diagscore_sortID[i]))
  }
  
  Cindx = seq(1,nrow(corr_m))
  Cindx[z1] = C$cluster
  Cindx[-z1] = 0
  CID = diagscore_sortID
  Clist = z1[inx_imporance]
  Clist = c(Clist, Cindx[-z1])
  
  #### end time counting
  t.overall2 = as.numeric(proc.time()[1])
  t.overall = t.overall2 - t.overall1
  
  return(list(Cindx = Cindx, CID = CID, Clist = Clist, K_selected = K_selected,  
              t.overall = t.overall, 
              Prp_net = Prp_net))
}

