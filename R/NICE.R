#' NICE
#'
#' @param corr_m correlation matrix
#' @param thres threshold
#'
#' @returns Cindx:    the cluster index of every non-isolated node
#' @returns CID:     the cluster index of every cluster in a power descending order. i.e. CID(1) will be the cluster index of the most concentrated cluster
#' @returns Clist:   the reordered node index, nodes in the same cluster are permuted together.
#' 
#' 
#' @export
#'
#' @examples
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
