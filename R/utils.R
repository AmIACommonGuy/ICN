#' NICE post processing
#' 
#' `NICER` process NICE object for step 2 analysis
#' 
#' @param NICE A list object; result from NICE
#' @return Returns a list of: K = number of clusters; 
#' K.ind = indices of clusters
#' @example
#' 
#' @export

NICER = function(NICE){
  # obtain number of communities (clusters with > 1 nodes) `K`
  lst.indx = NICE$Cindx
  counts.indx = as.vector(table(lst.indx))
  K = sum(counts.indx > 1)
  K.ind = which(counts.indx > 1)
  
  return(list(K = K, K.ind = K.ind))
}
