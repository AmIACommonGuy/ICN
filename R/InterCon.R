#' InterCon
#'
#' @param mat1 community Vc
#' @param mat2 community Vc'
#' @param mat12 community V(c, c')
#' @param lambda0 A hyper parameter in objective function
#' @param rseq A sequence of cutoffs
#' @param nodes1 nodes of the community c
#' @param nodes2 nodes of the community c'
#' @return result a list of 3 elements: rbest, interEdge and updated mat12
#' rbest is the best threshold.
#' InterEdge is a matrix including all significant edges
#' mat12 is the updated matrix after filtering by threshold.
#' @export
#' 
#' @import matlab
#'
InterCon = function(mat1, mat2, mat12, inter12, lambda0, rseq, nodes1, nodes2){
  # mat1 = mat1 - diag(mat1)
  # mat2 = mat2 - diag(mat2)
  # mat12 = mat12 - diag(mat12)
  diag(mat1) = 0
  diag(mat2) = 0
  diag(mat12) = 0
  obj = c()
  for(i in c(1: length(rseq))){
    vecW = mat12[abs(mat12) > rseq[i]]
    num = matlab::sum(abs(vecW)) + matlab::sum(abs(squareform(mat1))) + matlab::sum(abs(squareform(mat2)))
    denom = length(vecW) + length(squareform(mat1)) + length(squareform(mat2))
    obj[i] = num / (denom^lambda0)
  }
  rbest = rseq[which.max(obj)]
  result = list()
  result[[1]] = rbest
  interEdge = which(inter12 > rbest, arr.ind = TRUE)
  result[[2]] = interEdge
  inter12[inter12 < rbest] = 0
  n1 = nrow(mat1)
  n2 = nrow(mat2)
  mat12[c((n1+1):(n1+n2)), c(1:n1)] = inter12
  mat12[c(1:n1), c((n1+1):(n1+n2))] = t(inter12)
  result[[3]] = mat12
  return(result)
}

