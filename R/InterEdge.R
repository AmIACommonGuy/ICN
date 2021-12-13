#' Interedge
#'
#' @param mat1 community Vc
#' @param mat2 community Vc'
#' @param mat12 community V(c, c')
#' @lambda0 A hyper parameter in objective function
#' @rseq A sequence of cutoffs
#' @return rbest: best cutoff r  
#' 
#' 
#' @export
#'
#' @examples

#library(pracma)
InterCut = function(mat1, mat2, mat12, lambda0, rseq){
  obj = c()
  for(i in c(1: length(rseq))){
    vecW = mat12[abs(mat12) > rseq[i]]
    num = sum(abs(vecW)) + sum(abs(squareform(mat1))) + 
      sum(abs(squareform(mat2)))
    denom = length(vecW) + length(squareform(mat1)) + length(squareform(mat2))
    obj[i] = num / (denom^lambda0)
  }
  rbest = rseq(which.max(obj))
  return(rmbest)
}

