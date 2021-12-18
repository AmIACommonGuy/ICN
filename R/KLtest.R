#' KL Divergence Calculation
#' 
#' `KL` Calculates the KL divergence between two distributions y1 and y2
#' @param y1 distribution 1 (vector)
#' @param y2 distribution 2 (vector)
#' @param width bin width. Default is 0.001
#' 
#' @return KL distance between two distributions
#' @examples 
#' 
#' @export

KL = function(y1, y2, width){
  x = seq(min(c(min(y1), min(y2))), max(c(max(y1), max(y2))), length.out = 750)
  p = hist(y1, x, plot = F)$counts
  q = hist(y2, x, plot = F)$counts
  p = p/length(y1)
  q = q/length(y2)
  eps = 1e-16
  dist = sum(p*(log2(p+eps) - log2(q+eps)))
  return(dist)
}


#' KL test
#' 
#' `KLtest` Tests for the significance of KL divergence between 
#' G(V_c, V_c') and G_R matrices 
#' 
#' @param s G(Vc, Vc') matrix distribution to test
#' @param null vector to sample null distributions (singleton indices)
#' @param true_dist Singleton matrix
#' @param a Significant level alpha. Default is 0.05
#' @param width default is 0.001. Can take out if you don't know what it is
#' @return Returns `A` Boolean, 1 if KL-divergence significant (C and C' connected) 
#' and 0 if divergence not significant. 
#' @examples
#' 
#' @export

KLtest = function(s, null, true_dist, a = 0.05, width = 0.001){
  # All s, null, true are row vectors
  s_kl = KL(s, true_dist, width)
  t_kl = rep(NA, 1000)
  
  for (i in 1:1000){
    t = sample(null, length(s), TRUE)
    t_kl[i] = KL(t, true_dist, width)
  }
  
  P = sum(t_kl > s_kl) / length(t_kl)
  if (P < a){
    A = 1
  } else {
    A = 0
  }
  return(A)
}

