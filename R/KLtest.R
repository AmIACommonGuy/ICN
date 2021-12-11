#' KL test
#' 
#' `KLtest` Calculates the KL divergence between G(V_c, V_c') and G_R
#' matrices 
#' 
#' @param s G(Vc, Vc') matrix distribution to test
#' @param null vector to sample null distributions (singleton indices)
#' @param true_dist Singleton matrix
#' @param a Significant level alpha. Default is 0.05
#' @param width maxinum number of nearest neighbors to search (k)
#' @return Returns `A` Boolean, 1 if KL-divergence significant (C and C' connected) 
#' and 0 if divergence not significant. 
#' @examples
#' 
#' @export

KLtest = function(s, null, true_dist, a = 0.05, width){
  # All s, null, true are row vectors
  s_kl = KL.divergence(s, true_dist) #, width)
  t_kl = rep(NA, 1000)
  
  for (i in 1:1000){
    t = sample(null, dim(s)[2], TRUE)
    t_kl[i] = KL.divergence(t, true_dist[t, t], width)
  }
  
  P = sum(t_kl > s_kl) / length(t_kl)
  if (P < a){
    A = 1
  } else {
    A = 0
  }
  return(A)
}

