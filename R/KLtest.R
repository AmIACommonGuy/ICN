# s: distribution to test (dimension or vector? )
# null: vector to sample null distributions
# true: distribution to compare with
KLtest = function(s, null, true_dist, a, width){
  # All s, null, true are row vectors
  s_kl = KL.divergence(s,true_dist,width)
  t_kl = rep(NA, 1000)
  
  for (i in 1:1000){
    t = sample(null, length(s), true)
    t_kl[i] = KL.divergence(t, true_dist, width)
  }
  
  P = sum(t_kl > s_kl) / length(t_kl)
  if (P < a){
    A = 1
  } else {
    A = 0
  }
}

