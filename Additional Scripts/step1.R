# # execute NICE algorithm for interested correlation matrix
# input correlation matrix, self-determined threshold

#!/usr/bin/env Rscript
corr = commandArgs(trailingOnly=TRUE)
NICE = NICE_fast_cut(corr, threshold = 0.2, cutter = 1)
save(NICE,file = "data/cluster.rds")
save(corr,file = "data/shuffle.rds")
