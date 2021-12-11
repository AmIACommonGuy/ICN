# find pairwise connections
# input NICE result object 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

# default input is clust.rds
nice.res = readRDS(clust.rds)
K = nice.res$K_selected

# find community indices
comm.ind = vector(mode = "list", length = K)
comm.ind = 

# find all pairwise communities: 
pairs = ICN::lst_pairs(K)
n.conns = dim(pairs)[1]

# call KLtest function
connections = KLtest(s, null, true_dist, a = 0.05, width)
