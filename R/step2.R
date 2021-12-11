# find pairwise connections
# input NICE result object 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "data.rds"
}

# default input is clust.rds; perform post processing
NICE = readRDS(clust.rds)
NICER = NICER(NICE)
data = readRDS(data.rds)
n = dim(data)[1]

# find number of communities
K = NICER$K
K.ind = NICER$K.ind

# find community indices
comm.ind = vector(mode = "list", length = K)
for (ind in 1:K){
  comm.ind[[ind]] = as.vector(which(NICE$Cindx == K.ind[ind]))
}

# create a vector of singletons
ordered.ind = NICE$Cindx[NICE$Clist]
num.comm.ind = max(which(ordered.ind %in% K.ind))
# num.singles = dim(data)[1] - num.comm.ind
singles = NICE$Clist[(num.comm.ind + 1):n]

# create the singleton matrix
G.R = data[singles, singles]

# find all pairwise communities; create G(Vc, Vc')
pairs = lst_pairs(K)
n.pairs = dim(to_matrix(pairs))[1]
G.pairwise = vector(mode = "list", length = n.pairs)
for (p in 1:n.pairs){
  pairs = to_matrix(pairs)
  nodes.comm1 = comm.ind[[pairs[p,1]]]
  nodes.comm2 = comm.ind[[pairs[p,2]]]
  nodes = c(nodes.comm1, nodes.comm2)
  G.pairwise[[p]] = data[nodes, nodes]
}

# call KLtest function
connections = lapply(G.pairwise, KLtest, null = singles,
                     true_dist = G.R, a = 0.05)

