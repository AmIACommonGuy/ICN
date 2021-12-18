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
NICE = readRDS("data/cluster.rds")
NICER = NICER(NICE)
data = readRDS("data/shuffle.rds")
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
singles = NICE$Clist[(num.comm.ind + 1):n]   # equivalent to idx_left

# create the singleton matrix
G.R = data[singles, singles]   # equivalent to clu00
diag(G.R) = 0
true_dist = squareform(G.R)    # equivalent to true_dist. ha. 

# find all pairwise communities; create G(Vc, Vc')
pairs = lst_pairs(K)
n.pairs = dim(to_matrix(pairs))[1]
G.pairwise = vector(mode = "list", length = n.pairs)
null = vector()
for (p in 1:n.pairs){
  pairs = to_matrix(pairs)
  nodes.comm1 = comm.ind[[pairs[p,1]]]
  nodes.comm2 = comm.ind[[pairs[p,2]]]
  # nodes = c(nodes.comm1, nodes.comm2)
  CCp = data[nodes.comm1, nodes.comm2] # equivalent to CC
  G.pairwise[[p]] = as.vector(CCp)     # equivalent to Off{i,j}
  null = c(null, G.pairwise[[p]])
}

# write null distribution: 
Off_2 = data[NICE$Clist[1:num.comm.ind], singles]
null = c(null, as.vector(Off_2))

# call KLtest function
connections = lapply(G.pairwise, KLtest, null = null,
                     true_dist = true_dist, a = 0.05)

# transform list interconnection into a vector
interResult = unlist(connections)
# transform matrix of community pairs index into dataframe
inRes = as.data.frame(pairs)
colnames(inRes) = c('community_one', 'community_two')
inRes$whether_interconnected = interResult

# plot adjacency matrix

d = max(inRes$community_two)
A.mat = matrix(0, nrow = d, ncol = d)
diag(A.mat) = 1
for(i in c(1:nrow(inRes))){
  A.mat[inRes$community_one[i], inRes$community_two[i]] = 
    inRes$whether_interconnected[i]
  A.mat[inRes$community_two[i], inRes$community_one[i]] = 
    inRes$whether_interconnected[i]
}

par(cex.main=0.8, par(mar = c(1, 2, 2, 6)))
colors = rep(brewer.pal(9, "Blues"), each = 3)[1:17]
plot.A.mat = A.mat + diag(rep(2, nrow(A.mat)))    # emphasize diagonal in plot
adj.mat = heatmap.2(A.mat, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', 
                    key = F, col = as.vector(colors), main = "Communities Interconnection", 
                    colsep = 1:n, rowsep = 1:n, sepcolor = "grey", sepwidth=c(0.005,0.005))
image.plot(legend.only=T, zlim=range(0,1), col=colors, legend.mar = 3, 
           legend.shrink = 0.5, legend.width = 0.5, legend.cex = 0.3)
