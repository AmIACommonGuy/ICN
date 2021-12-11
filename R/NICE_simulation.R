### Data Simulation  ####
library(matlab)
library(gplots)
library(rARPACK)
## Simulate Two Communities (100 nodes)

# Simulate true corr matrix
true1 = matrix(0, nrow = 100, ncol = 100)
true1[1:20,1:20] = 0.5
true1[21:30,21:30] = 0.5
true1[1:30,1:30] = true1[1:30,1:30] + 0.5*diag(30)
true1[31:100,31:100] = true1[31:100,31:100] + diag(70);
heatmap.2(true1, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", breaks=seq(0,1,0.01),
          density.info = "none")

# Simulate samples based on true corr matrix, and then calculate the sample corr matrix.
nvars = 100
nobs = 50
R = chol(true1) # upper triangular
Z_vec = rnorm(nvars * nobs)
Z = matrix(Z_vec, nrow = 100)
samples =  crossprod(R,Z) # can use crossprod
samples = t(samples) # For each node (100 nodes in total), draw 50 samples.
simu1 = cor(samples) # Get corr between nodes
heatmap.2(simu1, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

# Shuffle the sample corr matrix
shuffle_index = sample(1:100)
simu1_shuf = simu1[shuffle_index, shuffle_index]
heatmap.2(simu1_shuf, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

result = NICE(simu1_shuf, 0.2)
system.time(NICE(simu1_shuf, 0.2))

# Run simulation again when encounter error
# When result$Clist has -1, code below will have error. -1 corresponds to zero row in W (Line 79).
simu1_reordered = simu1_shuf[result$Clist, result$Clist]

# Sample corr matrix under true order
heatmap.2(simu1, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

# Sample corr matrix odered by NICE, after shuffling the true order
heatmap.2(simu1_reordered, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100), 
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key", 
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

