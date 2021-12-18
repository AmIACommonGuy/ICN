# R script for applycation on TCGA

library(gplots)
library(matlab)
library(rARPACK)
library(parallel)
library(rARPACK)
library(stats)

# --------------------------------------------
#              Data Processing
# --------------------------------------------

# Data selecting
names(sort(table(cancerTypesVec), decreasing = TRUE)[1:17])
cancerTypesVec%in%names(sort(table(cancerTypesVec), decreasing = TRUE)[1:17])
cancerTypesVec_use = cancerTypesVec[cancerTypesVec%in%names(sort(table(cancerTypesVec), decreasing = TRUE)[1:17])]
allExprData2 = allExprData[,cancerTypesVec%in%names(sort(table(cancerTypesVec), decreasing = TRUE)[1:17])]

# Random sample 1000 cancer pacients
set.seed(1)
idx = sample(1:8214, 1000)
allExprData2_use = allExprData2[ ,idx]
cancerTypesVec_use = cancerTypesVec_use[idx]

# 17 Cancers, 1,000 Patients
cormat1000 = cor(allExprData2_use)
reordered1000 = cormat1000[NICE.apply$Clist, NICE.apply$Clist]
reordered1000.df = as.data.frame(reordered1000)

# Setting up 17 colors
load("~/biostat625/Project/17colors.Rdata")
Colors = cancerTypesVec_use[NICE.apply$Clist]
Colors[Colors == "BRCA"] = colors17[1]
Colors[Colors == "THCA"] = colors17[2]
Colors[Colors == "LUAD"] = colors17[3]
Colors[Colors == "LGG"] = colors17[4]
Colors[Colors == "CESC"] = colors17[5]
Colors[Colors == "OV"] = colors17[6]
Colors[Colors == "PRAD"] = colors17[7]
Colors[Colors == "KIRC"] = colors17[8]
Colors[Colors == "SKCM"] = colors17[9]
Colors[Colors == "LIHC"] = colors17[10]
Colors[Colors == "COAD"] = colors17[11]
Colors[Colors == "KIRP"] = colors17[12]
Colors[Colors == "HNSC"] = colors17[13]
Colors[Colors == "LUSC"] = colors17[14]
Colors[Colors == "SARC"] = colors17[15]
Colors[Colors == "BLCA"] = colors17[16]
Colors[Colors == "UCEC"] = colors17[17]

# --------------------------------------------
#              Making Graphs
# --------------------------------------------

# Graph heatmap for unordered correlation matrix
heatmap.2(cormat1000, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100),
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key",
          ColSideColors = Colors,
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

# Graph heatmap for ordered correlation matrix
heatmap.2(reordered1000, Rowv = FALSE, Colv = FALSE, margins = c(6,12), col = jet.colors(100),
          trace = "none", labRow = FALSE, labCol = FALSE, key.title = "Color Key",
          ColSideColors = Colors,
          keysize = 0.9, key.par = list(cex=0.5), key.xlab = "value", density.info = "none")

# Applying NICER
NICER(NICE.apply)

# Apply step 2 and 3





# Save unordered, reordered Correlation matrix and other data
save(Colors, allExprData2_use, cancerTypesVec_use, cormat1000, NICE.apply, reordered1000.df, file = "NICE_OUTPUT.RData")
