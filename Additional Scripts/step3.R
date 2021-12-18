# find interconnected edges within interconnected communities
# input Communities dataframe of interconnection detected result which is 
# 'inRes' in step 2.

#!/usr/bin/env Rscript
inRes = commandArgs(trailingOnly=TRUE)
interCon = inRes[inRes$whether_interconnected == 1,]
#kID1 = NICE$CID[interCon$community_one[1]]
#kID2 = NICE$CID[interCon$community_two[1]]
kID1 = NICE$CID[4] # an example, 4 here means the 4th cluster in heatmap.
kID2 = NICE$CID[5]
nodes1 = which(NICE$Cindx == kID1)
nodes2 = which(NICE$Cindx == kID2)
nodecom = c(nodes1, nodes2)
mat1 = data[nodes1, nodes1]
mat2 = data[nodes2, nodes2]
mat12 = data[nodecom, nodecom]
inter12 = data[nodes1, nodes2]

result3 = InterCon(mat1, mat2, mat12, inter12, 
                   0.6, c(0.3, 0.5, 0.6, 0.8, 0.85, 0.9), nodes1, nodes2)
print(result3)