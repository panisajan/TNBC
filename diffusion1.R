
## ppi construct from all genes, STRING combined score < 0.4 is filtered
ppi = read.csv('C:/Users/Asus/Desktop/AomamxAsus/breast cancer/ppi/edge.csv',header = T,sep = ",")
ppi_genename = ppi[,c(3,4,15)]

g = graph_from_data_frame(ppi_genename, directed = FALSE, vertices = NULL)
g <- set_edge_attr(g, "weight", value= ppi_genename$stringdb..score)

g3 = g
# K <- regularisedLaplacianKernel(g3,sigma2=1, add_diag=1,normalized=FALSE)
# save(K, file = "K_laplacian.Rdata")

load('K_laplacian.Rdata') ## variable name K
y_7 = ifelse(V(g3)$name %in% c('NOL4','STAR','C8G','NEIL1','SLC46A3','FRMD6','SCARF2'),1,0)

## diffusion score

F_all <- K %*% y_7
F_all = as.data.frame(F_all)
F_all$name = rownames(F_all)
# Sort the dataframe based on scores
sorted_scores_df <- F_all[order(F_all$V1,decreasing = TRUE), ]
# Calculate differences between adjacent scores
sorted_scores_df$differences <- c(NA, diff(sorted_scores_df$V1))
fil = subset(sorted_scores_df,sorted_scores_df$V1 > quantile(F_all$V1, 0.99))