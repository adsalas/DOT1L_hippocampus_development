# #######   Heatmap for the NSC and IPC selected markers  #################


# Define the cluster to keep
clustersTokeep <- c(5,9,3,8,12,11,18,23)
clustersVector <- vector()
for (i in 1:length(sc@cpart)) {
  if (sc@cpart[i] %in% clustersTokeep) {
    clustersVector <- c(clustersVector, sc@cpart[i])
  }
}

clusters <- as.data.frame(sc@cpart)
clusters$cell <- rownames(clusters)
colnames(clusters) <- c("cluster", "cell")
clusters <- clusters[rownames(clusters) %in% names(clustersVector), ]
library(dplyr)
clusters <- arrange_all(clusters, by_group = clusters)
rownames(clusters) <- clusters$cell

# Populate the data frame with the cluster names to which each cell belongs
for (i in 1:length(clusters$cluster)) {
  if (clusters$cluster[i] == 3) {
    clusters$cluster[i] <- "NSC1"
  }
  if (clusters$cluster[i] == 8) {
    clusters$cluster[i] <- "NSC2cyc"
  }
  if (clusters$cluster[i] == 11) {
    clusters$cluster[i] <- "NSC3"
  }
  if (clusters$cluster[i] == 12) {
    clusters$cluster[i] <- "NSC4cyc"
  }
  if (clusters$cluster[i] == 18) {
    clusters$cluster[i] <- "CortHem"
  }
  if (clusters$cluster[i] == 23) {
    clusters$cluster[i] <- "NSC5"
  }
  if (clusters$cluster[i] == 5) {
    clusters$cluster[i] <- "IPC1"
  }
  if (clusters$cluster[i] == 9) {
    clusters$cluster[i] <- "IPC2"
  }
}

library(circlize)
library(ComplexHeatmap)

# 1. Define the markers:
markers <- c("Wnt8b", "Cybrd1", "Sned1", "Adamts19", "Sox21", "Hes5", "Tac2")

# Subset the data
data3 <- as.matrix(sc@ndata)
dim(data3)
data3 <- data3[markers , ]
dim(data3)

data3_zsc <- t(as.data.frame((scale(t(data3)))))


column_ha2 = HeatmapAnnotation(cluster= clusters$cluster, col = list(cluster = c("NSC1" = sc@fcol[3], "NSC2cyc" = sc@fcol[8],
                                                                                 "NSC3" = sc@fcol[11], "NSC4cyc" = sc@fcol[12],
                                                                                 "CortHem" = sc@fcol[18], "NSC5" = sc@fcol[23],
                                                                                 "IPC1" = sc@fcol[5], "IPC2" = sc@fcol[9])))
# Plot the heatmap
set.seed(1500)
Heatmap(data3_zsc[ , rownames(clusters)], 
        show_row_names = T, 
        row_title = NULL, 
        column_title_side = "top", 
        show_heatmap_legend = T, 
        name = "zscore", 
        top_annotation = column_ha2,
        cluster_columns = F, 
        row_names_gp = gpar(fontsize = 11),
        border = T,
        show_column_names = F, 
        column_split =  clusters$cluster,
        cluster_column_slices = T, 
        show_row_dend = F,
        show_column_dend = T, 
        column_title_gp = gpar(fontsize = 11))

