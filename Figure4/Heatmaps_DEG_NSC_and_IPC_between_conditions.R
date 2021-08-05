##########   Heatmaps for DEGs between conditions for the NSC and IPC clusters  ###############################

# The heatmaps can be generated using as input the data frame (called "DEG_DF") created in the script from "Figure 5". 

# Load the libraries
library(ComplexHeatmap)
library(circlize)


# Subset the data for plotting the NSC clusters
set.seed(1500)
DEG_DF_NSCs <- DEG_DF[ , c(3,8,11,12, 20)]
DEG_DF_NSCs <- DEG_DF_NSCs[ rowSums(DEG_DF_NSCs) != 0, ]
Heatmap(DEG_DF_NSCs, show_row_names = T, row_title = "genes", column_title_side = "top", show_heatmap_legend = T, name = "log2FC", col = colorRamp2(c(-7, 0, 7), c("blue", "white", "red")), column_split = 2, column_title = "%s", border = T)

# Subset the data for plotting the IPC clusters
DEG_DF_IPCs <- DEG_DF[ , c(5,9)]
DEG_DF_IPCs <- DEG_DF_IPCs[ rowSums(DEG_DF_IPCs) != 0, ]
Heatmap(DEG_DF_IPCs, show_row_names = T, row_title = "genes", column_title_side = "top", show_heatmap_legend = T, name = "log2FC", col = colorRamp2(c(-7, 0, 7), c("blue", "white", "red")), border = T, row_names_gp = gpar(fontsize = 9), show_column_dend = F)
