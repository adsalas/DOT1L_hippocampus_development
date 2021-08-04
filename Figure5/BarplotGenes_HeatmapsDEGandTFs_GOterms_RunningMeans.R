

#Load the requiered libraries.
library(dplyr)
library(RaceID)
library(ComplexHeatmap)
library(circlize)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

##############   Computing the DEGs between control and cKO samples per cluster ##############################

### Differential gene expression analysis per cluster between conditions
# Create a list to store the information of the cells belonging to the different clusters divided by condition.
cellsInClCond <- list()
for (i in (1:max(sc@cpart))){
  cellsin <- names(sc@cpart[sc@cpart == i])
  cellsCont <- cellsin[grep("^C", cellsin)]
  cellscKO <- cellsin[grep("^K", cellsin)]
  cellsInClCond[[i]] <- list(cellsCont, cellscKO)  
} 

#Compute the DEGs per cluster between conditions
DEG_all_Cond <- list()
for (i in 1:max(sc@cpart)){  
  if (length(cellsInClCond[[i]][[1]]) != 0 & length(cellsInClCond[[i]][[2]]) != 0) {
    A <- cellsInClCond[[i]][[1]]
    B <- cellsInClCond[[i]][[2]]
    DEG_all_Cond[[i]] <- diffexpnb(getfdata(sc, n=c(A,B)), A=A, B=B)
  }
  else
    DEG_all_Cond[[i]] <- "No DEGs computed"
}


# Filter the results for each cluster based on the p value threshold (p < 0.05) and store all the different outputs as a list.
# Define the cell types
celltypes <- c("Granule", "CA1pyr", "NSC1", "CA3pyr","IPC1","Subic1", "CajRet1", "NSC2cyc", "IPC2","CajRet2", "NSC3", "NSC4cyc", "Subic2", "INs", "RedBl" , "Pericytes", "Unk1", "CortHem", "OPCs", "CajRet3", "Unk2", "Unk3", "NSC5", "ChPlexEp", "CajRet4",  "Microglia", "Unk4")
DEG_all_CondFTD <- DEG_all_Cond

for (i in 1:length(DEG_all_Cond)) {
  
  if (length(DEG_all_Cond[[i]]) > 1) {
  
    DEG_all_CondFTD[[i]] <-  DEG_all_Cond[[i]]$res
    DEG_all_CondFTD[[i]]$gene <- rownames(DEG_all_CondFTD[[i]])
    DEG_all_CondFTD[[i]] <- filter(DEG_all_CondFTD[[i]], padj < 0.05)
    DEG_all_CondFTD[[i]] <- DEG_all_CondFTD[[i]][ , c(4,5,6,7,8)]
    }
}

names(DEG_all_CondFTD) <- celltypes


# Keep the clusters that have cells for each of the samples
DEG_all_CondFTD <- DEG_all_CondFTD[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,20,21,23,25)]



#############  Barplot of the number of DE genes per cluster between the control and the cKOs  #####################

barplot(unlist(lapply(DEG_all_CondFTD, nrow)), names.arg = names(DEG_all_CondFTD),  
        cex.names = 1, col = sc@fcol[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,20,21,23,25)],
        ylab = "number of DEGs", ylim = c(0,100), las = 2)



############################################ Plotting all the DEGs between conditions as a heatmap  ##################################
# Create a data frames with the required information 
DEG_DF <- DEG_all_CondFTD[[1]][c(2,5)]
for (i in 2:length(DEG_all_CondFTD)) {
    DEG_DF <- full_join(DEG_DF, DEG_all_CondFTD[[i]][c(2,5)], by= "gene", keep=T)
} 

# Remove the NAs from the data frame
for (i in 1:nrow(DEG_DF)) {
  for (j in 1:ncol(DEG_DF)) {
    if (is.na(DEG_DF[i,j])) {
      DEG_DF[i,j] <- 0
    }
  }
}

# Simplify the data frame
# Set rownames
rownames(DEG_DF) <- DEG_DF$gene
# Remove gene column
DEG_DF <- DEG_DF[ , -2]
# Set the colnames
celltypes2 <- c("Granule", "CA1pyr", "NSC1", "CA3pyr","IPC1","Subic1", "CajRet1", "NSC2cyc", "IPC2","CajRet2", "NSC3", "NSC4cyc", "Subic2", "INs", "Pericytes", "Unk1", "CortHem", "CajRet3", "Unk2", "NSC5", "CajRet4")
colnames(DEG_DF) <- celltypes2

# Remove "Xist" gene
DEG_DF <- DEG_DF[!rownames(DEG_DF) == "Xist", ]

# Change the data frame to a matrix for the input
DEG_DF <- as.matrix(DEG_DF)

# Plot the heatmap
set.seed(1500)
Heatmap(DEG_DF, 
        show_row_names = F, 
        row_title = "genes", 
        column_title_side = "top", 
        show_heatmap_legend = T, 
        name = "log2FC", 
        col = colorRamp2(c(-7, 0, 7), c("blue", "white", "red")), 
        column_split = 2, 
        column_title = "%s", 
        border = T)



##################################################################################################################################
###########  Over-representation test for the genes increasing or decreasing in expression in each cluster in the cKOs ###########

# ORT using the clusterProfiler package 
# Granule cells (cluster 1)
Allgenes <- DEG_all_Cond[[1]]$res$log2FoldChange
names(Allgenes) <- rownames(DEG_all_Cond[[1]]$res)

genesClust1 <- DEG_all_CondFTD[[1]]$log2FoldChange
names(genesClust1) <- DEG_all_CondFTD[[1]]$gene

egoClust1 <- enrichGO(gene = names(genesClust1), universe = names(Allgenes), 
                      OrgDb = org.Mm.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05, minGSSize = 10, 
                      maxGSSize = 500, keyType = "SYMBOL")

genesClust1_up <- genesClust1[genesClust1 > 0]
egoClust1_up <- enrichGO(gene = names(genesClust1_up), universe = names(Allgenes), 
                         OrgDb = org.Mm.eg.db, ont = "BP", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, minGSSize = 10, 
                         maxGSSize = 500, keyType = "SYMBOL")

genesClust1_down <- genesClust1[genesClust1 < 0]
egoClust1_down <- enrichGO(gene = names(genesClust1_down), universe = names(Allgenes), 
                           OrgDb = org.Mm.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, minGSSize = 10, 
                           maxGSSize = 500, keyType = "SYMBOL")


# CA1pyr (cluster 2)
Allgenes <- DEG_all_Cond[[2]]$res$log2FoldChange
names(Allgenes) <- rownames(DEG_all_Cond[[2]]$res)

genesClust2 <- DEG_all_CondFTD[[2]]$log2FoldChange
names(genesClust2) <- DEG_all_CondFTD[[2]]$gene

egoClust2 <- enrichGO(gene = names(genesClust2), universe = names(Allgenes), 
                      OrgDb = org.Mm.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05, minGSSize = 10, 
                      maxGSSize = 500, keyType = "SYMBOL")

genesClust2_up <- genesClust2[genesClust2 > 0]
egoClust2_up <- enrichGO(gene = names(genesClust2_up), universe = names(Allgenes), 
                         OrgDb = org.Mm.eg.db, ont = "BP", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, minGSSize = 10, 
                         maxGSSize = 500, keyType = "SYMBOL")

genesClust2_down <- genesClust2[genesClust2 < 0]
egoClust2_down <- enrichGO(gene = names(genesClust2_down), universe = names(Allgenes), 
                           OrgDb = org.Mm.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, minGSSize = 10, 
                           maxGSSize = 500, keyType = "SYMBOL")


# CA3pyr (cluster 4)
Allgenes <- DEG_all_Cond[[4]]$res$log2FoldChange
names(Allgenes) <- rownames(DEG_all_Cond[[4]]$res)

genesClust4 <- DEG_all_CondFTD[[4]]$log2FoldChange
names(genesClust4) <- DEG_all_CondFTD[[4]]$gene

egoClust4 <- enrichGO(gene = names(genesClust4), universe = names(Allgenes), 
                      OrgDb = org.Mm.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05, minGSSize = 10, 
                      maxGSSize = 500, keyType = "SYMBOL")

genesClust4_up <- genesClust4[genesClust4 > 0]
egoClust4_up <- enrichGO(gene = names(genesClust4_up), universe = names(Allgenes), 
                         OrgDb = org.Mm.eg.db, ont = "BP", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, minGSSize = 10, 
                         maxGSSize = 500, keyType = "SYMBOL")

genesClust4_down <- genesClust4[genesClust4 < 0]
egoClust4_down <- enrichGO(gene = names(genesClust4_down), universe = names(Allgenes), 
                           OrgDb = org.Mm.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, minGSSize = 10, 
                           maxGSSize = 500, keyType = "SYMBOL")


# IPC1 (cluster 5)
Allgenes <- DEG_all_Cond[[5]]$res$log2FoldChange
names(Allgenes) <- rownames(DEG_all_Cond[[5]]$res)

genesClust5 <- DEG_all_CondFTD[[5]]$log2FoldChange
names(genesClust5) <- DEG_all_CondFTD[[5]]$gene

egoClust5 <- enrichGO(gene = names(genesClust5), universe = names(Allgenes), 
                      OrgDb = org.Mm.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05, minGSSize = 10, 
                      maxGSSize = 500, keyType = "SYMBOL")

genesClust5_up <- genesClust5[genesClust5 > 0]
egoClust5_up <- enrichGO(gene = names(genesClust5_up), universe = names(Allgenes), 
                         OrgDb = org.Mm.eg.db, ont = "BP", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, minGSSize = 10, 
                         maxGSSize = 500, keyType = "SYMBOL")

genesClust5_down <- genesClust5[genesClust5 < 0]
egoClust5_down <- enrichGO(gene = names(genesClust5_down), universe = names(Allgenes), 
                           OrgDb = org.Mm.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, minGSSize = 10, 
                           maxGSSize = 500, keyType = "SYMBOL")

# IPC2 (cluster 9)
Allgenes <- DEG_all_Cond[[9]]$res$log2FoldChange
names(Allgenes) <- rownames(DEG_all_Cond[[9]]$res)

genesClust9 <- DEG_all_CondFTD[[9]]$log2FoldChange
names(genesClust9) <- DEG_all_CondFTD[[9]]$gene

egoClust9 <- enrichGO(gene = names(genesClust9), universe = names(Allgenes), 
                      OrgDb = org.Mm.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05, minGSSize = 10, 
                      maxGSSize = 500, keyType = "SYMBOL")

genesClust9_up <- genesClust9[genesClust9 > 0]
egoClust9_up <- enrichGO(gene = names(genesClust9_up), universe = names(Allgenes), 
                         OrgDb = org.Mm.eg.db, ont = "BP", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, minGSSize = 10, 
                         maxGSSize = 500, keyType = "SYMBOL")

genesClust9_down <- genesClust9[genesClust9 < 0]
egoClust9_down <- enrichGO(gene = names(genesClust9_down), universe = names(Allgenes), 
                           OrgDb = org.Mm.eg.db, ont = "BP", 
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, minGSSize = 10, 
                           maxGSSize = 500, keyType = "SYMBOL")




#########################################################################################################################
#########  Barplot for the 6 main GO terms related to increased or decreased genes ######################################

# CLuster 1
# Increased genes
data1 <- egoClust1_up@result[1:6, ]
data_temp <- egoClust1_down@result[1:6, ]
data_temp$Count <- data_temp$Count * -1

data1 <- rbind(data1, data_temp)
# give the different terms the position for the text inside the bars
data1$textPos <- if_else(data1$Count > 0, 1,0)


# CLuster 2
# Increased genes
data2 <- egoClust2_up@result[1:6, ]
data_temp <- egoClust2_down@result[1:6, ]
data_temp$Count <- data_temp$Count * -1

data2 <- rbind(data2, data_temp)
# give the different terms the position for the text inside the bars
data2$textPos <- if_else(data2$Count > 0, 1,0)
data1 <- rbind(data1, data2)

# CLuster 4
# Increased genes
data2 <- egoClust4_up@result[1:6, ]
data_temp <- egoClust4_down@result[1:6, ]
data_temp$Count <- data_temp$Count * -1

data2 <- rbind(data2, data_temp)
# give the different terms the position for the text inside the bars
data2$textPos <- if_else(data2$Count > 0, 1,0)
data1 <- rbind(data1, data2)

# CLuster 9
# Increased genes
data2 <- egoClust9_up@result[1:6, ]
data_temp <- egoClust9_down@result[1:6, ]
data_temp$Count <- data_temp$Count * -1

data2 <- rbind(data2, data_temp)
# give the different terms the position for the text inside the bars
data2$textPos <- if_else(data2$Count > 0, 1,0)
data1 <- rbind(data1, data2)
rownames(data1) <- NULL


data1$Clusters <- c(rep("Granule", 12), rep("CA1pyr", 12),rep("CA3pyr", 12), rep("IPC2", 12))
data1$Clusters <- as.factor(data1$Clusters)
data1$order <- c(rep(1:12, 4))

ggplot(data1, aes(x = order, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Description, y = 0), hjust= data1$textPos, vjust= 0.5, size= 3) +
  coord_flip()  + 
  theme_classic() +
  facet_grid(rows = data1$Clusters) +
  scale_fill_gradient2(low = "lightblue", mid = "mediumblue", high = "darkblue", midpoint = 0.03) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(-20, 20), labels = c(20,10,0,10,20)) 



########################################################################################################################
############################### ANALYSIS OF TRANSCRIPTION FACTORS AFFECTED PER CLUSTER   ###############################

# Loading the selected databases
# TRRUST
trrustDB <- read.table("~/Downloads/trrust_rawdata.mouse.tsv", header = F)
colnames(trrustDB) <- c("TF", "Target", "Interaction", "V4")
# TF and co-factors DB 
tf_interactions <- read.table("~/Downloads/mouse_tf_interactions.csv", header = T, sep = ",")
TFs_both_DBs <- c(as.vector(unique(tf_interactions$Symbol1)), as.vector(unique(trrustDB$TF)))
TFs_both_DBs <- unique(TFs_both_DBs)
# Animal TF DB (Source: http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download)
tf_AnimalTFDB <- read.table("~/Downloads/Mus_musculus_TF.txt", header = F, sep = "", skip = 1)
colnames(tf_AnimalTFDB) <- c("Species", "Symbol", "Ensembl", "Family", "Protein", "Entrez_ID")
TF_BigDB <- c(TFs_both_DBs, as.vector(tf_AnimalTFDB$Symbol))
TF_BigDB <- unique(TF_BigDB)



# Remove the CortHem cluster 
dim(DEG_DF)
DEG_DF1 <- DEG_DF[ , -17]
dim(DEG_DF1)

# PLotting just the hippocmpal clusters included in the lineage trajectory analysis (DPT)
data2 <- DEG_DF1[rownames(DEG_DF1) %in% TF_BigDB , -c(6,7,10,13,14,15,16,17,18,20) ]
dim(data2)

data2 <- data2[ rowSums(data2) != 0, ]
dim(data2)
# Check min and max for defining the scale for the heatmap. parameter "col"
min(data2)
max(data2)

# Plot the heatmap
set.seed(1500)
Heatmap(data2, show_row_names = T, row_title = "genes", column_title_side = "top", 
        show_heatmap_legend = T, name = "log2FC", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")), 
        column_split = 4, 
        clustering_distance_columns = "pearson",
        column_labels = colnames(data2), 
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        border = T)


##########################################################################################################################
########   PLOTTING THE MAIN DE TFs BETWEEN CONTROL AND cKOs ALONG THE PSEUDOTIME ORDERING ###############################
##########################################################################################################################

library(ggplot2)

# Using the full data set 
data_subset <- as.data.frame(t(as.data.frame(fs[c("Tcf4", "Nr2f1", "Nfix", "Nfia", "Nfib", "Lhx2", "Zbtb20"), ])))
data_subset$cell <- rownames(data_subset)

data_subset$x <- c(1:length(data_subset$cell))
data_subset$sample2 <- substr(rownames(data_subset), start = 1, stop = 3)


rownames(data_colours) <- data_colours$cells_names
dat_temp <- data_colours[rownames(data_subset), ]
data_subset$colour <- dat_temp$colour
data_subset$lines <- "blue"
data_subset[data_subset$sample2 == "K81", ]$lines <- "red"
data_subset$colour <- as.character(data_subset$colour)

# Plotting the desired genes
# Confidence interval shown in black correspond to 0.95
ggplot(data_subset) + 
  geom_point(aes(x = x, y = Lhx2) , color= data_subset$colour, cex= 0.8, alpha= 0.5) +
  #geom_smooth(method = "loess", span= 0.35, se = TRUE, orientation = "x") +
  labs(title = "Lhx2", x= element_blank(), y= "Expression") +
  theme_classic() +
  geom_smooth(method = "loess", span= 0.35, size= 0.8, se = TRUE, fill= "black", show.legend = FALSE, orientation = "x",mapping = aes(x= x, y = Lhx2, colour= sample2))

ggplot(data_subset) + 
  geom_point(aes(x = x, y = Tcf4) , color= data_subset$colour, cex= 0.8, alpha= 0.5) +
  #geom_smooth(method = "loess", span= 0.35, se = TRUE, orientation = "x") +
  labs(title = "Tcf4", x= element_blank(), y= "Expression") +
  theme_classic() +
  geom_smooth(method = "loess", span= 0.35, size= 0.8, se = TRUE, fill= "black", show.legend = FALSE, orientation = "x",mapping = aes(x= x, y = Tcf4, colour= sample2))

ggplot(data_subset) + 
  geom_point(aes(x = x, y = Nr2f1) , color= data_subset$colour, cex= 0.8, alpha= 0.5) +
  #geom_smooth(method = "loess", span= 0.35, se = TRUE, orientation = "x") +
  labs(title = "Nr2f1", x= element_blank(), y= "Expression") +
  theme_classic() +
  geom_smooth(method = "loess", span= 0.35, size= 0.8, se = TRUE, fill= "black", show.legend = FALSE, orientation = "x",mapping = aes(x= x, y = Nr2f1, colour= sample2))

ggplot(data_subset) + 
  geom_point(aes(x = x, y = Nfix) , color= data_subset$colour, cex= 0.8, alpha= 0.5) +
  #geom_smooth(method = "loess", span= 0.35, se = TRUE, orientation = "x") +
  labs(title = "Nfix", x= element_blank(), y= "Expression") +
  theme_classic() +
  geom_smooth(method = "loess", span= 0.35, size= 0.8, se = TRUE, fill= "black", show.legend = FALSE, orientation = "x",mapping = aes(x= x, y = Nfix, colour= sample2))

ggplot(data_subset) + 
  geom_point(aes(x = x, y = Nfib) , color= data_subset$colour, cex= 0.8, alpha= 0.5) +
  #geom_smooth(method = "loess", span= 0.35, se = TRUE, orientation = "x") +
  labs(title = "Nfib", x= element_blank(), y= "Expression") +
  theme_classic() +
  geom_smooth(method = "loess", span= 0.35, size= 0.8, se = TRUE, fill= "black", show.legend = FALSE, orientation = "x",mapping = aes(x= x, y = Nfib, colour= sample2))

ggplot(data_subset) + 
  geom_point(aes(x = x, y = Nfia) , color= data_subset$colour, cex= 0.8, alpha= 0.5) +
  #geom_smooth(method = "loess", span= 0.35, se = TRUE, orientation = "x") +
  labs(title = "Nfia", x= element_blank(), y= "Expression") +
  theme_classic() +
  geom_smooth(method = "loess", span= 0.35, size= 0.8, se = TRUE, fill= "black", show.legend = FALSE, orientation = "x",mapping = aes(x= x, y = Nfia, colour= sample2))

ggplot(data_subset) + 
  geom_point(aes(x = x, y = Zbtb20) , color= data_subset$colour, cex= 0.8, alpha= 0.5) +
  #geom_smooth(method = "loess", span= 0.35, se = TRUE, orientation = "x") +
  labs(title = "Zbtb20", x= element_blank(), y= "Expression") +
  theme_classic() +
  geom_smooth(method = "loess", span= 0.35, size= 0.8, se = TRUE, fill= "black", show.legend = FALSE, orientation = "x",mapping = aes(x= x, y = Zbtb20, colour= sample2))

