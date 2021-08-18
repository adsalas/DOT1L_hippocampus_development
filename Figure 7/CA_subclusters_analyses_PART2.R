############################################################################################################################################
################################### SUBCLUSTERING ANALYSES FOR THE CA PYRAMIDAL CELL POPULATIONS PART 2  ###################################
############################################################################################################################################



###########################################################################################################################
######################################################  CA1  SUBCLUSTERS   ################################################


#################  FISHER's EXACT TEST FOR CA SUBCLUSTERS CONTROL vs CKOs   #############

sc <- sc2_2cls
cl <- sc@cpart
df_cl <- data.frame(names=names(cl), cluster= cl)
df_cl$condition <- substr(df_cl$names, start = 0, stop = 4)
table_clusters <- table(df_cl$condition, df_cl$cluster)

# We are interested in comparing the 2 main subclusters (the other ones are small and come from the outlier estimation)
table_clusters <- table_clusters[ , 1:2]

table_clusters1 <- table_clusters[1:2, ]
table_clusters2 <- table_clusters[c(1,3), ]
table_clusters3 <- table_clusters[2:3, ]


dataLIST1 <- list()
data <- data.frame("CellsIn"= c(1,1), "CellsOut"= c(1,1))

for (j in 1:dim(table_clusters1)[2]){
  for (i in 1:dim(table_clusters1)[1]) {
    data$CellsIn[i] <- table_clusters1[i,j]
    data$CellsOut[i] <- (sum(table_clusters1[i,]) - table_clusters1[i,j])
    dataLIST1[[j]] <- data
    rownames(dataLIST1[[j]]) <- rownames(table_clusters1)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST1)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST1[[i]])
}


FISHER_RESULTS1 <- FISHER_RESULTS


for (j in 1:dim(table_clusters2)[2]){
  for (i in 1:dim(table_clusters2)[1]) {
    data$CellsIn[i] <- table_clusters2[i,j]
    data$CellsOut[i] <- (sum(table_clusters2[i,]) - table_clusters2[i,j])
    dataLIST1[[j]] <- data
    rownames(dataLIST1[[j]]) <- rownames(table_clusters2)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST1)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST1[[i]])
}


FISHER_RESULTS2 <- FISHER_RESULTS


for (j in 1:dim(table_clusters3)[2]){
  for (i in 1:dim(table_clusters3)[1]) {
    data$CellsIn[i] <- table_clusters3[i,j]
    data$CellsOut[i] <- (sum(table_clusters3[i,]) - table_clusters3[i,j])
    dataLIST1[[j]] <- data
    rownames(dataLIST1[[j]]) <- rownames(table_clusters3)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST1)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST1[[i]])
}


FISHER_RESULTS3 <- FISHER_RESULTS


###########   Correcting the p values with the p.adj() built-in function and Bonferroni  ##############

# C817 vs K813
pvalues1 <- list()
for (i in 1:length(FISHER_RESULTS1)) {
  pvalues1[i] <- FISHER_RESULTS1[[i]]$p.value
}
pvalues1 <- unlist(pvalues1)
# Adjusting the p value
pvalues1_adj <- p.adjust(pvalues1, method = "bonferroni")
# Checkign the padj that passed the threshold (<0.05)
pvalues1_adj < 0.05

# C817 vs K814
pvalues2 <- list()
for (i in 1:length(FISHER_RESULTS2)) {
  pvalues2[i] <- FISHER_RESULTS2[[i]]$p.value
}
pvalues2 <- unlist(pvalues2)
# Adjusting the p value
pvalues2_adj <- p.adjust(pvalues2, method = "bonferroni") 
# Checkign the padj that passed the threshold (<0.05)
pvalues2_adj < 0.05

# K813 vs K814
pvalues3 <- list()
for (i in 1:length(FISHER_RESULTS3)) {
  pvalues3[i] <- FISHER_RESULTS3[[i]]$p.value
}
pvalues3 <- unlist(pvalues3)
# Adjusting the p value
pvalues3_adj <- p.adjust(pvalues3, method = "bonferroni")
# Checkign the padj that passed the threshold (<0.05)
pvalues3_adj < 0.05

# Determine in which clusters 1 vs 1 comparisons are significant between control and cKOs
pv_CA1 <- data.frame(cont_vs_cKO1 = pvalues1_adj, 
                     cont_vs_cKO2 = pvalues2_adj, 
                     cKO1_vs_cKO2 = pvalues3_adj)
pv_CA1_adj_bonf <- data.frame(cont_vs_cKO1 = pvalues1_adj < 0.05, 
                              cont_vs_cKO2 = pvalues2_adj < 0.05, 
                              cKO1_vs_cKO2 = pvalues3_adj < 0.05)



#############    Barplots cell proportions in CA subclusters per condition ##############################################

### Making a barplot for showing percentage of cells between control and cKOs for CA1 subclusters
library(ggplot2)
sc <- sc2_2cls
cl <- sc@cpart
df_cl <- data.frame(names=names(cl), cluster= cl)
df_cl$condition <- substr(df_cl$names, start = 0, stop = 4)
table_clusters <- table(df_cl$condition, df_cl$cluster)

subcluster <- rep(c("Less mature","More mature","3"), 3)
condition <- c(rep("Control", 3), rep("KO1", 3), rep("KO2", 3))
value <- c(table_clusters[1,], table_clusters[2,], table_clusters[3,])
total <- c(rep(sum(table_clusters[1,]),3), rep(sum(table_clusters[2,]),3), rep(sum(table_clusters[3,]),3))
data <- data.frame(subcluster,condition,value, total)


for (i in 1:dim(data)[1]) {
  data$percent <- data$value/data$total * 100
}

# Plotting just the 2 subclusters of interest
ggplot(data[c(1,2,4,5,7,8), ], aes(fill=subcluster, y=percent, x=condition)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x= "Sample", y= "Cells per cluster (%)") +
  scale_fill_manual(values=c("#0000FFFF", "#00FF00FF")) +
  ylim(0,100) + 
  theme(legend.position = "bottom")





########### ANALYSIS OF DEGS BETWEEN CONTROL AND cKOs FOR THE LESS MATURE CA SUBCLUSTERS  #######################

# Define the cells in each subcluster based in the condition
sc <- sc2_2cls

cellsInClCond <- list()
for (i in (1:max(sc@cpart))){
  cellsin <- names(sc@cpart[sc@cpart == i])
  cellsCont <- cellsin[grep("^C", cellsin)]
  cellscKO <- cellsin[grep("^K", cellsin)] 
  cellsInClCond[[i]] <- list(cellsCont, cellscKO)  
} 


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

DEGCA1_subclusters_All <- DEG_all_Cond


DEGCA1sub1condFTD <- DEG_all_Cond[[1]]$res[DEG_all_Cond[[1]]$res$padj < 0.05 , ]
DEGCA1sub2condFTD <- DEG_all_Cond[[2]]$res[DEG_all_Cond[[2]]$res$padj < 0.05 , ]


################ Over-representation test  #########################################

###########  Subcluster 1

DEGCA1sub1condFTD$gene <- rownames(DEGCA1sub1condFTD)

genesCA1sub1_cond <- DEGCA1sub1condFTD$gene
cellsCA1_1_cont <- unlist(cellsInClCond[[1]][1])
cellsCA1_1_cKO <- unlist(cellsInClCond[[1]][2])
cellsCA1_1_cKO813 <- grep("^K813", cellsCA1_1_cKO, value = T)
cellsCA1_1_cKO814 <- grep("^K814", cellsCA1_1_cKO, value = T)
list2 <- list(cellsCA1_1_cont, cellsCA1_1_cKO813, cellsCA1_1_cKO814)

data2 <- sc@ndata[genesCA1sub1_cond, unlist(list2)]
data2 <- as.matrix(data2)
data2_zsc <- t(as.data.frame((scale(t(data2)))))

# Filter out "Xist" from the z-score transformed data frame
data2_zsc <- data2_zsc[rownames(data2_zsc) != "Xist", ]

# Creating and plotting a heatmpa of the DEGs among control and cKOs for the less mature subcluster
# Define the annotation
column_ha = HeatmapAnnotation(sample= c(rep("Control", length(cellsCA1_1_cont)), rep("KO1", length(cellsCA1_1_cKO813)), rep("KO2", length(cellsCA1_1_cKO814))), col = list(sample = c("Control" = "blue", "KO1" = "green", "KO2" = "red")))
# Filter out 'Xist'
DEGCA1sub1condFTD <- DEGCA1sub1condFTD[rownames(DEGCA1sub1condFTD) != "Xist", ]
# Define the groups based on the direction of the genes: Up or Down
DEGCA1sub1condFTD$clust <- 1
DEGCA1sub1condFTD[DEGCA1sub1condFTD$log2FoldChange > 0, ]$clust <- 2
group <- DEGCA1sub1condFTD$clust

# Plot the heatmap
set.seed(100)
Heatmap(data2_zsc, cluster_columns = T, cluster_rows = T, column_names_rot = T, column_dend_reorder = F, column_names_centered = T, name = "z-score", top_annotation = column_ha, show_column_names = FALSE, row_split = group, row_names_gp = gpar(fontsize = 10))

# Run the over-representation test
##### Genes Up
AllgenesLog2FC <- as.numeric(DEGCA1_subclusters_All[[1]]$res$log2FoldChange)
names(AllgenesLog2FC) <- rownames(DEGCA1_subclusters_All[[1]]$res)

CA1c1 <- as.numeric(DEGCA1sub1condFTD$log2FoldChange[DEGCA1sub1condFTD$clust == 1])
names(CA1c1) <- rownames(DEGCA1sub1condFTD)[DEGCA1sub1condFTD$clust == 1]

# GO terms by Molecular function
egoCA1c1cond <- enrichGO(gene = names(CA1c1), universe = names(AllgenesLog2FC), 
                         OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

##### Genes Down
CA1c2 <- as.numeric(DEGCA1sub1condFTD$log2FoldChange[DEGCA1sub1condFTD$clust == 2])
names(CA1c2) <- rownames(DEGCA1sub1condFTD)[DEGCA1sub1condFTD$clust == 2]

# GO terms by Molecular function
egoCA1c2cond <- enrichGO(gene = names(CA1c2), universe = names(AllgenesLog2FC), 
                         OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

# Combine Go terms that share similarity of at least "cutoff"
egoSimplifiedCA1c1 <- simplify(egoCA1c1cond, cutoff =0.85, by="p.adjust")
egoSimplifiedCA1c2 <- simplify(egoCA1c2cond, cutoff =0.85, by="p.adjust")


#### Barplots
# Increased genes cKOs
data1 <- egoSimplifiedCA1c2@result[1:5, ]
data1$order <- c(1:length(data1$ID))
ggplot(data1, aes(x = order, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Description, y = 0.1), hjust=0, vjust= 0.5, size= 4) +
  coord_flip()  + 
  theme_classic() +
  scale_fill_gradient(low = "orange", high = "royalblue") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0,0)) 

# Decreased genes cKOs
data1 <- egoSimplifiedCA1c1@result[1:8, ]
data1$order <- c(1:length(data1$ID))
ggplot(data1, aes(x = order, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Description, y = 0.1), hjust=0, vjust= 0.5, size= 4) +
  coord_flip()  + 
  theme_classic() +
  scale_fill_gradient(low = "orange", high = "royalblue") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0,0)) 



########################################################################################################################
############################### ANALYSIS OF TRANSCRIPTION FACTORS    ###################################################

# Loading the selected databases
# TRRUST
trrustDB <- read.table("~/Downloads/trrust_rawdata.mouse.tsv", header = F)
colnames(trrustDB) <- c("TF", "Target", "Interaction", "V4")
# TF and co-factors DB 
tf_interactions <- read.table("~/Downloads/mouse_tf_interactions.csv", header = T, sep = ",")
# Animal TF DB (Source: http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download)
tf_AnimalTFDB <- read.table("~/Downloads/Mus_musculus_TF.txt", header = F, sep = "", skip = 1)
colnames(tf_AnimalTFDB) <- c("Species", "Symbol", "Ensembl", "Family", "Protein", "Entrez_ID")

# Combine the databases
TF_BigDB <- c(as.vector(unique(tf_interactions$Symbol1)), as.vector(unique(trrustDB$TF)), as.vector(tf_AnimalTFDB$Symbol))
TF_BigDB <- unique(TF_BigDB)


#####  VENN DIAGRAM #################################

# Check the TFs that are present in the DEGs between subcluster 1 and 2
TFs_sub1vs2 <- DEGCA1_1vs2FTD[DEGCA1_1vs2FTD$gene %in% TF_BigDB, c(8,5,7)]

# Check all the TFs present in the DEGs between control and cKO for subcluster 1
TFs_DE_sub1_cond <- DEGCA1sub1condFTD[DEGCA1sub1condFTD$gene %in% TF_BigDB, c(8,5,7)]

venn.diagram(x= list(TFs_sub1vs2$gene, TFs_DE_sub1_cond$gene),
             filename = "venn_CA1_TFs.png",
             category.names = c(" " , " "),
             output =T, 
             lwd = 2,
             lty = 'blank',
             fill= c("#FF6347","#FFD700"),
             cex = 3,
             fontfamily = "sans-serif")



# Check which of those TFs are DE between control and cKO in the subcluster 1
TFs_putative_reg <- DEGCA1sub1condFTD[DEGCA1sub1condFTD$gene %in% TFs_sub1vs2$gene, c(8,5,7)]
CA1_TFs_sub1vs2 <- TFs_sub1vs2


# Make a barplot for TFs DE between CA1 subcluster 1 and subcluster 2
ggplot(CA1_TFs_sub1vs2, aes(x= log2FoldChange, y= gene, fill= padj)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12))

# Make a barplot for TFs DE between cKO and control in the subcluster 1
ggplot(TFs_putative_reg, aes(x= log2FoldChange, y= gene, fill= padj)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12))



###########################################################################################################################
######################################################  CA3  SUBCLUSTERS   ################################################

#################  FISHER's EXACT TEST FOR CA SUBCLUSTERS CONTROL vs CKOs   #############

sc <- sc4_2cls
cl <- sc@cpart
df_cl <- data.frame(names=names(cl), cluster= cl)
df_cl$condition <- substr(df_cl$names, start = 0, stop = 4)
table_clusters <- table(df_cl$condition, df_cl$cluster)

# We are interested in comparing the 2 main subclusters (the other ones are small and come from the outlier estimation)
table_clusters <- table_clusters[ , 1:2]

table_clusters1 <- table_clusters[1:2, ]
table_clusters2 <- table_clusters[c(1,3), ]
table_clusters3 <- table_clusters[2:3, ]


dataLIST1 <- list()
data <- data.frame("CellsIn"= c(1,1), "CellsOut"= c(1,1))

for (j in 1:dim(table_clusters1)[2]){
  for (i in 1:dim(table_clusters1)[1]) {
    data$CellsIn[i] <- table_clusters1[i,j]
    data$CellsOut[i] <- (sum(table_clusters1[i,]) - table_clusters1[i,j])
    dataLIST1[[j]] <- data
    rownames(dataLIST1[[j]]) <- rownames(table_clusters1)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST1)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST1[[i]])
}


FISHER_RESULTS1 <- FISHER_RESULTS


for (j in 1:dim(table_clusters2)[2]){
  for (i in 1:dim(table_clusters2)[1]) {
    data$CellsIn[i] <- table_clusters2[i,j]
    data$CellsOut[i] <- (sum(table_clusters2[i,]) - table_clusters2[i,j])
    dataLIST1[[j]] <- data
    rownames(dataLIST1[[j]]) <- rownames(table_clusters2)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST1)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST1[[i]])
}


FISHER_RESULTS2 <- FISHER_RESULTS


for (j in 1:dim(table_clusters3)[2]){
  for (i in 1:dim(table_clusters3)[1]) {
    data$CellsIn[i] <- table_clusters3[i,j]
    data$CellsOut[i] <- (sum(table_clusters3[i,]) - table_clusters3[i,j])
    dataLIST1[[j]] <- data
    rownames(dataLIST1[[j]]) <- rownames(table_clusters3)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST1)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST1[[i]])
}


FISHER_RESULTS3 <- FISHER_RESULTS


###########   Correcting the p values with the p.adj() built-in function and Bonferroni  ##############

# C817 vs K813
pvalues1 <- list()
for (i in 1:length(FISHER_RESULTS1)) {
  pvalues1[i] <- FISHER_RESULTS1[[i]]$p.value
}
pvalues1 <- unlist(pvalues1)
# Adjusting the p value
pvalues1_adj <- p.adjust(pvalues1, method = "bonferroni")
# Checkign the padj that passed the threshold (<0.05)
pvalues1_adj < 0.05

# C817 vs K814
pvalues2 <- list()
for (i in 1:length(FISHER_RESULTS2)) {
  pvalues2[i] <- FISHER_RESULTS2[[i]]$p.value
}
pvalues2 <- unlist(pvalues2)
# Adjusting the p value
pvalues2_adj <- p.adjust(pvalues2, method = "bonferroni")
# Checkign the padj that passed the threshold (<0.05)
pvalues2_adj < 0.05

# K813 vs K814
pvalues3 <- list()
for (i in 1:length(FISHER_RESULTS3)) {
  pvalues3[i] <- FISHER_RESULTS3[[i]]$p.value
}
pvalues3 <- unlist(pvalues3)
# Adjusting the p value
pvalues3_adj <- p.adjust(pvalues3, method = "bonferroni")
# Checkign the padj that passed the threshold (<0.05)
pvalues3_adj < 0.05

# Determine in which clusters 1 vs 1 comparisons are significant between control and cKOs
pv_CA3 <- data.frame(cont_vs_cKO1 = pvalues1_adj, 
                     cont_vs_cKO2 = pvalues2_adj, 
                     cKO1_vs_cKO2 = pvalues3_adj)
pv_CA3_adj_bonf <- data.frame(cont_vs_cKO1 = pvalues1_adj < 0.05, 
                              cont_vs_cKO2 = pvalues2_adj < 0.05, 
                              cKO1_vs_cKO2 = pvalues3_adj < 0.05)



#############    Barplots cell proportions in CA subclusters per condition ##############################################

### Making a barplot for showing the percentage of cells between control and cKOs for CA3 subclusters
sc <- sc4_2cls
cl <- sc@cpart
df_cl <- data.frame(names=names(cl), cluster= cl)
df_cl$condition <- substr(df_cl$names, start = 0, stop = 4)
table_clusters <- table(df_cl$condition, df_cl$cluster)

subcluster <- rep(c("Less mature","More mature","3","4"), 3)
condition <- c(rep("Control", 4), rep("KO1", 4), rep("KO2", 4))
value <- c(table_clusters[1,], table_clusters[2,], table_clusters[3,])
total <- c(rep(sum(table_clusters[1,]),4), rep(sum(table_clusters[2,]),4), rep(sum(table_clusters[3,]),4))
data <- data.frame(subcluster,condition,value, total)


for (i in 1:dim(data)[1]) {
  data$percent <- data$value/data$total * 100
}


ggplot(data[c(1,2,5,6,9,10),], aes(fill=subcluster, y=percent, x=condition)) + 
  geom_bar(position = "dodge", stat = "identity") + theme_light() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x= "Sample", y= "Cells per cluster (%)") +
  scale_fill_manual(values=c("#FA8072FF", "#AFEEEEFF")) + 
  ylim(0,100) + 
  theme(legend.position = "bottom")



########### ANALYSIS OF DEGS BETWEEN CONTROL AND cKOs FOR THE LESS MATURE CA CLUSTERS AND SELECTION OF THE TRANSCRIPTION FACTORS THAT ARE PUTATIVE CANDIDATES   #######################

# Define the cells in each subcluster based in the condition
sc <- sc4_2cls

cellsInClCond <- list()
for (i in (1:max(sc@cpart))){
  cellsin <- names(sc@cpart[sc@cpart == i])
  cellsCont <- cellsin[grep("^C", cellsin)]
  cellscKO <- cellsin[grep("^K", cellsin)] 
  cellsInClCond[[i]] <- list(cellsCont, cellscKO)  
}


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

DEGCA3_subclusters_All <- DEG_all_Cond

DEGCA3sub1condFTD <- DEG_all_Cond[[1]]$res[DEG_all_Cond[[1]]$res$padj < 0.05 , ]
DEGCA3sub2condFTD <- DEG_all_Cond[[2]]$res[DEG_all_Cond[[2]]$res$padj < 0.05 , ]

################ Over-representation test  #########################################


###########  Subcluster 1

DEGCA3sub1condFTD$gene <- rownames(DEGCA3sub1condFTD)

genesCA3sub1_cond <- DEGCA3sub1condFTD$gene
cellsCA3_1_cont <- unlist(cellsInClCond[[1]][1])
cellsCA3_1_cKO <- unlist(cellsInClCond[[1]][2])
cellsCA3_1_cKO813 <- grep("^K813", cellsCA3_1_cKO, value = T)
cellsCA3_1_cKO814 <- grep("^K814", cellsCA3_1_cKO, value = T)
list2 <- list(cellsCA3_1_cont, cellsCA3_1_cKO813, cellsCA3_1_cKO814)

data2 <- sc@ndata[genesCA3sub1_cond, unlist(list2)]
data2 <- as.matrix(data2)
data2_zsc <- t(as.data.frame((scale(t(data2)))))

# Filter out "Xist" from the z-score transformed data frame
data2_zsc <- data2_zsc[rownames(data2_zsc) != "Xist", ]

# Creating and plotting a heatmpa of the DEGs among control and cKOs for the less mature subcluster
# Define the annotation
column_ha = HeatmapAnnotation(sample= c(rep("Control", length(cellsCA3_1_cont)), rep("KO1", length(cellsCA3_1_cKO813)), rep("KO2", length(cellsCA3_1_cKO814))), col = list(sample = c("Control" = "blue", "KO1" = "green", "KO2" = "red")))
# Filter out 'Xist'
DEGCA3sub1condFTD <- DEGCA3sub1condFTD[rownames(DEGCA3sub1condFTD) != "Xist", ]
# Define the groups based on the direction of the genes: Up or Down
DEGCA3sub1condFTD$clust <- 1
DEGCA3sub1condFTD[DEGCA3sub1condFTD$log2FoldChange > 0, ]$clust <- 2
group <- DEGCA3sub1condFTD$clust

# Plot the heatmap
set.seed(100)
Heatmap(data2_zsc, cluster_columns = T, cluster_rows = T, column_names_rot = T, column_dend_reorder = F, column_names_centered = T, name = "z-score", top_annotation = column_ha, show_column_names = FALSE, row_split = group, row_names_gp = gpar(fontsize = 10))

# Run the over-representation test
##### Genes Up
AllgenesLog2FC <- as.numeric(DEGCA3_subclusters_All[[1]]$res$log2FoldChange)
names(AllgenesLog2FC) <- rownames(DEGCA3_subclusters_All[[1]]$res)

CA3c1 <- as.numeric(DEGCA3sub1condFTD$log2FoldChange[DEGCA3sub1condFTD$clust == 1])
names(CA3c1) <- rownames(DEGCA3sub1condFTD)[DEGCA3sub1condFTD$clust == 1]

# GO terms by Molecular function
egoCA3c1cond <- enrichGO(gene = names(CA3c1), universe = names(AllgenesLog2FC), 
                         OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

##### Genes Down

CA3c2 <- as.numeric(DEGCA3sub1condFTD$log2FoldChange[DEGCA3sub1condFTD$clust == 2])
names(CA3c2) <- rownames(DEGCA3sub1condFTD)[DEGCA3sub1condFTD$clust == 2]

# GO terms by Molecular function
egoCA3c2cond <- enrichGO(gene = names(CA3c2), universe = names(AllgenesLog2FC), 
                         OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")


### Combine Go terms that share similarity of at least "cutoff"
egoSimplifiedCA3c1 <- simplify(egoCA3c1cond, cutoff =0.85, by="p.adjust")
egoSimplifiedCA3c2 <- simplify(egoCA3c2cond, cutoff =0.85, by="p.adjust")


#### Barplots
# Increased genes cKOs
data1 <- egoSimplifiedCA3c2@result[1:8, ]
data1$order <- c(1:length(data1$ID))
ggplot(data1, aes(x = order, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Description, y = 0.1), hjust=0, vjust= 0.5, size= 4) +
  coord_flip()  + 
  theme_classic() +
  scale_fill_gradient(low = "orange", high = "royalblue") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0,0)) 

# Decreased genes cKOs
data1 <- egoSimplifiedCA3c1@result[1:8, ]
data1$order <- c(1:length(data1$ID))
ggplot(data1, aes(x = order, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Description, y = 0.1), hjust=0, vjust= 0.5, size= 4) +
  coord_flip()  + 
  theme_classic() +
  scale_fill_gradient(low = "orange", high = "royalblue") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0,0)) 


########################################################################################################################
############################### ANALYSIS OF TRANSCRIPTION FACTORS    ###################################################

# NOTE: the TFs databases were loaded in the sections above


#####  VENN DIAGRAM #################################

# Check the TFs that are present in the DEGs between subcluster 1 and 2
TFs_sub1vs2 <- DEGCA3_1vs2FTD[DEGCA3_1vs2FTD$gene %in% TF_BigDB, c(8,5,7)]

# Check all the TFs present in the DEGs between control and cKO for subcluster 1
TFs_DE_sub1_cond <- DEGCA3sub1condFTD[DEGCA3sub1condFTD$gene %in% TF_BigDB, c(8,5,7)]

venn.diagram(x= list(TFs_sub1vs2$gene, TFs_DE_sub1_cond$gene),
             filename = "venn_CA3_TFs.png",
             category.names = c(" " , " "),
             output =T, 
             lwd = 2,
             lty = 'blank',
             fill= c("#FF6347","#FFD700"),
             cex = 3,
             fontfamily = "sans-serif")


# Check which of those TFs are DE between control and cKO in the subcluster 1
TFs_putative_reg <- DEGCA3sub1condFTD[DEGCA3sub1condFTD$gene %in% TFs_sub1vs2$gene, c(8,5,7)]
CA3_TFs_sub1vs2 <- TFs_sub1vs2

# Make a barplot for TFs DE between CA3 subcluster 1 and subcluster 2
ggplot(CA3_TFs_sub1vs2, aes(x= log2FoldChange, y= gene, fill= padj)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12))

# Make a barplot for TFs DE between cKO and control in the subcluster 1
ggplot(TFs_putative_reg, aes(x= log2FoldChange, y= gene, fill= padj)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12))
