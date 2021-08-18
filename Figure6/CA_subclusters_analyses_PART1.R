#############################################################################################################################################
################################### SUBCLUSTERING ANALYSES FOR THE CA PYRAMIDAL CELL POPULATIONS PART 1   ###################################
#############################################################################################################################################

# Load the libraries
library(RaceID) # NOTE: the clustering has to be performed using RaceID version 0.1.3
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(VennDiagram)
library(ggplot2)
library(EnhancedVolcano)

###########################################################################################
################################### CA1 PYRAMIDAL CELLS ###################################
###########################################################################################

# Subset the data corresponding to the specific CA clusters
# CA1 cells
cells_cl2 <- names(sc@cpart[sc@cpart == 2])
prdata_cl2 <- prdata[ , which(names(prdata) %in% cells_cl2)]

############# PERFORMING THE SUBCLUSTERING ANALYSIS DEFINING 2 CLUSTERS AND LOGPEARSON AS METRIC ###############

# Filter the data based on total counts
cs <- colSums(prdata_cl2)
prdata_cl2 <- prdata_cl2[,cs > 500]

# Initialize SC object with prdata
sc2_2cls <- SCseq(prdata_cl2)

# Filter cells with at least 'mintotal' number of  total transcripts
sc2_2cls <- filterdata(sc2_2cls, 
                       mintotal     =  3000, 
                       minexpr      =  5, 
                       minnumber    =  1,
                       LBatch       =  NULL, 
                       knn          =  10, 
                       CGenes       =  c("Mki67", "Pcna"), 
                       FGenes       =  c("Kcnq1ot1"),
                       ccor         =  0.4,
                       bmode        =  "RaceID")

# Compute cell-to-cell distance matrix for clustering and outlier identification
sc2_2cls <- compdist(sc2_2cls, 
                     metric        =  "logpearson", 
                     FSelect       =  T, 
                     knn           =  NULL)

# Perform clustering
sc2_2cls <- clustexp(sc2_2cls, 
                     sat          =  TRUE, 
                     #samp         =  1000, 
                     cln          =  2, 
                     clustnr      =  30,
                     bootnr       =  50, 
                     rseed        =  17000, 
                     FUNcluster   =  "kmedoids")

# Compute the t-SNE map
sc2_2cls <- comptsne(sc2_2cls, initial_cmd = TRUE, perplexity = 30, rseed = 15555)

# k-nearest neighbour graph layout utilizing the Fruchterman-Rheingold algorithm
sc2_2cls <- compfr(sc2_2cls,knn=10)

# Find outliers in clusters
sc2_2cls <- findoutliers(sc2_2cls, 
                         probthr      =  1e-5, 
                         outminc      =  3, 
                         outlg        =  3,
                         outdistquant = 0.95)



############## MAKING SOME REPRESENTATIVE PLOTS ########################################

# Plot the t-SNE representation for the markers expression in the original sc object
# CA1 marker
plotexpmap(sc, "Pou3f1", logsc = T)

# Plot the t-SNE representations for the CA subclusters
# CA1 subclusters
plotmap(sc2_2cls, final = F)

# Plot the t-SNE heatmaps for the CA markers
# CA1
plotexpmap(sc2_2cls, "Pou3f1", logsc = T)


#############  Computing the DEG between the subclusters of interest  ############################

## DEGs subclusters 1 and 2 of CA1 ##

A <- names(sc2_2cls@cpart)[sc2_2cls@cpart == 1] 
B <- names(sc2_2cls@cpart)[sc2_2cls@cpart == 2]  # CLuster 2 is the one with higher expression of Pou3f1

DEGCA1_1vs2 <- diffexpnb(getfdata(sc2_2cls, n= c(A,B)), A=A, B=B)
DEGCA1_1vs2FTD <- DEGCA1_1vs2$res
DEGCA1_1vs2FTD$gene <- rownames(DEGCA1_1vs2FTD)
DEGCA1_1vs2FTD <- filter(DEGCA1_1vs2FTD, padj < 0.05)


####################  Volcano plot for the DEGs between the 2 CA1 sub-clusters of interest
data <- DEGCA1_1vs2$res
data <- data[rownames(data) != "Xist", ]
EnhancedVolcano(data, x= "log2FoldChange", y= "padj", lab = rownames(data), FCcutoff = 1.025339, pCutoff = 0.05, 
                transcriptPointSize = 2, transcriptLabSize = 3.5, transcriptLabhjust = -0.4, selectLab =  c("Tac2", "Pcp4", "Luzp2", "Pou3f1"), 
                ylab = bquote(~-Log[10]~adjusted~italic(P)), xlim = c(-4,4), title = NULL, subtitle = NULL, transcriptLabvjust = 0.7, boxedlabels = T, 
                drawConnectors = T, axisLabSize = 12, caption = NULL, legendLabSize = 12) 


### Over-representation test for genes increased in expression in CA1 subcluster 2 compared to subcluster 1

DEGCA1_1vs2FTD_Up <- filter(DEGCA1_1vs2FTD, foldChange >= 1)

AllgenesLog2FC <- as.numeric(DEGCA1_1vs2$res$log2FoldChange)
names(AllgenesLog2FC) <- rownames(DEGCA1_1vs2$res)

DEGCA1_1vs2_genes_Up <- as.numeric(DEGCA1_1vs2FTD_Up$log2FoldChange)
names(DEGCA1_1vs2_genes_Up) <- DEGCA1_1vs2FTD_Up$gene

egoCA1_1vs2Up <- enrichGO(gene = names(DEGCA1_1vs2_genes_Up), universe = names(AllgenesLog2FC), 
                          OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")


# Generate the barplot of the GO terms
data1 <- egoCA1_1vs2Up@result[1:8, ]
data1$order <- c(1:length(data1$ID))
ggplot(data1, aes(x = order, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Description, y = 0.1), hjust=0, vjust= 0.5, size= 4) +# nudge_x = c(0.22, -0.22)) + 
  coord_flip()  + 
  theme_classic() +
  scale_fill_gradient(low = "orange", high = "royalblue") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0,0)) 


###########################################################################################
################################### CA3 PYRAMIDAL CELLS ###################################
###########################################################################################

# Subset the data corresponding to the specific CA clusters
# CA3 cells
cells_cl4 <- names(sc@cpart[sc@cpart == 4])
prdata_cl4 <- prdata[ , which(names(prdata) %in% cells_cl4)] 


############# PERFORMING THE SUBCLUSTERING ANALYSIS DEFINING 2 CLUSTERS AND LOGPEARSON AS METRIC ###############

# Filter the data based on total counts
cs <- colSums(prdata_cl4)
prdata_cl4 <- prdata_cl4[,cs > 500]

# Initialize SC object with prdata
sc4_2cls <- SCseq(prdata_cl4)

# Filter cells with at least 'mintotal' number of  total transcripts
sc4_2cls <- filterdata(sc4_2cls, 
                       mintotal     =  3000, 
                       minexpr      =  5, 
                       minnumber    =  1,
                       LBatch       =  NULL, 
                       knn          =  10, 
                       CGenes       =  c("Mki67", "Pcna"), 
                       FGenes       =  c("Kcnq1ot1"),
                       ccor         =  0.4,
                       bmode        =  "RaceID")

# Compute cell-to-cell distance matrix for clustering and outlier identification
sc4_2cls <- compdist(sc4_2cls, 
                     metric        =  "logpearson", 
                     FSelect       =  T, 
                     knn           =  NULL)

# Perform clustering
sc4_2cls <- clustexp(sc4_2cls, 
                     sat          =  TRUE, 
                     #samp         =  1000, 
                     cln          =  2, 
                     clustnr      =  30,
                     bootnr       =  50, 
                     rseed        =  17000, 
                     FUNcluster   =  "kmedoids")

# Compute the t-SNE map
sc4_2cls <- comptsne(sc4_2cls, initial_cmd = TRUE, perplexity = 30, rseed = 15555)

# k-nearest neighbour graph layout utilizing the Fruchterman-Rheingold algorithm
sc4_2cls <- compfr(sc4_2cls,knn=10)

# Find outliers in clusters
sc4_2cls <- findoutliers(sc4_2cls, 
                         probthr      =  1e-5, 
                         outminc      =  3, 
                         outlg        =  3,
                         outdistquant = 0.95)



############## MAKING SOME REPRESENTATIVE PLOTS ########################################

# Plot the t-SNE representation for the markers expression in the original sc object
# CA3 marker
plotexpmap(sc, "Grik4", logsc = T)

# Plot the t-SNE representations for the CA subclusters
# CA3 subclusters
plotmap(sc4_2cls, final = F)

# Plot the t-SNE heatmaps for the CA markers
# CA3 subclusters
plotexpmap(sc4_2cls, "Grik4", logsc = T)


#############  Computing the DEG between the subclusters of interest  ############################

## DEGs subclusters 1 and 2 of CA3 ##

A <- names(sc4_2cls@cpart)[sc4_2cls@cpart == 1]
B <- names(sc4_2cls@cpart)[sc4_2cls@cpart == 2]

DEGCA3_1vs2 <- diffexpnb(getfdata(sc4_2cls, n= c(A,B)), A=A, B=B)
DEGCA3_1vs2FTD <- DEGCA3_1vs2$res
DEGCA3_1vs2FTD$gene <- rownames(DEGCA3_1vs2FTD)
DEGCA3_1vs2FTD <- filter(DEGCA3_1vs2FTD, padj < 0.05)


####################  Volcano plot for the DEGs between the 2 CA3 sub-clusters of interest

data <- DEGCA3_1vs2$res
data <- data[rownames(data) != "Xist", ]
EnhancedVolcano(data, x= "log2FoldChange", y= "padj", lab = rownames(data), FCcutoff = 1.025339, pCutoff = 0.05, 
                transcriptPointSize = 2, transcriptLabSize = 3.5, transcriptLabhjust = -0.4, selectLab =  c("Cpne4", "Grp", "Grik4"), 
                ylab = bquote(~-Log[10]~adjusted~italic(P)), xlim = c(-4,4), title = NULL, subtitle = NULL, transcriptLabvjust = 0.7, boxedlabels = T, 
                drawConnectors = T, axisLabSize = 12, caption = NULL, legendLabSize = 12) 


### Over-representation test for genes increased in expression in CA3 subcluster 2 compared to subcluster 1
DEGCA3_1vs2FTD_Up <- filter(DEGCA3_1vs2FTD, foldChange >= 1)


AllgenesLog2FC <- as.numeric(DEGCA3_1vs2$res$log2FoldChange)
names(AllgenesLog2FC) <- rownames(DEGCA3_1vs2$res)

DEGCA3_1vs2_genes_Up <- as.numeric(DEGCA3_1vs2FTD_Up$log2FoldChange)
names(DEGCA3_1vs2_genes_Up) <- DEGCA3_1vs2FTD_Up$gene

egoCA3_1vs2Up <- enrichGO(gene = names(DEGCA3_1vs2_genes_Up), universe = names(AllgenesLog2FC), 
                          OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

# Generate the barplot of the GO terms
data1 <- egoCA3_1vs2Up@result[1:8, ]
data1$order <- c(1:length(data1$ID))
ggplot(data1, aes(x = order, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Description, y = 0.1), hjust=0, vjust= 0.5, size= 4) +# nudge_x = c(0.22, -0.22)) + 
  coord_flip()  + 
  theme_classic() +
  scale_fill_gradient(low = "orange", high = "royalblue") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0,0)) 

