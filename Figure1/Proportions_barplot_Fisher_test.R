
######  Plot a barplot with the cell proportions for each sample   ################################

cl <- sc@cpart 
df_cl <- data.frame(names= names(cl), cluster=cl)
df_cl$condition <- substr(df_cl$names, start = 0, stop = 4)
table(df_cl$condition, df_cl$cluster)
table_clusters <-table(df_cl$condition, df_cl$cluster)
table_clusters

#### Transforming the cell numbers to percentages for plotting a "corrected" barplot
sum_rows <- apply(table_clusters, 1, sum)
sum_rows <- as.table(sum_rows)
clus_cellsPerc <- apply(table_clusters, 2, function(x){x/sum_rows*100})

### 1. Filtering the clusters with just one or 2 cells
clus_cellsPercFTD <- clus_cellsPerc[ , -c(19,22,24,26,27)]
### 2. Horizontal Barplot
barplot(clus_cellsPerc, horiz = T, col = c("#bababa","#ca0020", "#f4a582"), border = "white", space = 0.3, legend= c("Control", "cKO1", "cKO2"), font.axis=2, xlab = "Percentage of cells")
# Filtered barplot
barplot(clus_cellsPerc, horiz = T, col = c("#bababa","#ca0020", "#f4a582"), border = "white", space = 0.3, legend= c("Control", "cKO1", "cKO2"), font.axis=2, xlab = "Percentage of cells")


#######   Perform the Fisher's exact test for comparing the changes in proportions per cluster between samples   #####################
# Subset the table
table_clusters1 <- table_clusters[1:2, ]
table_clusters2 <- table_clusters[c(1,3), ]
table_clusters3 <- table_clusters[2:3, ]


# Create a list for storing the information of the number of cells inside or outside the clusters
dataLIST <- list()
data <- data.frame("CellsIn"= c(1,1), "CellsOut"= c(1,1))

# Control vs cKO1

for (j in 1:dim(table_clusters1)[2]){
  for (i in 1:dim(table_clusters1)[1]) {
    data$CellsIn[i] <- table_clusters1[i,j]
    data$CellsOut[i] <- (sum(table_clusters1[i,]) - table_clusters1[i,j])
    dataLIST[[j]] <- data
    rownames(dataLIST[[j]]) <- rownames(table_clusters1)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST[[i]])
}


FISHER_RESULTS1 <- FISHER_RESULTS

# Control vs cKO2

for (j in 1:dim(table_clusters2)[2]){
  for (i in 1:dim(table_clusters2)[1]) {
    data$CellsIn[i] <- table_clusters2[i,j]
    data$CellsOut[i] <- (sum(table_clusters2[i,]) - table_clusters2[i,j])
    dataLIST[[j]] <- data
    rownames(dataLIST[[j]]) <- rownames(table_clusters2)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST[[i]])
}


FISHER_RESULTS2 <- FISHER_RESULTS

# cKO1 vs cKO2

for (j in 1:dim(table_clusters3)[2]){
  for (i in 1:dim(table_clusters3)[1]) {
    data$CellsIn[i] <- table_clusters3[i,j]
    data$CellsOut[i] <- (sum(table_clusters3[i,]) - table_clusters3[i,j])
    dataLIST[[j]] <- data
    rownames(dataLIST[[j]]) <- rownames(table_clusters3)
  }
}

FISHER_RESULTS <- list()

for (i in 1:length(dataLIST)) {
  FISHER_RESULTS[[i]] <- fisher.test(dataLIST[[i]])
}

FISHER_RESULTS3 <- FISHER_RESULTS


#######################################################################################################
###########   Correcting the p values with the p.adj() built-in function and Bonferroni  ##############
#######################################################################################################

# Control vs cKO1 (C817 vs K813)
pvalues1 <- list()
for (i in 1:length(FISHER_RESULTS1)) {
  pvalues1[i] <- FISHER_RESULTS1[[i]]$p.value
}
pvalues1 <- unlist(pvalues1)
# Adjusting the p value
pvalues1_adj <- p.adjust(pvalues1, method = "bonferroni")
# Checkign the padj that passed the threshold (<0.05)
pvalues1_adj < 0.05

# # Control vs cKO2 (C817 vs K814)
pvalues2 <- list()
for (i in 1:length(FISHER_RESULTS2)) {
  pvalues2[i] <- FISHER_RESULTS2[[i]]$p.value
}
pvalues2 <- unlist(pvalues2)
# Adjusting the p value
pvalues2_adj <- p.adjust(pvalues2, method = "bonferroni")
# Checkign the padj that passed the threshold (<0.05)
pvalues2_adj < 0.05

# cKO1 vs scKO2 (K813 vs K814)
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
pv_adj_bonf <- data.frame(cont_vs_cKO1 = pvalues1_adj < 0.05, 
                          cont_vs_cKO2 = pvalues2_adj < 0.05, 
                          cKO1_vs_cKO2 = pvalues3_adj < 0.05)
