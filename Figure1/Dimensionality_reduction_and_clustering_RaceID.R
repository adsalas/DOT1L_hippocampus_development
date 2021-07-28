
# Load required libraries

library(RaceID)

# Load the data
prdata <- read.table("~/Documents/Paper_Dot1l_Hippocampus/Data_submission_GEO/counts_matrix.txt", sep = "")
#prdata           <- as.data.frame(merged_data_tib)
#rownames(prdata) <- prdata$GENEID
#prdata$GENEID    <- NULL
# Remove ERCC and mitochondrial genes
#prdata           <- prdata[grep("^(ERCC|mt|Gm|Rik)",row.names(prdata),invert=TRUE),]

# Remove cells with less than or equal to 500 total transcripts
cs <- colSums(prdata)
prdata <- prdata[,cs > 500]

# Initialize SC object with prdata
sc <- SCseq(prdata)

# Filter cells with at least 'mintotal' number of  total transcripts
sc <- filterdata(sc, 
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
sc <- compdist(sc, 
               metric        =  "pearson", 
               FSelect       =  T, 
               knn           =  NULL)

# Perform clustering
sc <- clustexp(sc, 
               sat          =  TRUE, 
               #samp         =  1000, 
               cln          =  16, # Number of clusters is adjusted based on inspection of the saturation plot
               clustnr      =  30,
               bootnr       =  50, 
               rseed        =  17000, 
               FUNcluster   =  "kmedoids")

# Compute the t-SNE map
sc <- comptsne(sc, initial_cmd = TRUE, perplexity = 30, rseed = 15555)

# k-nearest neighbour graph layout utilizing the Fruchterman-Rheingold algorithm
sc <- compfr(sc,knn=10)

# Find outliers in clusters
sc <- findoutliers(sc, 
                   probthr      =  1e-5, 
                   outminc      =  3, 
                   outlg        =  3,
                   outdistquant = 0.95)

# Plot of change in log within cluster dispersion at different values of k. The results of this plot is used as a reference for defining the “cln” parameter in the clustexp() function.
plotsaturation(sc)


# Change the colors to get a little bit more of contrast between clusters 
sc@fcol <- sc@fcol[c(1,2,3,14,20,13,7,8,9,10,11,12,6,4,18,16,17,15,19,5,21,22,23,24,25,26,27)]

# tSNE representation of the data coloured by cluster
plotmap(sc)

# tSNE representation of the data coloured by sample
types <- substr(colnames(sc@ndata), start = 1, stop = 4)
plotsymbolsmap(sc, types, samples_col = c("#bababa","#ca0020","#f4a582"), cex = 0.5)
