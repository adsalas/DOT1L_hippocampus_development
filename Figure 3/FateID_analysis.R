#####################  FATEID ANALYSIS FOR HIPPOCAMPAL LINEAGES  ###################################

# Load the required packages
library(RaceID)
library(FateID)

# Filter out the cells that are not part of the hippocampal lineage
 
fdClust <- c( cellsC14 <- names(sc@cpart[sc@cpart == 14]), # interneurons
              cellsC15 <- names(sc@cpart[sc@cpart == 15]), # red blood cells
              cellsC19 <- names(sc@cpart[sc@cpart == 19]), # OPCs
              cellsC16 <- names(sc@cpart[sc@cpart == 16]), # Perycites
              cellsC17 <- names(sc@cpart[sc@cpart == 17]), # Unknown
              cellsC24 <- names(sc@cpart[sc@cpart == 24]), # Choroid plexus
              cellsC26 <- names(sc@cpart[sc@cpart == 26]), # Microglia
              cellsC27 <- names(sc@cpart[sc@cpart == 27]), # Unknown 1 cell
              cellsC22 <- names(sc@cpart[sc@cpart == 22]), # NSC5Unk_22 1 cell
              cellsC21 <- names(sc@cpart[sc@cpart == 21]), # Unkn2_21 10 cells
              cellsC7 <- names(sc@cpart[sc@cpart == 7]), # CR1
              cellsC25 <- names(sc@cpart[sc@cpart == 25]), # CR2
              cellsC10 <- names(sc@cpart[sc@cpart == 10]), # CR3
              cellsC20 <- names(sc@cpart[sc@cpart == 20])#, # CR4
              
)

# Filtering the matrix
mat <- sc@expdata[ , colnames(sc@ndata)]
mat <- mat[ , !colnames(sc@ndata) %in% fdClust]
dim(mat)

# Filter the tSNE coordinates
tSNEcoor <- sc@tsne[!colnames(sc@ndata) %in% fdClust, ]

# Filter the UMAP coordinates
UMAPcoor <- sc@umap[!colnames(sc@ndata) %in% fdClust, ]

# Filter the cpart
cpart <- sc@cpart[!names(sc@cpart) %in% fdClust]

# Provide the expression matrix as input (transcript counts)
x <- mat
x <- as.matrix(x)

# Filter the matrix based on selected features
x <- x[sc@cluster$features, ]

# Provide the partition of cells
y <- cpart

# Define the most mature populations in the data
tar <- c(1,2,4,6,13)

# Computing fate bias
fb  <- fateBias(x, 
                y, 
                tar, 
                z=NULL, 
                minnr=5, 
                minnrh=30, 
                adapt=TRUE, 
                nbfactor=5, 
                use.dist=FALSE, 
                seed=12345, 
                nbtree=NULL)

# Store the precomputed coordinates in the dimmensionality reduction step
dr <- list()
dr$tsne <- list("D2"= tSNEcoor)
dr$umap <- list("D2"= UMAPcoor)

# PLot the dimensionality reduced representation
plotFateMap(y,dr,k=2,m="tsne")

# Plot the fate bias towards specific clusters
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t1", n = "to DG granule")
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t2", n = "to CA1pyr")
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t4", n = "to CA3pyr")
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t6", n = "to Subic1")
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t13", n = "to Subic2")

