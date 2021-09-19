
################## RNA velocity analysis with scvelo PART1: generate a "sce" object using Alevin/Salmon's outputs #############################

# Written by Michael Stadler and Charlotte Soneson, FMI, Basel, Switzerland (2020)

#################################################################################################################
####################### Read Alevin quantifications and generate SCE object  ####################################
#################################################################################################################

# Set root directory and load packages

## Set top directory. 
topdir <- "/path/to/your/files"   # adjust this to the path containing Alevin and Salmon outputs, and the metadata files

ncpu <- 32

## Load packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(umap)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(reshape2)
  library(tximeta)
  library(Matrix)
  library(parallel)
  library(org.Mm.eg.db)
  library(GenomicFeatures)
})


## helper function
read_alevin_cdna_introns <- function(alevindir, sampleid, allgenes) {
  cdna_introns <- tximeta(coldata = data.frame(
    names = sampleid,
    files = file.path(alevindir, "quants_mat.gz"),
    stringsAsFactors = FALSE
  ), type = "alevin")
  ## Get indexes of unspliced and spliced targets
  uidx <- grep("\\.*I\\.*$", rownames(cdna_introns))
  sidx <- grep("\\.*I\\.*$", rownames(cdna_introns), invert = TRUE)
  ## Subset count matrix to unspliced and spliced targets, respectively
  ucounts <- assay(cdna_introns, "counts")[uidx, ]
  scounts <- assay(cdna_introns, "counts")[sidx, ]
  rownames(ucounts) <- gsub("\\.*I\\.*$", "", rownames(ucounts))
  rownames(scounts) <- gsub("\\.*$", "", rownames(scounts))
  ssum <- sum(scounts)
  usum <- sum(ucounts)
  ## Match with allgenes
  ucounts <- ucounts[rownames(ucounts) %in% allgenes, ]
  scounts <- scounts[rownames(scounts) %in% allgenes, ]
  utmp <- stmp <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                               dims = c(length(allgenes), ncol(scounts)),
                               dimnames = list(allgenes, colnames(scounts)))
  utmp[rownames(ucounts), ] <- ucounts
  stmp[rownames(scounts), ] <- scounts
  ucounts <- utmp
  scounts <- stmp
  stopifnot(all(rownames(ucounts) == rownames(scounts)))
  stopifnot(all(colnames(ucounts) == colnames(scounts)))
  message(sampleid, ":")
  message("Spliced counts excluded since gene is missing from allgenes: ",
          round((ssum - sum(scounts)) / 1e6, 1), " Mio. (",
          round(100 * (ssum - sum(scounts)) / ssum, 1), "%)")
  message("Unspliced counts excluded since gene is missing from allgenes: ",
          round((usum - sum(ucounts)) / 1e6, 1), " Mio. (",
          round(100 * (usum - sum(ucounts)) / usum, 1), "%)\n")
  SingleCellExperiment(
    assays = list(spliced = scounts,
                  unspliced = ucounts)
  )
}

#################################################################################################################
############################### Read gene annotations and sample metadata  ######################################
#################################################################################################################

## Read gene annotations. Note: provide the respective path to the files
tx2gene <- read.delim(file.path(topdir, "annot/refsalmon_spliced_intron_separate_k23_fl50/tx2gene.tsv"),
                      as.is = TRUE, col.names = c("tx_id", "gene_id"), header = FALSE)

gnsymbols <- mapIds(x = org.Mm.eg.db,
                    keys = unique(tx2gene$gene_id),
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")

txdb <- loadDb(file.path(topdir, "RData/TxDb.sqlite"))

gnchrs <- mapIds(x = txdb,
                 keys = unique(tx2gene$gene_id),
                 column = "TXCHROM",
                 keytype = "GENEID",
                 multiVals = "first")

## Read metadata
meta <- read.csv(file.path(topdir, "samples.csv"), as.is = TRUE)

## Read known barcodes. Note: provide the respective path to the files
cellbc <- readLines(file.path(topdir, "metadata/celseq_barcodes.192.tabular"))


## Read cell clusters as defined in the metadata. Note: provide the respective path to the files
cellcluster <- readRDS(file.path(topdir, "metadata/cells_barcode_clusters.RDS"))
cellcluster <- cbind(cellcluster, cellid = paste0(cellcluster$barcode, "__", cellcluster$library))
head(cellcluster)



#################################################################################################################
############################  Read quantifications and add gene annotations  ####################################
#################################################################################################################

## define samples and genes
tmp <- read.csv(file.path(topdir, "samples.csv"), as.is = TRUE)
#sids <- c("C81721", "C81722")
sids <- tmp$SampleName
sids

names(sids) <- sids
allgenes <- unique(tx2gene$gene_id)
length(allgenes)

## read counts
setTximetaBFC("/path/to/your/local/BiocFileCache/if/different/from/default") # optional (only needed if you want to set a non-default location)

tmpL <- mclapply(sids, function(sid) {
  read_alevin_cdna_introns(alevindir = file.path(topdir, "data_processed",
                                                 "alevin_spliced_introns_separate", sid, "alevin"),
                           sampleid = sid,
                           allgenes = allgenes)
}, mc.cores = ncpu)

## keep expected cells only
sapply(tmpL, ncol)

sapply(tmpL, function(x) sum(colnames(x) %in% cellbc))

tmpL <- mclapply(tmpL, function(x) x[, colnames(x) %in% cellbc], mc.cores = ncpu)
sapply(tmpL, ncol)

## add experiment to cell barcode to make them unique
tmpL <- mclapply(names(tmpL), function(nm) {
  x <- tmpL[[nm]]
  colnames(x) <- paste0(colnames(x), "__", nm)
  x
}, mc.cores = ncpu)
sce <- do.call(cbind, tmpL)

## map cells to clusters as defined in the metadata
summary(colnames(sce) %in% cellcluster$cellid)

ncol(sce)

nrow(cellcluster)

clustername <- ifelse(colnames(sce) %in% cellcluster$cellid,
                      as.character(cellcluster[match(colnames(sce), cellcluster$cellid), "cluster"]),
                      "unknown")
table(clustername)

## add row/columns data
cd <- DataFrame(barcode = sapply(strsplit(colnames(sce), "__"), "[", 1),
                experiment = sapply(strsplit(colnames(sce), "__"), "[", 2),
                cluster = clustername)
rd <- DataFrame(gene_id = rownames(sce),
                symbol = gnsymbols[rownames(sce)],
                chr = gnchrs[rownames(sce)])
colData(sce) <- cd
rowData(sce) <- rd


#################################################################################################################
#################################### Calculate QC metrics with scater  ##########################################
#################################################################################################################

(Mt <- rownames(sce)[rowData(sce)$chr == "MT"])

sce <- scater::addPerCellQC(sce,
                            subsets = list(Mt = Mt),
                            BPPARAM = MulticoreParam(48),
                            exprs_values = "spliced")
sce <- scater::addPerFeatureQC(sce,
                               BPPARAM = MulticoreParam(48),
                               exprs_values = "spliced")

#################################################################################################################
######################################  Initial QC and filtering   ##############################################
#################################################################################################################

detected_features_threshold <- 1000
Mt_percent_threshold <- 25

sce$retain <- sce$detected > detected_features_threshold &
  sce$subsets_Mt_percent < Mt_percent_threshold
table(sce$experiment, sce$retain)

ggplot(data.frame(colData(sce)), aes(x = experiment, fill = retain)) +
  geom_bar() + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("") + ylab("Number of cells")


ggplot(data.frame(colData(sce)),
       aes(x = sum, fill = retain)) +
  geom_histogram(bins = 40, alpha = 0.5) + theme_cowplot() +
  xlab("Total UMI count per cell") + ylab("Number of cells")


ggplot(data.frame(colData(sce)),
       aes(x = detected, fill = retain)) +
  geom_histogram(bins = 40, alpha = 0.5) + theme_cowplot() +
  geom_vline(xintercept = detected_features_threshold, linetype = 3) +
  xlab("Number of detected genes per cell") + ylab("Number of cells")


ggplot(data.frame(colData(sce)),
       aes(x = sum, y = detected,
           color = retain)) +
  geom_point(size = 2, alpha = 0.4) + theme_cowplot() +
  geom_hline(yintercept = detected_features_threshold, linetype = 3) +
  xlab("Total UMI count per cell") + ylab("Number of detected genes per cell") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

ggplot(data.frame(colData(sce)) %>% dplyr::arrange(desc(detected)),
       aes(x = seq_along(detected), y = detected, color = retain)) +
  geom_line(size = 2) + theme_cowplot() +
  geom_hline(yintercept = detected_features_threshold, linetype = 3) +
  xlab("Rank") + ylab("Number of detected genes per cell")


## Mitochondrial content

ggplot(data.frame(colData(sce)), aes(x = detected, y = subsets_Mt_percent, colour = retain)) +
  geom_point(size = 2, alpha = 0.4) + theme_cowplot() +
  geom_hline(yintercept = Mt_percent_threshold, linetype = 3) +
  geom_vline(xintercept = detected_features_threshold, linetype = 3) +
  xlab("Number of detected genes per cell") + ylab("Percentage of mitochondrial counts") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

## Filter
sce <- sce[, sce$retain]
sce

table(sce$cluster)

## Most highly expressed genes
scater::plotHighestExprs(sce[!is.na(rowData(sce)$symbol), ],
                         n = 25, exprs_values = "spliced",
                         feature_names_to_plot = "symbol")


#################################################################################################################
################################  Calculate normalization factors  ##############################################
#################################################################################################################

size_factor_threshold <- 0.1

# first quickly find rough clusters
# ... using scran::quickCluster
blockfact <- if (ncol(sce) > 1000) cut(seq_len(ncol(sce)), round(ncol(sce) / 500)) else NULL
clusters <- scran::quickCluster(sce, method = "igraph", min.size = 100, use.ranks = FALSE,
                                BPPARAM = MulticoreParam(workers = round(60 / nlevels(blockfact))),
                                assay.type = "spliced",
                                block = blockfact, block.BPPARAM = MulticoreParam(workers = nlevels(blockfact)))
print(table(clusters))

# ... make sure there are not singletons
ncl <- table(clusters)
if (any(f <- ncl == 1)) { # throw singletons into largest clusters
  clusters[clusters %in% which(f)] <- which.max(ncl)
}
print(table(clusters))

# cacluate norm factors using deconvolution within cluster
sce <- scran::computeSumFactors(sce, min.mean = 0.5, cluster = clusters,
                                assay.type = "spliced",
                                BPPARAM = MulticoreParam(workers = 20))
print(summary(sizeFactors(sce)))

print(ggplot(data.frame(libSize = sce$sum,
                        sizeFactor = sizeFactors(sce),
                        stringsAsFactors = FALSE),
             aes(x = libSize, y = sizeFactor)) +
        geom_point(color = "lightgrey") + theme_cowplot() +
        geom_hline(yintercept = size_factor_threshold, linetype = "dashed") +
        xlab("Total UMI count") + ylab("Size factor") + ggtitle(metadata(sce)$dataset))

summary(sizeFactors(sce) > size_factor_threshold)

sce <- sce[, sizeFactors(sce) > size_factor_threshold]

#################################################################################################################
############################  Normalize (calculate logcounts)  ##################################################
#################################################################################################################

sce <- scater::logNormCounts(sce, exprs_values = "spliced")


#################################################################################################################
#################################### Dimension reduction  #######################################################
#################################################################################################################

set.seed(42)
sce <- scater::runPCA(sce, exprs_values = "logcounts", ncomponents = 30)
set.seed(43)
sce <- scater::runTSNE(sce, dimred = "PCA")
set.seed(44)
sce <- scater::runUMAP(sce, dimred = "PCA")

plotReducedDim(sce, dimred = "PCA",  colour_by = "experiment") + ggtitle("PCA")

plotReducedDim(sce, dimred = "TSNE", colour_by = "experiment") + ggtitle("T-SNE")

plotReducedDim(sce, dimred = "UMAP", colour_by = "experiment") + ggtitle("UMAP")

plotReducedDim(sce, dimred = "PCA",  colour_by = "cluster") + ggtitle("PCA")

plotReducedDim(sce, dimred = "TSNE", colour_by = "cluster") + ggtitle("T-SNE")

plotReducedDim(sce, dimred = "UMAP", colour_by = "cluster") + ggtitle("UMAP")


#################################################################################################################
################################## Save SingleCellExperiment objects  ###########################################
#################################################################################################################

# make "spliced" also the "counts" assay, as anndata2r will make the first assay X and thus loose "spliced"
assays(sce) <- c(SimpleList(counts = assay(sce, "spliced")), assays(sce))
saveRDS(sce, file = "sce.rds")


#################################################################################################################
#########################################    Session info     ###################################################
#################################################################################################################

## R version 3.6.1 (2019-07-05)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /tungstenfs/groups/gbioinfo/Appz/easybuild/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_skylakex-r0.3.7.so
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] GenomicFeatures_1.38.2      org.Mm.eg.db_3.10.0         AnnotationDbi_1.48.0        Matrix_1.2-18              
##  [5] tximeta_1.4.3               reshape2_1.4.3              cowplot_1.0.0               tidyr_1.0.2                
##  [9] dplyr_0.8.4                 tibble_2.1.3                umap_0.2.4.1                scran_1.14.6               
## [13] scater_1.14.6               ggplot2_3.2.1               SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
## [17] DelayedArray_0.12.2         BiocParallel_1.20.1         matrixStats_0.55.0          Biobase_2.46.0             
## [21] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0         IRanges_2.20.2              S4Vectors_0.24.3           
## [25] BiocGenerics_0.32.0         RColorBrewer_1.1-2         
## 
## loaded via a namespace (and not attached):
##  [1] Rtsne_0.15               ggbeeswarm_0.6.0         colorspace_1.4-1         XVector_0.26.0           BiocNeighbors_1.4.1     
##  [6] farver_2.0.3             bit64_0.9-7              RSpectra_0.16-0          tximport_1.14.0          knitr_1.28              
## [11] jsonlite_1.6.1           Rsamtools_2.2.3          dbplyr_1.4.2             uwot_0.1.5               compiler_3.6.1          
## [16] httr_1.4.1               dqrng_0.2.1              assertthat_0.2.1         lazyeval_0.2.2           limma_3.42.2            
## [21] BiocSingular_1.2.2       htmltools_0.4.0          prettyunits_1.1.1        tools_3.6.1              rsvd_1.0.3              
## [26] igraph_1.2.4.2           gtable_0.3.0             glue_1.3.1               GenomeInfoDbData_1.2.2   rappdirs_0.3.1          
## [31] Rcpp_1.0.3               vctrs_0.2.3              Biostrings_2.54.0        rtracklayer_1.46.0       DelayedMatrixStats_1.8.0
## [36] xfun_0.12                stringr_1.4.0            lifecycle_0.1.0          irlba_2.3.3              ensembldb_2.10.2        
## [41] statmod_1.4.34           XML_3.99-0.3             edgeR_3.28.0             zlibbioc_1.32.0          scales_1.1.0            
## [46] hms_0.5.3                ProtGenerics_1.18.0      AnnotationFilter_1.10.0  yaml_2.2.1               curl_4.3                
## [51] memoise_1.1.0            reticulate_1.14          gridExtra_2.3            biomaRt_2.42.0           stringi_1.4.6           
## [56] RSQLite_2.2.0            rlang_0.4.4              pkgconfig_2.0.3          bitops_1.0-6             evaluate_0.14           
## [61] lattice_0.20-38          purrr_0.3.3              GenomicAlignments_1.22.1 labeling_0.3             bit_1.1-15.2            
## [66] tidyselect_1.0.0         plyr_1.8.5               magrittr_1.5             R6_2.4.1                 DBI_1.1.0               
## [71] pillar_1.4.3             withr_2.1.2              RCurl_1.98-1.1           crayon_1.3.4             BiocFileCache_1.10.2    
## [76] rmarkdown_2.1            viridis_0.5.1            progress_1.2.2           locfit_1.5-9.1           grid_3.6.1              
## [81] FNN_1.1.3                blob_1.2.1               digest_0.6.25            RcppParallel_4.4.4       openssl_1.4.1           
## [86] munsell_0.5.0            beeswarm_0.2.3           viridisLite_0.3.0        vipor_0.4.5              askpass_1.1