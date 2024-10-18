#####################################
#####################################

# Single Cell-RNA-seq preprocessing

#####################################
#####################################



#####################################
# 1. Libraries
#####################################
library(Matrix)
library(plyr)
library(dplyr)
library(Seurat)
library(sctransform)
library(igraph)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
require(Hmisc)
require(dplyr)
require(openxlsx)
require(ggplot2)
library(ggpubr)
require(cowplot)
library(data.table)
library(RColorBrewer)
library(SingleR)
library(scater)
library(pheatmap)
library(nichenetr)
library(patchwork)
library(DoubletFinder)
library(SoupX)
library(scDblFinder)
library(glmGamPoi)
library(scuttle)
library(SingleR)
library(SingleCellExperiment)
library(CellBench)
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(msigdbi)
library("cerebroApp")
library("presto")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
options(future.globals.maxSize = 1e10)

set.seed(49506)

#####################################
# 1. Functions
#####################################
args = commandArgs(trailingOnly=TRUE)
samples <- c(args[1], args[2], args[3], args[4], args[5])
#samples <- c("TB1", "TB2", "TB3", "LTBI1", "LTBI2")
samples <- gsub("/home/workspace/jogrady/heQTL/work/scRNA_seq/", "", samples)
pre_process_metrics <- function(sample_id) {
  path <- paste0("/home/workspace/jogrady/heQTL/work/scRNA_seq/",sample_id, "/outs/filtered_feature_bc_matrix")
  print(path)
  sobj <- Read10X(data.dir = path)
  sobj <- CreateSeuratObject(counts = sobj, min.cells = 0, min.features = 500, project = sample_id)
  sobj$sample_id <- sample_id


  # QC metrics <- very basic
  sobj$log1p_total_counts <- log1p(sobj@meta.data$nCount_RNA)
  sobj$log1p_n_genes_by_counts <- log1p(sobj@meta.data$nFeature_RNA)
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")
  sobj[["percent.globin"]] <- PercentageFeatureSet(sobj, pattern = "^HB[^(P)]")
  sobj[["percent.rrna"]] <- PercentageFeatureSet(sobj, pattern = "^RNA\\d8S5")
  
  if (sample_id == "LTBI1" | sample_id == "LTBI2") {
    sobj$Condition <- "LTBI"
    sobj$Sex <- "Female"
  }
  else if (sample_id == "TB1" | sample_id == "TB2"){
    sobj$Condition <- "TB"
    sobj$Sex <- "Female"
  }
  else {
    sobj$Condition <- "TB"
    sobj$Sex <- "Male"
  }
  
  if (sample_id == "LTBI1") {
    sobj$Age <- 51
  }
  else if (sample_id == "LTBI2") {
    
    sobj$Age <- 51
  }
  
  else if (sample_id == "TB1") {
    sobj$Age <- 33
    
  }
  else if (sample_id == "TB2") {
    
    sobj$Age <- 35
  }
  
  else {
    
    sobj$Age <- 34
  }
  return(sobj)
}
                 
##### Generating soup groups
get_soup_groups <- function(sobj) {
  
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:18, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  
  return(sobj@meta.data[['seurat_clusters']])
  
}

#### helper function for soup clusters
add_soup_groups <- function(sobj){
  sobj$soup_group <- get_soup_groups(sobj)
  return(sobj)
}


######### Actually make the soub group using the original counts
make_soup <- function(sobj){
  sample_id <- as.character(sobj$sample_id[1])
  print(sample_id)
  path <- paste0("/home/workspace/jogrady/heQTL/work/scRNA_seq/", sample_id, "/outs/raw_feature_bc_matrix/")
  raw <- Read10X(data.dir = path)
  sc = SoupChannel(raw,sobj@assays$RNA$counts) # need raw matrix and counts matrix from object
  sc = setClusters(sc, sobj$soup_group) # set the clusters based on soup group
  sc = autoEstCont(sc, doPlot=FALSE) # estimate contamination/soup channel
  out = adjustCounts(sc, roundToInt = TRUE) # get adjusted counts and round to nearest integer
  
  sobj[["original.counts"]] <- CreateAssayObject(counts = sobj@assays$RNA$counts)
  sobj@assays$RNA$counts <- out # redefine normalised counts
  
  return(sobj)
  
}





# Annotate the doublets
find_doublets <- function(sobj) {
  
  # sclDblFinder needs normalised data
  sobj2 <- NormalizeData(sobj, verbose = FALSE)
  sobj2 <- FindVariableFeatures(object = sobj2, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  sobj2 <- ScaleData(sobj2, verbose = FALSE)
  sobj2 <- RunPCA(sobj2, npcs = 20, verbose = FALSE)
  sobj2 <- FindNeighbors(sobj2, dims = 1:18, verbose = FALSE)
  sobj2 <- FindClusters(sobj2, resolution = 0.5, verbose = FALSE)
  
  # Make sure to use layer = data to use normalised counts
  # use clusters identified by seurat
  sce <- scDblFinder(GetAssayData(sobj2, layer = "counts"), dbr=0.01, clusters = Idents(sobj2))
  
  # return scores and classes
  sobj$scDblFinder.score <- sce$scDblFinder.score
  sobj$scDblFinder.class <- sce$scDblFinder.class
  return(sobj)
  
}



####### Actually filter the doublets
filter_doublets <- function(sobj) {
  
  sobj <- subset(sobj, cell = which(sobj$scDblFinder.class == "singlet"))
  
  return(sobj)
  
}


filter_seurat <- function(sobj) {
  
  # find outliers, subset and remove
  sobj <- subset(sobj, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 7)
  
  return(sobj)
}

final_preprocess <- function(sobj) {
  sobj2 <- NormalizeData(sobj, verbose = FALSE)
  sobj2 <- FindVariableFeatures(object = sobj2, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  return(sobj2)
}



# Preprocess and filter mitochondrial genes
data_list <- sapply(samples, pre_process_metrics)
data_list <- sapply(data_list, add_soup_groups)
# Ambient RNA removal
data_list <- sapply(data_list,make_soup)
head(data_list)

# View how much RNA was removed
sum(data_list[1]$TB1@assays$RNA$counts)/sum(data_list[1]$TB1@assays$original.counts@counts)
sum(data_list[2]$TB2@assays$RNA$counts)/sum(data_list[2]$TB2@assays$original.counts@counts)
sum(data_list[3]$TB3@assays$RNA$counts)/sum(data_list[3]$TB3@assays$original.counts@counts)
sum(data_list[4]$LTBI1@assays$RNA$counts)/sum(data_list[4]$LTBI1@assays$original.counts@counts)
sum(data_list[5]$LTBI2@assays$RNA$counts)/sum(data_list[5]$LTBI2@assays$original.counts@counts)


# Doublet detection
data_list <- sapply(data_list,find_doublets)


# filter doublets
data_list <- sapply(data_list, filter_doublets)


# Perform final filtering
data_list <- sapply(data_list, filter_seurat)


#####################################
# 3. Integration
#####################################

# Merge the data togehter
#TB.big <- merge(data_list[1]$LTBI1, y = c(data_list[2]$LTBI2,data_list[3]$TB1,data_list[4]$TB2,data_list[5]$TB3), add.cell.ids = c("LTBI1","LTBI2", "TB1", "TB2", "TB3"), project = "TB_singleCell")

# normalize and identify variable features for each dataset independently
TB.list <- lapply(X = data_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = TB.list)

TB.anchors <- FindIntegrationAnchors(object.list = TB.list, anchor.features = features, reduction = "cca")
# this command creates an 'integrated' data assay
TB.combined <- IntegrateData(anchorset = TB.anchors)



DefaultAssay(TB.combined) <- "integrated"


# Run standard pipeline
# Run the standard workflow for visualization and clustering
TB.combined <- ScaleData(TB.combined, verbose = F, vars.to.regress = c("nCount_RNA", "percent.mt", "nFeature_RNA", "percent.ribo"))
TB.combined <- RunPCA(TB.combined, npcs = 50, verbose = F)


#Select number of pcs
pct <- TB.combined[['pca']]@stdev / sum(TB.combined[['pca']]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_point() + 
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

pcs

DefaultAssay(TB.combined) <- "integrated"
TB.combined <- RunUMAP(TB.combined, reduction = "pca", dims = 1:pcs, verbose = F)
TB.combined <- FindNeighbors(TB.combined, reduction = "pca", dims = 1:pcs)
TB.combined <- FindClusters(TB.combined, resolution = 0.4)
p1 <- DimPlot(TB.combined, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(TB.combined, reduction = "umap", label = TRUE)
p2
dev.off()
p1
#ggsave("/home/workspace/jogrady/heQTL/work/scRNA_seq/UMAP_TB.no.batch.correction.pdf", width = 15, height = 15, dpi = 600)





#TB.big.harmony <- RunHarmony(TB.big, plot_convergence = T, c("sample_id"))

# Get coordinates of cells
#harmony_embeddings <- Embeddings(TB.big.harmony, 'harmony')
#harmony_embeddings[1:5, 1:5]

# VIsualise new cluster points
#p1 <- DimPlot(object = TB.big.harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
#p2 <- VlnPlot(object = TB.big.harmony, features = "harmony_1", group.by = "orig.ident", pt.size = .05) + NoLegend()
#plot_grid(p1,p2)

pcs

# Cluster using the harmony embeddings
#TB.big.harmony <- TB.big.harmony %>% 
 # RunUMAP(reduction = "harmony", verbose = F, dims = 1:pcs) %>% # Note these are the harmony dimensions not PCS
  #FindNeighbors(reduction = "harmony", dims = 1:pcs) %>%  # We find the neighbours on the harmony reduction but we plot the umap reduction
  #FindClusters(resolution = c(0.2,0.25,0.3,0.4,0.5,0.6))






# First rejoin layers
#TB.big.harmony[["RNA"]] <- JoinLayers(TB.big.harmony[["RNA"]])

# Conserved markers - for cell types that do not have enough cells (17, can use FindALL markers)

#**Determine the number of clusters **
num_clusters <- max(as.numeric(as.character(
  TB.combined@meta.data$integrated_snn_res.0.4)))
library(multtest)
library(metap)
num_clusters
levels(TB.combined)
DefaultAssay(TB.combined) <- "RNA"
TB.combined[["RNA"]] <- JoinLayers(TB.combined[["RNA"]])
all.markers_0 <- FindConservedMarkers(TB.combined, ident.1 = 0, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_1 <- FindConservedMarkers(TB.combined, ident.1 = 1, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_2 <- FindConservedMarkers(TB.combined, ident.1 = 2, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_3 <- FindConservedMarkers(TB.combined, ident.1 = 3, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_4 <- FindConservedMarkers(TB.combined, ident.1 = 4, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_5 <- FindConservedMarkers(TB.combined, ident.1 = 5, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_6 <- FindConservedMarkers(TB.combined, ident.1 = 6, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_7 <- FindConservedMarkers(TB.combined, ident.1 = 7, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_8 <- FindConservedMarkers(TB.combined, ident.1 = 8, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_9 <- FindConservedMarkers(TB.combined, ident.1 = 9, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_10 <- FindConservedMarkers(TB.combined, ident.1 =10, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_11<- FindConservedMarkers(TB.combined, ident.1 = 11, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_12 <- FindConservedMarkers(TB.combined, ident.1 = 12, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_13<- FindConservedMarkers(TB.combined, ident.1 = 13, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_14<- FindConservedMarkers(TB.combined, ident.1 = 14, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
all.markers_15<- FindConservedMarkers(TB.combined, ident.1 = 15, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)
#all.markers_16<- FindConservedMarkers(TB.combined, ident.1 = 16, grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE)


#saveRDS(TB.big.harmony, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.big.harmony.rds")



library(celldex)

####################################################################
####################################################################
####################################################################
############## Reference based annotation ##########################
####################################################################
####################################################################
####################################################################




# Load internal reference of human immune cells for SingleR
# https://bioconductor.org/packages/3.12/data/experiment/vignettes/celldex/inst/doc/userguide.html

references = list(
  HPCA =celldex::HumanPrimaryCellAtlasData(),
  BPEN = celldex::BlueprintEncodeData(),
  IMGEN =celldex::ImmGenData(),
  MONAC = celldex::MonacoImmuneData())


# Annotation is performed on cluster-level rather than default single-cell level
rst_HPCA<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$HPCA, labels = references$HPCA$label.main)

#Assign predicted cell labels to seurat object
TB.combined[["HPCA"]] <-rst_HPCA$labels

# Annotation is performed on cluster-level rather than default single-cell level
rst_BPEN<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$BPEN, labels = references$BPEN$label.main)

#Assign predicted cell labels to seurat object
TB.combined[["BPEN"]] <-rst_BPEN$labels

# Annotation is performed on cluster-level rather than default single-cell level
rst_IMGEN<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$IMGEN, labels = references$IMGEN$label.main)

#Assign predicted cell labels to seurat object
TB.combined[["IMGEN"]] <-rst_IMGEN$labels

# Annotation is performed on cluster-level rather than default single-cell level
rst_MONAC<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$MONAC, labels = references$MONAC$label.main)

#Assign predicted cell labels to seurat object
TB.combined[["MONAC"]] <-rst_MONAC$labels

#TB.combined@meta.data$R

singler_HPCA.results <- merge(data.frame(cell = rownames(rst_HPCA), singler = rst_HPCA$labels), 
                         data.frame(cell = rownames(TB.combined@meta.data), 
                                    cluster = TB.combined@meta.data$integrated_snn_res.0.4), 
                         by = "cell", 
                         all.y = FALSE)
singler_HPCA.results$cell <- NULL
singler_HPCA.results$count <- 1
singler_HPCA.results <- aggregate(count ~ ., singler_HPCA.results, FUN = sum)
singler_HPCA.final <- singler_HPCA.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)


singler_BPEN.results <- merge(data.frame(cell = rownames(rst_BPEN), singler = rst_BPEN$labels), 
                              data.frame(cell = rownames(TB.combined@meta.data), 
                                         cluster = TB.combined@meta.data$integrated_snn_res.0.4), 
                              by = "cell", 
                              all.y = FALSE)
singler_BPEN.results$cell <- NULL
singler_BPEN.results$count <- 1
singler_BPEN.results <- aggregate(count ~ ., singler_BPEN.results, FUN = sum)
singler_BPEN.final <- singler_BPEN.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)


singler_IMGEN.results <- merge(data.frame(cell = rownames(rst_IMGEN), singler = rst_IMGEN$labels), 
                              data.frame(cell = rownames(TB.combined@meta.data), 
                                         cluster = TB.combined@meta.data$integrated_snn_res.0.4), 
                              by = "cell", 
                              all.y = FALSE)
singler_IMGEN.results$cell <- NULL
singler_IMGEN.results$count <- 1
singler_IMGEN.results <- aggregate(count ~ ., singler_IMGEN.results, FUN = sum)
singler_IMGEN.final <- singler_IMGEN.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)

singler_IMGEN.final
singler_MONAC.results <- merge(data.frame(cell = rownames(rst_MONAC), singler = rst_MONAC$labels), 
                              data.frame(cell = rownames(TB.combined@meta.data), 
                                         cluster = TB.combined@meta.data$integrated_snn_res.0.4), 
                              by = "cell", 
                              all.y = FALSE)
singler_MONAC.results$cell <- NULL
singler_MONAC.results$count <- 1
singler_MONAC.results <- aggregate(count ~ ., singler_MONAC.results, FUN = sum)
singler_MONAC.final <- singler_MONAC.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)


singler_MONAC.final
singler_HPCA.final
singler_MONAC.final
dev.off()
#DimPlot(TB.combined, group.by = c("seurat_clusters", "HPCA", "BPEN", "MONAC", "IMGEN"), reduction = "umap", label = T)
DimPlot(TB.combined, group.by = c("integrated_snn_res.0.4", "BPEN", "HPCA"), reduction = "umap", label = T)

DimPlot(TB.combined, reduction = "umap", label = T) 

No_annotation <- DimPlot(TB.combined, group.by = c("integrated_snn_res.0.4"), reduction = "umap", label = T) 


# Finer annoation
# Check finer annotations for memory B cell
# Annotation is performed on cluster-level rather than default single-cell level
rst_HPCA.fine<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$HPCA, labels = references$HPCA$label.fine)

#Assign predicted cell labels to seurat object
TB.combined[["HPCA.fine"]] <-rst_HPCA.fine$pruned.labels

# Annotation is performed on cluster-level rather than default single-cell level
rst_BPEN.fine<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$BPEN, labels = references$BPEN$label.fine)

#Assign predicted cell labels to seurat object
TB.combined[["BPEN.fine"]] <-rst_BPEN.fine$pruned.labels

# Annotation is performed on cluster-level rather than default single-cell level
rst_IMGEN.fine<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$IMGEN, labels = references$IMGEN$label.fine)

#Assign predicted cell labels to seurat object
TB.combined[["IMGEN.fine"]] <-rst_IMGEN.fine$pruned.labels

# Annotation is performed on cluster-level rather than default single-cell level
rst_MONAC.fine<-SingleR(test = as.SingleCellExperiment(TB.combined), ref = references$MONAC, labels = references$MONAC$label.fine)

#Assign predicted cell labels to seurat object
TB.combined[["MONAC.fine"]] <-rst_MONAC.fine$pruned.labels
table(rst_HPCA.fine$pruned.labels)

BPEN <- DimPlot(TB.combined, reduction = "umap", group.by = c("BPEN.fine"), label = T, label.size = 7) + NoLegend()
MONAC <- DimPlot(TB.combined, reduction = "umap", group.by = c("MONAC.fine"), label = T, label.size = 4) 
HPCA <- DimPlot(TB.combined, reduction = "umap", group.by = c("HPCA"), label = T, label.size = 2)

MONAC
dev.off()
plot_grid(No_annotation, MONAC, BPEN)





####################################################################
####################################################################
####################################################################
############## Marker assisted annotation ##########################
####################################################################
####################################################################
####################################################################

# Here we will also use scmap for marker assisted annotation
library("SCINA")
# Import the marker genes as a GMT file and store as a variable
#
markers <- msigdbi::read.gmt(args[6])
# Convert the expression data from Seurat object into a matrix data structure
exprMatrix <- as.matrix(Seurat::GetAssayData(TB.combined))

# Run SCINA on the query data using the marker genes to identify cell types
# Specifying rm_overlap = FALSE allows the same marker gene to specify multiple cell types which
# may be useful if identifying cell subtypes or other similar types of cells
# Specifying allow_unknown = TRUE allows cells to be labeled as "unknown" instead of being
# assigned a low-confident label
predictions.scina = SCINA::SCINA(exp = exprMatrix, signatures = markers$genesets,
                                 rm_overlap = FALSE, allow_unknown = TRUE)
# Add SCINA annotation information to each cell in Seurat object
TB.combined[["SCINA"]] <- predictions.scina$cell_labels
DimPlot(TB.combined, reduction = "umap", group.by = c("SCINA"), label = T, label.size = 2)



# Get to here

###########################################
# At this stage, have rough annotations for cell types
# Can annotate clusters fully to improve this
# Use marker genes
##########################################
No_annotation <- DimPlot(TB.combined, reduction = "umap", label = T)
No_annotation
ggsave(args[7], width = 12, height = 12, dpi = 600)

save.image(file = args[8])

load(args[8])
#################################################
## Marker selection #############################
#################################################

# https://www.nature.com/articles/s41467-022-30893-5#MOESM7

######################################
# Platelets/ Megakaryocytes ##########
######################################
# Genes from https://www.science.org/doi/10.1126/sciadv.abh2169

No_annotation

# Can identify these clusters
DotPlot(TB.combined, features = c("GP9", "ITGA2B", "PF4", "PPBP")) # Cluster 11

#######################################
# Dendritic Cells
#######################################
# Look for marker genes in this cluster
# CLEC10A = marker for CDC1+ Dendritic cells

CLECL4 <- FeaturePlot(TB.combined, reduction = "umap", features = c("CLEC4C")) # Tiny little cluster there is plasmacytoid dendritic cells
CLECL4 # cluster 15
# Plasmacycoit dendritic cells = CLEC4C - on their own so will have to edit manually 


# Dendritic cells # Myeloid and plasmacytoid
# Genes from https://www.science.org/doi/10.1126/sciadv.abh2169


No_annotation
FeaturePlot(TB.combined, reduction = "umap", features = c("CD1C", "ITGAX", "CLEC4C", "CD14")) 
FeaturePlot(subset(TB.combined, seurat_clusters == "0"), features = c("CD1C"), min.cutoff = 0) # may need to rename as our own cluster

#######################################
# Monocytes
#######################################
# Genes from https://www.science.org/doi/10.1126/sciadv.abh2169

DotPlot(TB.combined, features = c("CD14", "LYZ", "S100A9", "CSF3R", "FCGR3A", "LYN", "CSF1R","IFITM1", "IFITM2", "IFITM3"))
# Classical monocytes = cluster = 0
# Non-classical monocyte = 12


#######################################
# CD4 + T cells
#######################################
# Genes from https://www.science.org/doi/10.1126/sciadv.abh2169
# Naive = Cluster 1 
# Memory = Cluster 2

No_annotation
DotPlot(TB.combined, features = c("CD3D", "CD3E", "CD3G", "CCR7", "SELL", "LTB", "IL7R", "S100A4", "CD8A", "CD8B", "GZMB", "PRF1"))
DotPlot(TB.combined, features = c("CD3D", "CD3E", "CD3G", "LTB", "IL7R", "S100A4"))
DotPlot(TB.combined, features = c("CD8A", "CD8B", "GZMB", "PRF1"))
#######################################
# CD8+ T cells
#######################################
# Genes from https://www.science.org/doi/10.1126/sciadv.abh2169
# 6 = Naieve CD8 T cell
# 5 + 8 = Memory CD8 T cell (a bit rough - likely has some NKT cells in there)



DotPlot(TB.combined, features = c("CD3D", "CD3E", "CD3G", "CCR7", "SELL", "CD8A", "CD8B", "GZMB", "PRF1"))
DotPlot(TB.combined, features = c("CD3D", "CD3E", "CD3G", "LTB","IL7R", "S100A4", "CD8A", "CD8B", "GZMB", "PRF1"))


#######################################
# NK cells
#######################################
# Genes from https://www.science.org/doi/10.1126/sciadv.abh2169
# Cluster 4+ 9 are NK cells



VlnPlot(subset(TB.combined, seurat_clusters %in% c("4", "9")), features = c("NKG7", "GNLY", "KLRC1","CD8A", "CD8B", "FCGR3A","GZMB", "PRF1"))
DotPlot(subset(TB.combined, seurat_clusters %in% c("4", "9")), features = c("NKG7", "GNLY", "KLRC1","CD8A", "CD8B", "FCGR3A","GZMB", "PRF1"))
FeaturePlot(subset(TB.combined, seurat_clusters %in% c("3", "10")), reduction = "umap", features = c("NKG7", "GNLY", "KLRC1","CD8A", "CD8B", "FCGR3A","GZMB", "PRF1"))




#######################################
# B Cells
#######################################
# Cluster 3 are Naieve B cells
# Cluster 7 are memory B cells
No_annotation
BPEN
MONAC # Bit of cluster 8 appears to be MAIT T cells
# Cluster 15 are plasma blasts
# 11 are platelets


# Cluster 10 are T reg cells

DotPlot(TB.combined, features = c("SLC4A10"))
DotPlot(TB.combined, features = c("CD79A", "MS4A1"))

# Get final IDS
No_annotation
New_ids = c("Classical Monocyte", #0 
            "Naïve CD4+ T cell", #1
            "Memory CD4+ T cell", #2
            "Naïve B cell", # 3
            "NK cell", #4
            "Memory CD8+ T cell", #5 
            "Naïve CD8+ T cell", #6
            "Memory B cell", #7
            "Memory CD8+ T cell", #8
            "NK cell", #9
            "T regulatory cell", #10, 
            "Megakaryocyte",#11
            "Non-classical Monocyte", # 12
            "Undefined T cell", # 13
            "Plasma B cell", # 14
            "Dendritic cell", #15,
            "MAIT cell" # 16 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8363247/ (SLC4A10 conserved marker)
           
)
New_ids
CellsMeta = TB.combined@meta.data
head(CellsMeta)
CellsMeta <- as.data.frame(CellsMeta)
CellsMeta <- CellsMeta %>% mutate(seurat_clusters = case_when(
                                                              MONAC.fine == "MAIT cells" ~ "16", 
                                                              .default = as.character(seurat_clusters)))
table(CellsMeta$seurat_clusters)

CellsMeta$seurat_clusters <- factor(CellsMeta$seurat_clusters, levels = as.character(c(0:16)))



TB.combined <- AddMetaData(TB.combined, CellsMeta)
table(TB.combined$seurat_clusters)


No_annotation_new <- DimPlot(TB.combined, group.by = c("seurat_clusters"), reduction = "umap", label = T) 
No_annotation_new



Idents(object = TB.combined) <- "seurat_clusters"

levels(TB.combined)

names(New_ids) <- levels(TB.combined)

#Idents(TB.combined)
TB.combined <- RenameIdents(TB.combined, New_ids)










DimPlot(TB.combined, reduction = "umap", label = TRUE, pt.size = 0.2, label.size = 4.5) + NoLegend() + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) 
ggsave(args[9], width = 10, height = 10, dpi = 600)


saveRDS(TB.combined, file = args[10])

TB.combined <- readRDS(args[10])


# Now onto generating the signature matrix
# find conserved markers for each cluster
levels(TB.combined)
all.markers_0_b <- FindConservedMarkers(TB.combined, ident.1 = "Memory CD8+ T cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_1_b <- FindConservedMarkers(TB.combined, ident.1 = "NK cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_2_b <- FindConservedMarkers(TB.combined, ident.1 = "Memory B cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_3_b <- FindConservedMarkers(TB.combined, ident.1 = "Memory CD4+ T cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_4_b <- FindConservedMarkers(TB.combined, ident.1 = "Naïve CD8+ T cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_5_b <- FindConservedMarkers(TB.combined, ident.1 = "Naïve B cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_6_b <- FindConservedMarkers(TB.combined, ident.1 = "Classical Monocyte", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_7_b <- FindConservedMarkers(TB.combined, ident.1 = "Naïve CD4+ T cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_8_b <- FindConservedMarkers(TB.combined, ident.1 = "T regulatory cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_9_b <- FindConservedMarkers(TB.combined, ident.1 = "MAIT cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_10_b <- FindConservedMarkers(TB.combined, ident.1 ="Undefined T cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.15, min.pct = 0.25,logfc.threshold = 0.5)
all.markers_11_b <- FindConservedMarkers(TB.combined, ident.1 = "Non-classical Monocyte", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_12_b <- FindConservedMarkers(TB.combined, ident.1 = "Megakaryocyte", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_13_b <- FindConservedMarkers(TB.combined, ident.1 = "Plasma B cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)
all.markers_14_b <- FindConservedMarkers(TB.combined, ident.1 = "Dendritic cell", grouping.var = "Condition", verbose = TRUE,  only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25,logfc.threshold = 0.75)

# add in cell types for denoting
all.markers_0_b$celltype = "Memory CD8+ T cell"
all.markers_1_b$celltype  = "NK cell"
all.markers_2_b$celltype  = "Memory B cell"
all.markers_3_b$celltype  = "Memory CD4+ T cell"
all.markers_4_b$celltype  = "Naïve CD8+ T cell"
all.markers_5_b$celltype  = "Naïve B cell"
all.markers_6_b$celltype  = "Classical Monocyte"
all.markers_7_b$celltype  = "Naïve CD4+ T cell"
all.markers_8_b$celltype  = "T regulatory cell"
all.markers_9_b$celltype  = "MAIT cell"
all.markers_10_b$celltype  ="Undefined T cell"
all.markers_11_b$celltype  = "Non-classical Monocyte"
all.markers_12_b$celltype  = "Megakaryocyte"
all.markers_13_b$celltype  = "Plasma B cell"
all.markers_14_b$celltype  = "Dendritic cell"

all.markers_0_b$gene <- rownames(all.markers_0_b)
all.markers_1_b$gene <- rownames(all.markers_1_b)
all.markers_2_b$gene <- rownames(all.markers_2_b)
all.markers_3_b$gene <- rownames(all.markers_3_b)
all.markers_4_b$gene <- rownames(all.markers_4_b)
all.markers_5_b$gene <- rownames(all.markers_5_b)
all.markers_6_b$gene <- rownames(all.markers_6_b)
all.markers_7_b$gene <- rownames(all.markers_7_b)
all.markers_8_b$gene <- rownames(all.markers_8_b)
all.markers_9_b$gene <- rownames(all.markers_9_b)
all.markers_10_b$gene <- rownames(all.markers_10_b)
all.markers_11_b$gene <- rownames(all.markers_11_b)
all.markers_12_b$gene <- rownames(all.markers_12_b)
all.markers_13_b$gene <- rownames(all.markers_13_b)
all.markers_14_b$gene <- rownames(all.markers_14_b)

all.conserved <- rbind(all.markers_0_b,
                       all.markers_1_b,
                       all.markers_2_b,
                       all.markers_3_b,
                       all.markers_4_b,
                       all.markers_5_b,
                       all.markers_6_b,
                       all.markers_7_b,
                       all.markers_8_b,
                       all.markers_9_b,
                       all.markers_10_b,
                       all.markers_11_b,
                       all.markers_12_b,
                       all.markers_13_b,
                       all.markers_14_b)




head(all.conserved)
dim(all.conserved)

all.conserved <- all.conserved %>% group_by(celltype) %>%
  arrange(as.numeric(minimump_p_val)) %>%
  slice_head(n = 100) %>%
  ungroup()



# Remove duplicates
all.conserved.unique <- all.conserved[!duplicated(all.conserved$gene),]

dim(all.conserved.unique)

write.table(all.conserved.unique, file=args[11], quote=FALSE, row.names = T, sep='\t', col.names = TRUE)

# Print the result
#TB.combined@active.ident


# write cell type to file for deconvolution
write.table(TB.combined@active.ident, file=args[12], quote=FALSE, sep='\t', col.names = TRUE)


# write gene counts per cell
write.table(TB.combined@assays$RNA$counts, file=args[13], quote=FALSE, sep='\t', col.names = TRUE)


# Perform pseudobulk
# Group by donor
TB.combined.pseudo <- AggregateExpression(TB.combined, assays = "RNA", group.by = "sample_id")

write.table(as.data.frame(TB.combined.pseudo$RNA), file = args[14], sep = "\t", quote = F, row.names = T, col.names = T)

write.table(as.data.frame(TB.combined@meta.data), file =args[15], sep = "\t", quote = F, row.names = T, col.names = T )
