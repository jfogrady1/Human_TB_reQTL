####################################
#
# Deconvolution Benchmark
####################################


## Install



# load
library("BisqueRNA")
library(MuSiC)
library(Biobase)
library(data.table)
library(Seurat)
library(tidyr)
library(SummarizedExperiment)
library("SimBu")


################################################################################################################
#### Goal to Generate this matrix for Pearson and RMSE #########################################################
##          Music_marker_raw   X      0      0
##          Music_raw_raw      X      0      0                                         # MUSIC requires raw dataframe at least for sc data
#           NNLS_marker_raw    X      0      0
## sc       NNLS_raw_raw       X      0      0
## data     Ciber_raw          X      X      X
##          Ciber_TPM          X      X      X
##          Bisque_raw         X      0      0                                         # Note internal normalisation of CPM so need raw
###########                   RAW....TPM....CPM
############################ Bulk data############################################################################





# Functions
TPM_sc = function(counts,lengths) {
  # counts = seurat object
  # lenghts = vector of gene lenghts with names corresponding to gene names 
  
  # Intersect for common - cannot TPM normalise genes with no length information
  A = intersect(rownames(counts),names(lengths))
  counts = counts[A,]
  lengths = lengths[A]
  
  print(length(A))
  # Perform TPM normalisation
  rate = counts@assays$RNA$counts / lengths
  apply(rate,2,function(x) 1e6*x/sum(x))
  counts@assays$RNA$TPM <- rate
  # return
  return(counts)
}

TPM_bulk = function(counts,lengths) {
  # counts = count matrix of bulk data
  # lenghts = vector of gene lenghts with names corresponding to gene names 
  
  # Intersect for common - cannot TPM normalise genes with no length information
  A = intersect(rownames(counts),names(lengths))
  counts = counts[A,]
  lengths = lengths[A]
  print(all(rownames(counts) == names(lengths)))
  # Perform TPM normalisation
  rate = counts/ lengths
  apply(rate,2,function(x) 1e6*x/sum(x))
  # return
  return(rate)
}



# Necessary files
# Need to convert gene ids to gene symbols
# Get in gene symbol information
symbols <- fread("/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf") %>% filter(V3 == "gene") # select(V9)
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)
symbols$length <- abs(symbols$V4 - symbols$V5)
symbols <- symbols %>% dplyr::select(gene_id, gene_name, length)


# Set up lengths data frame
lengths = symbols$length
names(lengths) = symbols$gene_name



# Read in the single cell data

# Read in marker genes)
markers <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Conserved_Marker_genes.txt")

TB.sc.data = readRDS("~/heQTL/work/scRNA_seq/TB.combined.final.rds")


# We will first normalise by TPM and make a signature matrix for cibersort
TB.sc.data.tpm = TPM_sc(counts = TB.sc.data, lengths = lengths)
TB.sc.data.tpm@assays$RNA$counts <- TB.sc.data.tpm@assays$RNA$TPM # overrule the counts when converting to SCE object

TB.sc.data.tpm@assays$RNA$counts
TB.sc.data$cell_type <- Idents(TB.sc.data)
# Average expression for all cell types
Average_raw_expression = AverageExpression(
  TB.sc.data,
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "cell_type",
  add.ident = NULL,
  layer = "counts",
  verbose = TRUE)

RAW_matrix <- as.matrix(Average_raw_expression$RNA)
RAW_matrix <- RAW_matrix[rownames(RAW_matrix) %in% markers$gene,]


TB.sc.data.tpm$cell_type <- Idents(TB.sc.data.tpm)
Average_TPM_expression = AverageExpression(
  TB.sc.data.tpm,
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "cell_type",
  add.ident = NULL,
  layer = "counts",
  verbose = TRUE)
TPM_matrix <- as.matrix(Average_TPM_expression$RNA)
TPM_matrix <- TPM_matrix[rownames(TPM_matrix) %in% markers$gene,]


write.table(TPM_matrix, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/TPM_signature_matrix_Cibersort.txt", sep = "\t", quote = F, row.names = T)
write.table(RAW_matrix, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/RAW_signature_matrix_Cibersort.txt", sep = "\t", quote = F, row.names = T)

# Remove for space
rm(TB.sc.data.tpm)


dim(markers)
TB.sc.data$ID <- colnames(TB.sc.data)

annotation = as.data.frame(TB.sc.data@meta.data)




# Create data structure for simulation
ds <- SimBu::dataset(
  annotation = annotation,
  count_matrix = GetAssayData(object = TB.sc.data, assay = "RNA", layer = "counts"), # need raw counts for the simulation
  name = annotation$cell_type
)

set.seed(4589)

# Perform simulation of bulk data
simulation <- SimBu::simulate_bulk(
  data = ds,
  scenario = "mirror_db",
  scaling_factor = "NONE",
  ncells = 10000,
  nsamples = 100,
  run_parallel = FALSE,
  balance_even_mirror_scenario = 0.5
) # multi-threading to TRUE
#> Using parallel generation of simulations.
#> Finished simulation.


# Structure manipulation
colnames(simulation$bulk) <- c(paste0("Sample_", 1:100))
rownames(simulation$cell_fractions) <- c(paste0("Sample_", 1:100))

# Look at plot
SimBu::plot_simulation(simulation = simulation) + theme_bw() + theme(axis.text.y  = element_text(size = 0, angle = 0))
ggsave("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Simbu_Simulation_Cell_Components.pdf", width = 12, height = 12, dpi = 600)


# 3 
#pseudo bulk data

#pseudo <- read.table('/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.pseudobulk.txt')

# convert to sce
as.matrix(assays(simulation$bulk)[["bulk_counts"]])
bulk_sim <- as.matrix(assays(simulation$bulk)[["bulk_counts"]])

# These are the raw counts
bulk_sim

# These are the TPM counts
bulk_TPM <- TPM_bulk(bulk_sim, lengths = lengths)


# These are the CPM counts
library(edgeR) 
bulk_CPM <- edgeR::cpm(bulk_sim)


# Write files
write.table(bulk_CPM, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Pseudobulk_CPM.txt", sep = "\t", quote = F)
write.table(bulk_TPM, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Pseudobulk_TPM.txt", sep = "\t", quote = F)
write.table(bulk_sim, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Pseudobulk_RAW.txt", sep = "\t", quote = F)



# Adjust our simulation to make a dataframe
simulation$TPM <- bulk_TPM



pheno_pseudo <- as.data.frame(matrix(nrow = 100, ncol = 1))
rownames(pheno_pseudo) <- colnames(simulation$TPM)
pheno_pseudo$V1 <- rownames(pheno_pseudo)
pseudo <- ExpressionSet(assayData = simulation$TPM, phenoData = AnnotatedDataFrame(pheno_pseudo))
pseudo@assayData$exprs
pseudo.tpm.mtx <- exprs(pseudo)

simulation$CPM <- bulk_CPM
pheno_pseudo <- as.data.frame(matrix(nrow = 100, ncol = 1))
rownames(pheno_pseudo) <- colnames(simulation$CPM)
pheno_pseudo$V1 <- rownames(pheno_pseudo)
pseudo <- ExpressionSet(assayData = simulation$CPM, phenoData = AnnotatedDataFrame(pheno_pseudo))
pseudo@assayData$exprs
pseudo.cpm.mtx <- exprs(pseudo)

pheno_pseudo <- as.data.frame(matrix(nrow = 100, ncol = 1))
rownames(pheno_pseudo) <- colnames(simulation$bulk)
pheno_pseudo$V1 <- rownames(pheno_pseudo)
pseudo <- ExpressionSet(assayData = as.matrix(assays(simulation$bulk)[["bulk_counts"]]), phenoData = AnnotatedDataFrame(pheno_pseudo))
pseudo@assayData$exprs
pseudo.raw.mtx <- exprs(pseudo) # last pseudo needed for bisque


pseudo.raw.mtx
pseudo.cpm.mtx
pseudo.tpm.mtx

# Read in marker genes)
markers <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Conserved_Marker_genes.txt")


# set up results

RESULTS <- simulation$cell_fractions
RESULTS$Mixture <- rownames(RESULTS)
RESULTS
RESULTS <- RESULTS %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "ground_truth")
RESULTS$Mixture <- factor(RESULTS$Mixture, levels = c(paste0("Sample_", 1:100)))
RESULTS <- RESULTS[order(RESULTS$cell_types, RESULTS$Mixture),]
head(RESULTS)

# RMSE - cibersort will be added later
# RMSE 

RMSE = matrix(nrow = 8, ncol = 3) # music, music + marker, nnls + marker, nnls, Bisque, Bisque + marker x TPM, Raw , CPM
PEARSON = matrix(nrow = 8, ncol = 3)
RMSE
rownames(RMSE) <- c("MuSIC_marker_raw", "NNLS_marker_raw", "music_raw", "NNLS_raw", "Bisque_marker", "Bisque_raw", "Cibersort_Raw", "Cibersort_TPM")
colnames(RMSE) <- c("RAW", "TPM", "CPM")

rownames(PEARSON) <- c("MuSIC_marker_raw", "NNLS_marker_raw", "music_raw", "NNLS_raw", "Bisque_marker", "Bisque_raw","Cibersort_Raw", "Cibersort_TPM")
colnames(PEARSON) <- c("RAW", "TPM", "CPM")


####################################################################

####################################################################

############## Benchmarking deconvolution #########################


####################################################################

####################################################################
####################################################################


###########################################################################
###########################################################################
###########################################################################
############## MUSIC + NNLS Benchmark with marker #########################
###########################################################################
###########################################################################
###########################################################################


# First try it with raw
TB.sc.data.sce <- as.SingleCellExperiment(TB.sc.data)
# Estimate cell type proportions
Est.prop.human.TB = music_prop(bulk.mtx = pseudo.raw.mtx, sc.sce  = TB.sc.data.sce, markers = markers$gene, clusters = 'ident', samples = 'sample_id', verbose = T)




music_results <- as.data.frame(Est.prop.human.TB$Est.prop.weighted)
music_results$Mixture <- rownames(music_results)
music_results_long <- music_results %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "observed")

music_results_long <- music_results_long %>% group_by(cell_types) %>% reframe(observed = observed,
                                                                                Mixture = Mixture)

nnls_results <- as.data.frame(Est.prop.human.TB$Est.prop.allgene)
nnls_results$Mixture <- rownames(nnls_results)
nnls_results_long <- nnls_results %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "observed")

nnls_results_long <- nnls_results_long %>% group_by(cell_types) %>% reframe(observed = observed,
                                                                              Mixture = Mixture)

head(RESULTS)

music_results_long <- music_results_long[order(as.character(music_results_long$cell_types)),]
nnls_results_long <- nnls_results_long[order(as.character(nnls_results_long$cell_types)),]

RESULTS
head(music_results_long)
RESULTS$Music_marker_observed <- as.vector(music_results_long$observed)
RESULTS$NNLS_marker_observed <- as.vector(nnls_results_long$observed)


RESULTS_music_marker_raw_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Music_marker_observed - ground_truth)^2)) %>% round(.,3), 
                             Pearson = cor(Music_marker_observed, ground_truth) %>% round(.,3),
                             MAD = mean(abs(Music_marker_observed - ground_truth)) %>% round(.,3))


RESULTS_music_marker_raw_df$Pearson <- if_else(is.na(RESULTS_music_marker_raw_df$Pearson),true = -1, false = RESULTS_music_marker_raw_df$Pearson)
# Mean for final matrix
RMSE[1,1] <- mean(RESULTS_music_marker_raw_df$RMSE)
PEARSON[1,1] <- mean(RESULTS_music_marker_raw_df$Pearson)




RESULTS_nnls_marker_raw_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((NNLS_marker_observed - ground_truth)^2)) %>% round(.,3), 
                                                        Pearson = cor(NNLS_marker_observed, ground_truth) %>% round(.,3),
                                                        MAD = mean(abs(NNLS_marker_observed - ground_truth)) %>% round(.,3))

RMSE[2,1] <- mean(RESULTS_nnls_marker_raw_df$RMSE)
PEARSON[2,1] <- mean(RESULTS_nnls_marker_raw_df$Pearson)


label = paste(" Pearson = ", round(PEARSON[1,1],3), "\n RMSE = ", round(RMSE[1,1],3), sep = "")

RESULTS$cell_types <- as.factor(RESULTS$cell_types)
P_Music_marker <- ggplot(RESULTS, aes(x = Music_marker_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label, size = 4)+
  ggtitle("MuSiC / markers / Raw counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "bottom")
P_Music_marker

label2 = paste(" Pearson = ", round(PEARSON[2,1],3), "\n RMSE = ", round(RMSE[2,1],3), sep = "")

P_NNLS_marker <-  ggplot(RESULTS, aes(x = NNLS_marker_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.55, label = label2, size = 4)+
  ggtitle("NNLS / markers / Raw counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "bottom")
P_NNLS_marker


####################################################################
####################################################################
####################################################################
############## MUSIC annotation with no marker #####################
####################################################################
####################################################################
####################################################################

# Estimate cell type proportions
Est.prop.human.TB = music_prop(bulk.mtx = pseudo.raw.mtx, sc.sce  = TB.sc.data.sce, markers = NULL, clusters = 'ident', samples = 'sample_id', verbose = T)


music_no_marker_results <- as.data.frame(Est.prop.human.TB$Est.prop.weighted)
music_no_marker_results$Mixture <- rownames(music_no_marker_results)
music_no_marker_results_long <- music_no_marker_results %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "observed")
music_no_marker_results_long <- music_no_marker_results_long %>% group_by(cell_types) %>% reframe(observed = observed,
                                                                              Mixture = Mixture)


nnls_no_marker_results <- as.data.frame(Est.prop.human.TB$Est.prop.allgene)
nnls_no_marker_results$Mixture <- rownames(nnls_no_marker_results)
nnls_no_marker_results_long <- nnls_no_marker_results %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "observed")

nnls_no_marker_results_long <- nnls_no_marker_results_long %>% group_by(cell_types) %>% reframe(observed = observed,
                                                                            Mixture = Mixture)


music_no_marker_results_long <- music_no_marker_results_long[order(as.character(music_no_marker_results_long$cell_types)),]
nnls_no_marker_results_long <- nnls_no_marker_results_long[order(as.character(nnls_no_marker_results_long$cell_types)),]


RESULTS$MUSIC_no_marker_observed <- as.vector(music_no_marker_results_long$observed)
RESULTS$NNLS_no_marker_observed <- as.vector(nnls_no_marker_results_long$observed)


RESULTS_music_raw_raw_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((MUSIC_no_marker_observed - ground_truth)^2)) %>% round(.,3), 
                                               Pearson = cor(MUSIC_no_marker_observed, ground_truth) %>% round(.,3),
                                               MAD = mean(abs(MUSIC_no_marker_observed - ground_truth)) %>% round(.,3))



RESULTS_music_raw_raw_df$Pearson <- if_else(is.na(RESULTS_music_raw_raw_df$Pearson),true = -1, false = RESULTS_music_raw_raw_df$Pearson)
RMSE[3,1] <- mean(RESULTS_music_raw_raw_df$RMSE)
PEARSON[3,1] <- mean(RESULTS_music_raw_raw_df$Pearson)

RESULTS_nnls_raw_raw_df = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((NNLS_no_marker_observed - ground_truth)^2)) %>% round(.,3), 
                                                Pearson = cor(NNLS_no_marker_observed, ground_truth) %>% round(.,3),
                                                MAD = mean(abs(NNLS_no_marker_observed - ground_truth)) %>% round(.,3))

RMSE[4,1] <- mean(RESULTS_nnls_raw_raw_df$RMSE)
PEARSON[4,1] <- mean(RESULTS_nnls_raw_raw_df$Pearson)


label3 =  paste(" Pearson = ", round(PEARSON[3,1],3), "\n RMSE = ", round(RMSE[3,1],3), sep = "")

label4 =  paste(" Pearson = ", round(PEARSON[4,1],3), "\n RMSE = ", round(RMSE[4,1],3), sep = "")


P_Music_no_marker <- ggplot(RESULTS, aes(x = MUSIC_no_marker_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label3, size = 4)+
  ggtitle("MuSiC / No markers / Raw counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")

P_NNLS_no_marker <- ggplot(RESULTS, aes(x = NNLS_no_marker_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label4, size = 4)+
  ggtitle("NNLS / no markers / Raw counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")


P_NNLS_no_marker
P_Music_no_marker

####################################################################
####################################################################
####################################################################
############## Bisque annotation with marker #######################
####################################################################
####################################################################
####################################################################


markers_bisque <- markers[,c(15,14,7)]
head(markers)
colnames(markers_bisque) <- c("gene", "cluster", "avg_logFC")

sc.eset <- BisqueRNA::SeuratToExpressionSet(TB.sc.data, delimiter="_", position=2, version="v3") # USE TPM again


res <- BisqueRNA::ReferenceBasedDecomposition(pseudo, sc.eset, markers = markers_bisque, verbose = T, use.overlap=FALSE) # Keep false as in real life we don't have matching so unfair comapirson

res$bulk.props
Bisque_results <- as.data.frame(t(res$bulk.props))
Bisque_results$Mixture <- rownames(Bisque_results)
Bisque_results_long <- Bisque_results %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")

Bisque_results_long <- Bisque_results_long[order(as.character(Bisque_results_long$cell_types)),]

RESULTS$Bisque_marker_observed <- Bisque_results_long$proportion

Bisque_marker_raw_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Bisque_marker_observed - ground_truth)^2)) %>% round(.,3), 
                                               Pearson = cor(Bisque_marker_observed, ground_truth) %>% round(.,3),
                                               MAD = mean(abs(Bisque_marker_observed - ground_truth)) %>% round(.,3))


RMSE[5,1] <- mean(Bisque_marker_raw_df$RMSE)
PEARSON[5,1] <- mean(Bisque_marker_raw_df$Pearson)

label5 = paste(" Pearson = ", round(PEARSON[5,1],3), "\n RMSE = ", round(RMSE[5,1],3))

P_bisque_marker <- ggplot(RESULTS, aes(x = Bisque_marker_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label5, size = 4) +
  ggtitle("Bisque / Markers / Raw counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")

P_bisque_marker

####################################################################
####################################################################
####################################################################
############## Bisque annotation with no marker ####################
####################################################################
####################################################################
####################################################################

res <- BisqueRNA::ReferenceBasedDecomposition(pseudo, sc.eset,verbose = T, use.overlap=FALSE, )

Bisque_results <- as.data.frame(t(res$bulk.props))
Bisque_results$Mixture <- rownames(Bisque_results)
Bisque_results_long <- Bisque_results %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")

Bisque_results_long <- Bisque_results_long[order(as.character(Bisque_results_long$cell_types)),]


RESULTS$Bisque_no_marker <- Bisque_results_long$proportion

Bisque_raw_raw_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Bisque_no_marker - ground_truth)^2)) %>% round(.,3), 
                                               Pearson = cor(Bisque_no_marker, ground_truth) %>% round(.,3),
                                               MAD = mean(abs(Bisque_no_marker - ground_truth)) %>% round(.,3))


RMSE[6,1] <- mean(Bisque_raw_raw_df$RMSE)
PEARSON[6,1] <- mean(Bisque_raw_raw_df$Pearson)
RMSE

label6 = paste(" Pearson = ", round(PEARSON[6,1],3), "\n RMSE = ", round(RMSE[6,1],3))
label6
P_bisque_no_marker <- ggplot(RESULTS, aes(x = Bisque_no_marker , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label6, size = 4)+
  ggtitle(" Bisque / No markers / Raw counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")


tail(Bisque_results_long,15)
tail(nnls_results_long,15)
tail(music_results_long,15)
tail(music_no_marker_results_long,15)
tail(nnls_no_marker_results_long, 15)

# Look at cibersort results
head(Ciber_RAW_RAW)
Ciber_RAW_RAW <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_RAW_RAW.txt") %>% as.data.frame() %>% dplyr::select(-c("RMSE", "Correlation", "P-value"))
Ciber_RAW_TPM <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_RAW_TPM.txt") %>% as.data.frame() %>% dplyr::select(-c("RMSE", "Correlation", "P-value"))
Ciber_RAW_CPM <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_RAW_CPM.txt") %>% as.data.frame() %>% dplyr::select(-c("RMSE", "Correlation", "P-value"))
Ciber_TPM_RAW <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_TPM_RAW.txt") %>% as.data.frame() %>% dplyr::select(-c("RMSE", "Correlation", "P-value"))
Ciber_TPM_TPM <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_TPM_TPM.txt") %>% as.data.frame() %>% dplyr::select(-c("RMSE", "Correlation", "P-value"))
Ciber_TPM_CPM <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_TPM_CPM.txt") %>% as.data.frame() %>% dplyr::select(-c("RMSE", "Correlation", "P-value"))

head(Ciber_RAW_RAW)
head(Ciber_RAW_CPM)
Ciber_RAW_RAW_long <- Ciber_RAW_RAW %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")
Ciber_RAW_RAW_long <- Ciber_RAW_RAW_long[order(as.character(Ciber_RAW_RAW_long$cell_types)),]

Ciber_RAW_TPM_long <- Ciber_RAW_TPM %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")
Ciber_RAW_TPM_long <- Ciber_RAW_TPM_long[order(as.character(Ciber_RAW_TPM_long$cell_types)),]

Ciber_RAW_CPM_long <- Ciber_RAW_CPM %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")
Ciber_RAW_CPM_long <- Ciber_RAW_CPM_long[order(as.character(Ciber_RAW_CPM_long$cell_types)),]

Ciber_TPM_RAW_long <- Ciber_TPM_RAW %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")
Ciber_TPM_RAW_long <- Ciber_TPM_RAW_long[order(as.character(Ciber_TPM_RAW_long$cell_types)),]

Ciber_TPM_TPM_long <- Ciber_TPM_TPM %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")
Ciber_TPM_TPM_long <- Ciber_TPM_TPM_long[order(as.character(Ciber_TPM_TPM_long$cell_types)),]

Ciber_TPM_CPM_long <- Ciber_TPM_CPM %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")
Ciber_TPM_CPM_long <- Ciber_TPM_CPM_long[order(as.character(Ciber_TPM_CPM_long$cell_types)),]


RESULTS$Ciber_RAW_RAW_observed <- Ciber_RAW_RAW_long$proportion
RESULTS$Ciber_RAW_TPM_observed <- Ciber_RAW_TPM_long$proportion
RESULTS$Ciber_RAW_CPM_observed <- Ciber_RAW_CPM_long$proportion
RESULTS$Ciber_TPM_RAW_observed <- Ciber_TPM_RAW_long$proportion
RESULTS$Ciber_TPM_TPM_observed <- Ciber_TPM_TPM_long$proportion
RESULTS$Ciber_TPM_CPM_observed <- Ciber_TPM_RAW_long$proportion

Ciber_RAW_RAW_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Ciber_RAW_RAW_observed - ground_truth)^2)) %>% round(.,7), 
                                                                           Pearson = cor(Ciber_RAW_RAW_observed, ground_truth) %>% round(.,3),
                                                                           MAD = mean(abs(Ciber_RAW_RAW_observed - ground_truth)) %>% round(.,3))

Ciber_RAW_TPM_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Ciber_RAW_TPM_observed - ground_truth)^2)) %>% round(.,3), 
                                                                          Pearson = cor(Ciber_RAW_TPM_observed, ground_truth) %>% round(.,3),
                                                                          MAD = mean(abs(Ciber_RAW_TPM_observed - ground_truth)) %>% round(.,3))


Ciber_RAW_TPM_df$Pearson <- if_else(is.na(Ciber_RAW_TPM_df$Pearson),true = -1, false = Ciber_RAW_TPM_df$Pearson)

Ciber_RAW_CPM_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Ciber_RAW_CPM_observed - ground_truth)^2)) %>% round(.,7), 
                                                                          Pearson = cor(Ciber_RAW_CPM_observed, ground_truth) %>% round(.,3),
                                                                          MAD = mean(abs(Ciber_RAW_CPM_observed - ground_truth)) %>% round(.,3))
Ciber_TPM_RAW_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Ciber_TPM_RAW_observed - ground_truth)^2)) %>% round(.,3), 
                                                                          Pearson = cor(Ciber_TPM_RAW_observed, ground_truth) %>% round(.,3),
                                                                          MAD = mean(abs(Ciber_TPM_RAW_observed - ground_truth)) %>% round(.,3))
Ciber_TPM_TPM_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Ciber_TPM_TPM_observed - ground_truth)^2)) %>% round(.,3), 
                                                                          Pearson = cor(Ciber_TPM_TPM_observed, ground_truth) %>% round(.,3),
                                                                          MAD = mean(abs(Ciber_TPM_TPM_observed - ground_truth)) %>% round(.,3))
Ciber_TPM_CPM_df <- RESULTS %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Ciber_TPM_CPM_observed - ground_truth)^2)) %>% round(.,3), 
                                                                          Pearson = cor(Ciber_TPM_CPM_observed, ground_truth) %>% round(.,3),
                                                                          MAD = mean(abs(Ciber_TPM_CPM_observed - ground_truth)) %>% round(.,3))

RMSE
RMSE <- RMSE[1:6,]
Ciber_RAW_CPM_df
Ciber_RAW_RAW_df
RMSE <- rbind(RMSE, c(mean(Ciber_RAW_RAW_df$RMSE), mean(Ciber_RAW_TPM_df$RMSE), mean(Ciber_RAW_CPM_df$RMSE)))
RMSE <- rbind(RMSE, c(mean(Ciber_TPM_RAW_df$RMSE), mean(Ciber_TPM_TPM_df$RMSE), mean(Ciber_TPM_CPM_df$RMSE)))
RMSE
Ciber_RAW_CPM_long
Ciber_RAW_RAW_long


rownames(RMSE)[7] <- "Cibersort_RAW"
rownames(RMSE)[8] <- "Cibersort_TPM"

RMSE_inverse = 1/RMSE
RMSE_inverse = round(1/RMSE)
          
PEARSON <- PEARSON[1:6,]

PEARSON <- rbind(PEARSON, c(mean(Ciber_RAW_RAW_df$Pearson), mean(Ciber_RAW_TPM_df$Pearson), mean(Ciber_RAW_CPM_df$Pearson)))
PEARSON <- rbind(PEARSON, c(mean(Ciber_TPM_RAW_df$Pearson), mean(Ciber_TPM_TPM_df$Pearson), mean(Ciber_TPM_CPM_df$Pearson)))
rownames(PEARSON)[7] <- "Cibersort_RAW"
rownames(PEARSON)[8] <- "Cibersort_TPM"

PEARSON
pdf(file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Benchmark_heatmap.pdf", width = 12, height = 12)
ComplexHeatmap::pheatmap(RMSE_inverse, display_numbers = T, color = colorRampPalette(c('white','red'))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15)
dev.off()
pdf(file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Benchmark_heatmap_pearson.pdf", width = 12, height = 12)
ComplexHeatmap::pheatmap(PEARSON, display_numbers = T, color = colorRampPalette(c('blue', 'white', "red"))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15)
dev.off()

P_bisque_no_marker
P_NNLS_marker
P_NNLS_no_marker
P_Music_marker
P_Music_no_marker
P_bisque_marker
P_bisque_no_marker

prow <- plot_grid(
  P_NNLS_marker + theme(legend.position="none"),
  P_NNLS_no_marker + theme(legend.position="none"),
  P_Music_marker + theme(legend.position="none"),
  P_Music_no_marker + theme(legend.position="none"),
  P_bisque_marker + theme(legend.position="none"),
  P_bisque_no_marker + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C", "D","E","F"),
  hjust = -1,
  nrow = 3
)


# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  P_Music_marker +
    guides(color = guide_legend(nrow = 3, override.aes = list(size = 5))) +
    labs(colour = "Cell Type") +
    theme(legend.position = "bottom")
)



# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
p_tpm <- plot_grid(prow, legend, ncol = 1, rel_heights = c(1, 0.1))
p_tpm

ggsave(plot = p_tpm, filename = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Pseudo_Benchmark_raw.png", width = 12, height = 12, dpi = 600)



# Now do the cibersort
RMSE
label7 = paste(" Pearson = ", round(PEARSON[7,1],3), "\n RMSE = ", round(RMSE[7,1],3))
label7
P_Ciber_RAW_RAW <- ggplot(RESULTS, aes(x = Ciber_RAW_RAW_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label7, size = 4)+
  ggtitle(" Cibersort / Raw / Raw counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")


RMSE
label8 = paste(" Pearson = ", round(PEARSON[7,2],3), "\n RMSE = ", round(RMSE[7,2],3))
label8
P_Ciber_RAW_TPM <- ggplot(RESULTS, aes(x = Ciber_RAW_TPM_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label8, size = 4)+
  ggtitle(" Cibersort / Raw / TPM counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")

label9= paste(" Pearson = ", round(PEARSON[7,3],3), "\n RMSE = ", round(RMSE[7,3],3))
label9
P_Ciber_RAW_CPM <- ggplot(RESULTS, aes(x = Ciber_RAW_CPM_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label9, size = 4)+
  ggtitle(" Cibersort / Raw / CPM counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")

label10= paste(" Pearson = ", round(PEARSON[8,1],3), "\n RMSE = ", round(RMSE[8,1],3))
label10
P_Ciber_TPM_RAW <- ggplot(RESULTS, aes(x = Ciber_TPM_RAW_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label10, size = 4)+
  ggtitle(" Cibersort / TPM / RAW counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")


label11= paste(" Pearson = ", round(PEARSON[8,2],3), "\n RMSE = ", round(RMSE[8,2],3))
label11
P_Ciber_TPM_TPM <- ggplot(RESULTS, aes(x = Ciber_TPM_TPM_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label11, size = 4)+
  ggtitle(" Cibersort / TPM / TPM counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")

label12= paste(" Pearson = ", round(PEARSON[8,3],3), "\n RMSE = ", round(RMSE[8,3],3))
label12
P_Ciber_TPM_CPM <- ggplot(RESULTS, aes(x = Ciber_TPM_CPM_observed , y = ground_truth, colour = cell_types)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
  annotate("text", x = 0.3, y = 0.75, label = label12, size = 4)+
  ggtitle(" Cibersort / TPM / CPM counts") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("observed fraction") +
  ylab("expected fraction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")


prow <- plot_grid(
  P_Ciber_RAW_RAW + theme(legend.position="none"),
  P_Ciber_TPM_RAW + theme(legend.position="none"),
  P_Ciber_RAW_TPM + theme(legend.position="none"),
  P_Ciber_TPM_TPM + theme(legend.position="none"),
  P_Ciber_RAW_CPM + theme(legend.position="none"),
  P_Ciber_TPM_CPM + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C", "D","E","F"),
  hjust = -1,
  nrow = 3
)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  P_Music_marker +
    guides(color = guide_legend(nrow = 3, override.aes = list(size = 5))) +
    labs(colour = "Cell Type") +
    theme(legend.position = "bottom")
)



# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
p_tpm <- plot_grid(prow, legend, ncol = 1, rel_heights = c(1, 0.1))

p_tpm
ggsave(plot = p_tpm, filename = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Deconvolution_Pseudo_Benchmark_Cibersort_raw.pdf", width = 12, height = 12, dpi = 600)



# Now check cibersorts performance on each cell type
colnames(RESULTS)
RESULTS_cibersort_RAW_CPM <- RESULTS %>% dplyr::select(1,2,3,12)

RESULTS_cibersort_RAW_CPM_values <- RESULTS_cibersort_RAW_CPM %>% group_by(cell_types) %>% dplyr::summarise(RMSE = sqrt(mean((Ciber_RAW_CPM_observed - ground_truth)^2)) %>% round(.,3), 
                                                                                                        Pearson = cor(Ciber_RAW_CPM_observed, ground_truth) %>% round(.,3),
                                                                                                        MAD = mean(abs(Ciber_RAW_CPM_observed - ground_truth)) %>% round(.,3))

RESULTS_cibersort_RAW_CPM %>% ggplot(., aes(x = Ciber_RAW_CPM_observed, y = ground_truth, col = cell_types)) + 
  geom_point() + 
  xlim(0,0.4) +
  ylim(0,0.4) +
  labs(x = "Observed fraction",
       y = "Expected fraction") +
  facet_wrap(~cell_types) +
  geom_text(data = RESULTS_cibersort_RAW_CPM_values, aes(x = 0.2, y = 0.37, label = paste0("RMSE = ",RMSE), color = NULL,group= NULL)) +
  geom_text(data = RESULTS_cibersort_RAW_CPM_values, aes(x = 0.2, y = 0.32, label = paste0(" Pearson = ",Pearson), color = NULL,group= NULL)) +
  theme(legend.position="none")
  
ggsave(filename = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Benchmark_Cibersort_all_celltypes_RAW_CPM.pdf", width = 12, height = 12, dpi = 600)


write.table(RESULTS, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/ALL_benchmark_results.text", sep = "\t", quote = FALSE)

