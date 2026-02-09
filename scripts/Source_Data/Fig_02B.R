# Script to plot mean expression along chromosome and also
# imputed variants with respective DR2 values.

# Packages
library(tidyverse)
library(ggplot2)
library(matrixStats)
library(karyoploteR)

args = commandArgs(trailingOnly=TRUE)

# Note
#args[1] = count_matrix
#args[2] = annotation file
#args[3] = SNP data for chromosome 14
#args[4] = Output PNG file.

# Expression and annotation file
# Note using all samples together as variants derived from this set
Exp <- read.table('work/RNA_seq/quantification/count_matrix_ordered_ALL.txt', sep = "\t", header = T)
annotation <- read.table('data/ref_genome/gencode.v43.annotation.gtf', header = F, sep = "\t")

# Modify annotation file
annotation <- annotation %>% filter(V3 == "gene")
annotation <- annotation %>% separate(., V9, into = c("gene_id", "gene_type", "gene_name", "level"), sep = ";") %>% select(1:9)
annotation$gene_id <- gsub("gene_id ", "", annotation$gene_id)
annotation$length <- abs(annotation$V4 - annotation$V5)
annotation_simple <- annotation %>% select(gene_id, length)

# Join together to get coordinates of genes for plotting
Exp <- left_join(Exp, annotation_simple, by = c("Geneid" = "gene_id"))


# TPM normalise to make it easier for viewing
x <- Exp[,c(-1,-length(colnames(Exp)))] / Exp$length
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
tpm.mat <- as.data.frame(tpm.mat)
tpm.mat$Mean <- rowMeans(as.matrix(tpm.mat))

# Rejoin the annotation
merged = cbind(tpm.mat, annotation)

# Select for chromosome 14
merged <- merged %>%filter(V1 == "chr14")

# Need to have expression > 0 for a variant to be detected
merged_filtered <- merged %>% filter(merged$Mean >0)
merged_filtered$Scale <- round(log2(merged_filtered$Mean + 1))


rm(annotation)
rm(annotation_simple)


# Now read in the SNP data
# This will come from the rule above where we extract data for chromosome 14
SNP_data <- read.table('work/DNA_seq/imputation/final/chr14_dr2.txt')
colnames(SNP_data) <- c("CHR", "POS", "META")
SNP_data <- SNP_data %>% separate("META", into = c("DR2", "AF", "IMP"), sep = ";")

# Modify to select IMPUTED and TYPED variants
SNP_data <- SNP_data %>% mutate(IMP = case_when(is.na(IMP) ~ "Typed",
                                                IMP == "IMP" ~ "Imputed")) %>% select(1,2,3,4,5)

# Modify to convert to numeric
SNP_data$DR2 = gsub("DR2=", "", SNP_data$DR2)
SNP_data$DR2 <- as.numeric(SNP_data$DR2)
SNP_data$AF = gsub("AF=", "", SNP_data$AF)
SNP_data$AF = as.numeric(SNP_data$AF)


SNP_data_typed <- SNP_data %>% filter(IMP == "Typed")
SNP_data_IMP <- SNP_data %>% filter(IMP == "Imputed")
SNP_data_IMP<- SNP_data_IMP %>% filter (DR2 >= 0.25)
print(dim(SNP_data_IMP))
print(dim(SNP_data_typed))
dim(SNP_data)

SNP_data <- SNP_data %>% filter(DR2 >= 0.25)
dim(SNP_data)
write.table(SNP_data, file = '/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_02_SNP_data.txt', row.names = FALSE, col.names = T, sep = "\t", quote = F)
write.table(merged_filtered, file = '/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_02_expression_data.txt', row.names = FALSE, col.names = T, sep = "\t", quote = F)
