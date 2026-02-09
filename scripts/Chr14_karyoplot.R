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
Exp <- read.table(args[1], sep = "\t", header = T)
annotation <- read.table(args[2], header = F, sep = "\t")

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
SNP_data <- read.table(args[3])
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
# Draw karyoplot with expression and DR2 values.
# Note output file here will be the picture in PNG format
pdf(args[4], width = 15, height = 12)
kp <- plotKaryotype(genome = "hg38", plot.type = 2, chromosomes=c("chr14"))
kpDataBackground(kp, data.panel = 1, col="white")
kpDataBackground(kp, data.panel = 2, col="white")
kpAxis(kp, data.panel = 1, tick.pos = c(0, 5, 10), ymin=0, ymax=10, col = "black")
kpAddLabels(kp, labels = "Mean log2 TPM expression", srt=90, pos=1, label.margin = 0.045, ymax=15, ymin=0, col = "black")
kpRect(kp, chr="chr14", data.panel = 1, x0=merged_filtered$V4, x1=merged_filtered$V5, y0 =0, y1 = merged_filtered$Scale, r1 = 0.1, col="darkred",border="darkred")
kpPoints(kp, chr = "chr14", data.panel = 2, x = SNP_data_IMP$POS, y = SNP_data_IMP$DR2, col = "#92c5de", alpha = 0.4, r0 =1, r1 = -0)
kpPoints(kp, chr = "chr14", data.panel = 2, x = SNP_data_typed$POS, y = SNP_data_typed$DR2, col = "violet", r1 = 0)
kpAxis(kp, data.panel = 2, tick.pos = c(1, 0.8, 0.6, 0.3, 0.25), ymin=1, ymax=0, col = "black")
kpAddLabels(kp, data.panel = 2, labels = "Dosage R2", srt=90, pos=1, label.margin = 0.045, ymax=1, ymin=0, col = "black")


