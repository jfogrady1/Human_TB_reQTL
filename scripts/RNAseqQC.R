library(tidyverse)
library(DESeq2)
library(ggplot2)
library(tidyquant)
library(reshape2)
library(patchwork)
# Read in the data for each timepoint

# Read in all of the files
counts0 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T0.txt", sep = "\t", header = T, row.names = 1))
counts1 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T1.txt", sep = "\t", header = T, row.names = 1))
counts2 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T2.txt", sep = "\t", header = T, row.names = 1))
counts3 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T3.txt", sep = "\t", header = T, row.names = 1))
counts4 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T4.txt", sep = "\t", header = T, row.names = 1))

# Read in covariate data
coldata_0 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T0/T0_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_1 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T1/T1_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_2 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T2/T2_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_3 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T3/T3_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_4 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T4/T4_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_3 <- coldata_3[,-3]
coldata_4 <- coldata_4[,-c(2,4,6)]

# make column for merged dataset
Time0 <- "T0"
Time1 <- "T1"
Time2 <- "T2"
Time3 <- "T3"
Time4 <- "T4"


coldata_0 <- cbind(coldata_0,Time0)
coldata_1 <- cbind(coldata_1,Time1)
coldata_2 <- cbind(coldata_2,Time2)
coldata_3 <- cbind(coldata_3,Time3)
coldata_4 <- cbind(coldata_4,Time4)


rownames(coldata_1) <- paste0(rownames(coldata_1),"_T1")
rownames(coldata_2) <- paste0(rownames(coldata_2),"_T2")
rownames(coldata_3) <- paste0(rownames(coldata_3),"_T3")
rownames(coldata_4) <- paste0(rownames(coldata_4),"_T4")
head(coldata_1)
colnames(counts1) <- paste0(colnames(counts1),"_T1")
colnames(counts2) <- paste0(colnames(counts2),"_T2")
colnames(counts3) <- paste0(colnames(counts3),"_T3")
colnames(counts4) <- paste0(colnames(counts4),"_T4")


edata <- cbind(counts0, counts1, counts2, counts3, counts4)


####
# Sex check
####


# arg 2 = annotation with gene length
annotation <- read.table("../data/ref_genome/gencode.v43.annotation.gtf", header = F, sep = "\t")


# Modify annotation file
# Need to extract gene length and gene information etc if want to convert to symbols
annotation <- annotation %>% filter(V3 == "gene")
annotation <- annotation %>% separate(., V9, into = c("gene_id", "gene_type", "gene_name", "level"), sep = ";") %>% select(1:11)
annotation$gene_name <- gsub("gene_name ", "", annotation$gene_name)
annotation$gene_id <- gsub("gene_id ", "", annotation$gene_id)
annotation$length <- abs(annotation$V4 - annotation$V5)



annotation_simple <- annotation %>% select(gene_name, gene_id, length)
rownames(annotation_simple) <- annotation_simple$gene_id
head(annotation_simple)
GeneSym <- annotation_simple$gene_name
annotation_simple <- annotation_simple %>% select(-c(1,2)) %>% as.matrix()

head(annotation_simple)
counts.table = cbind(edata, annotation_simple)
colnames(counts.table)
x <- counts.table[,c(-length(colnames(counts.table)))] / counts.table[,length(colnames(counts.table))]
tpm.mat <- t( t(x) * 1e6 / colSums(x) )

tpm.sex <- data.frame(t(tpm.mat[c("ENSG00000129824.16", "ENSG00000229807.13"),]))

tpm.sex$Sample <- rownames(tpm.sex)

tpm.sex$Sex <- coldata_4[,"gender"]

colnames(tpm.sex)[1] <- "RPS4Y1"
colnames(tpm.sex)[2] <- "XIST"

ggplot(data = tpm.sex, aes(x = XIST, y = RPS4Y1, col = Sex)) + geom_point() +
  theme_bw() +
  labs(x = "XIST (TPM normalised expression)",
       y = "RPS4Y1 (TPM normalised expression)",
       colour = "Reported Sex") +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
ggsave("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Sex_check.pdf", width = 12, height = 12, dpi = 600)


# Now onto boxplot

tpm.long <- data.frame(tpm.mat)
tpm.long$Gene_ID <- rownames(tpm.long)
dim(tpm.long)
tpm.long <- pivot_longer(tpm.long, cols = 1:240,names_to = "Sample_ID", values_to = "TPM_counts", )
tpm.long$Timepoint <- rep(c(rep("T0",48),rep("T1",48) ,rep("T2", 48),rep("T3", 48),rep("T4", 48)),
View(tpm.long)




















