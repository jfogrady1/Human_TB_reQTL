# TPM normalisation

# John O'Grady

library("tidyverse")
args <- commandArgs(trailingOnly = TRUE)
# arg 1 = counts
counts0 <- as.matrix(read.table(args[1], sep = "\t", header = T, row.names = 1))

# arg 2 = annotation with gene length
annotation <- read.table(args[2], header = F, sep = "\t")


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


counts.table = cbind(counts0, annotation_simple)

x <- counts.table[,c(-length(colnames(counts.table)))] / counts.table[,length(colnames(counts.table))]
tpm.mat <- t( t(x) * 1e6 / colSums(x) )



tpm.mat <- as.data.frame(tpm.mat)

# Write to file
write.table(tpm.mat, args[3], sep = "\t", row.names = T, col.names = T, quote = F)
