# Script to CPM normalise the bulk rnaseq data
library(edgeR)
# Modify annotation file
# Need to extract gene information etc if want to convert to symbols
annotation <- read.table("/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf", header = F, sep = "\t")
annotation <- annotation %>% filter(V3 == "gene")
annotation <- annotation %>% separate(., V9, into = c("gene_id", "gene_type", "gene_name", "level"), sep = ";") %>% select(1:11)
annotation$gene_name <- gsub("gene_name ", "", annotation$gene_name)
annotation$gene_id <- gsub("gene_id ", "", annotation$gene_id)
annotation$length <- abs(annotation$V4 - annotation$V5)



annotation_simple <- annotation %>% select(gene_name, gene_id)

T0 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T0.txt") %>% 
  left_join(., annotation_simple, by = c("Geneid" = "gene_id")) %>% 
  mutate(duplicated = duplicated(gene_name)) %>% 
  mutate(gene_name = if_else(duplicated == TRUE, Geneid, gene_name)) %>%
  select(-c(Geneid, duplicated))
         
gene_ids <- T0$gene_name

T0 <- T0 %>% select(-c(gene_name)) %>% as.matrix()
rownames(T0) <- gene_ids


T1 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T1.txt") %>% select(-c(Geneid)) %>% as.matrix()
T2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T2.txt") %>% select(-c(Geneid)) %>% as.matrix()
T3 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T3.txt") %>% select(-c(Geneid)) %>% as.matrix()
T4 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T4.txt") %>% select(-c(Geneid)) %>% as.matrix()
rownames(T1) <- gene_ids
rownames(T2) <- gene_ids
rownames(T3) <- gene_ids
rownames(T4) <- gene_ids

colnames(T1) <- paste0(colnames(T1), "_T1")
colnames(T2) <- paste0(colnames(T2), "_T2")
colnames(T3) <- paste0(colnames(T3), "_T3")
colnames(T4) <- paste0(colnames(T4), "_T4")

bulk <- cbind(T0, T1, T2, T3, T4)
head(bulk)
bulk_cpm <- edgeR::cpm(bulk)
rownames(bulk_cpm) <- gsub(" ", "", rownames(bulk_cpm))

write.table(bulk_cpm, file = "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/Real_bulk_CPM_normalised.txt", sep = "\t", quote = F)
rownames(bulk_cpm)
