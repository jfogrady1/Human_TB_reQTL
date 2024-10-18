library(data.table)
library(tidyverse)



args = commandArgs(trailingOnly=TRUE)
counts0 <- fread(args[1])
counts1 <- fread(args[2])
counts2 <- fread(args[3]) 
counts3 <- fread(args[4])
counts4 <- fread(args[5])


colnames(counts0)[2:49] <- paste0(colnames(counts0)[2:49], "_T0")
colnames(counts1)[2:49] <- paste0(colnames(counts1)[2:49], "_T1")
colnames(counts2)[2:49] <- paste0(colnames(counts2)[2:49], "_T2")
colnames(counts3)[2:49] <- paste0(colnames(counts3)[2:49], "_T3")
colnames(counts4)[2:49] <- paste0(colnames(counts4)[2:49], "_T4")

colnames(counts1)


all(counts0$Geneid == counts3$Geneid)

data_final <- cbind(counts0, counts1[,2:49], counts2[,2:49], counts3[,2:49], counts4[,2:49])


colnames(data_final)[1] <- "GeneSymbol" # Needed for cibersort

symbols <- fread(args[6]) %>% filter(V3 == "gene") %>% select(V9)
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols <- symbols %>% select(1,3)
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)

head(symbols)

data_final <- left_join(data_final, symbols, by = c("GeneSymbol" = "gene_id"))

data_final <- data_final %>% select(-1)
data_final <- data_final %>% select(c(241, 1:240))



colnames(data_final)[1] <- "GeneSymbol" # Needed for cibersort

write.table(data_final, file =args[7], sep = "\t", quote = FALSE, row.names = FALSE)
