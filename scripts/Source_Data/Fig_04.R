data <- readRDS("/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.final.rds")
library(Seurat)
library(ggplot2)
library(cowplot)
library(tidyverse)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.2, label.size = 4.5) + NoLegend()  +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10))) 


cell_embeddings <- data[["umap"]]@cell.embeddings

cell_annotations <- as.data.frame(data@active.ident)

all(rownames(cell_annotations) == rownames(cell_embeddings))

source_data <- cbind(cell_annotations, cell_embeddings)
head(source_data)
colnames(source_data)[1] <- "cell_type"

write.table(source_data, file = "/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_04A_source_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)



#### Deconvolution stuff

library(tidyverse)
library(ggplot2)
library(data.table)
library(speckle)
library(limma)
library(readr)
library(vctrs)
library(MuSiC)
library(Seurat)
library(SummarizedExperiment)
set.seed(4589)
library(Biobase)

args = commandArgs(trailingOnly=TRUE)
data <- fread("/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Results_Real_deconvolution.txt")
data <- data %>% select(1:16)

data_long <- data %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")

dim(data)
15*48
times <- c(rep("T0", 720 ), rep("T1", 720), rep ("T2", 720), rep("T3", 720), rep("T4", 720))
data_plot <- data_long
data_plot$Time <- times
my_palette = c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026")
data_plot %>% group_by(cell_types) %>% ggplot(., aes(x = Time, y = proportion, fill = Time)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.15, alpha = 0.3) +
  scale_fill_manual(values = my_palette) +
  facet_wrap(~ cell_types, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size= 12, colour = "black"),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 10)) +
  labs(y = "Cell type proportion",
       x = "Time")
head(data_plot)
write.table(data_plot, file = "/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_04B_source_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)
