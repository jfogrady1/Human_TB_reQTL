# Script to get eQTL numbers which are significant
# Also to get a violin plot of gene-wise FDR cutoffs


library(data.table)
library(tidyverse)
library(ggplot2)
library(UpSetR)
args = commandArgs(trailingOnly = T)
args[1] <- "/home/workspace/jogrady/heQTL/results/eQTL/T0.50000.cis_qtl_fdr0.1.txt"
args[2] <- "/home/workspace/jogrady/heQTL/results/eQTL/T1.50000.cis_qtl_fdr0.1.txt"  
args[3] <- "/home/workspace/jogrady/heQTL/results/eQTL/T2.50000.cis_qtl_fdr0.1.txt"
args[4] <- "/home/workspace/jogrady/heQTL/results/eQTL/T3.50000.cis_qtl_fdr0.1.txt"
args[5] <- "/home/workspace/jogrady/heQTL/results/eQTL/T4.50000.cis_qtl_fdr0.1.txt"
args[6] <- "/home/workspace/jogrady/heQTL/results/eQTL/T0.50000.cis_qtl_pairs.ALL.txt"  
args[7] <- "/home/workspace/jogrady/heQTL/results/eQTL/T1.50000.cis_qtl_pairs.ALL.txt"   
args[8] <- "/home/workspace/jogrady/heQTL/results/eQTL/T2.50000.cis_qtl_pairs.ALL.txt"   
args[9] <- "/home/workspace/jogrady/heQTL/results/eQTL/T3.50000.cis_qtl_pairs.ALL.txt"
args[10] <- "/home/workspace/jogrady/heQTL/results/eQTL/T4.50000.cis_qtl_pairs.ALL.txt"
  
  
data_T0 <- fread(args[1]) %>% filter(is_eGene == TRUE)
data_T1 <- fread(args[2]) %>% filter(is_eGene == TRUE)
data_T2 <- fread(args[3]) %>% filter(is_eGene == TRUE)
data_T3 <- fread(args[4]) %>% filter(is_eGene == TRUE)
data_T4 <- fread(args[5]) %>% filter(is_eGene == TRUE)


dim(data_T0)
dim(data_T1)
dim(data_T2)
dim(data_T3)
dim(data_T4)

T0_nominal <- fread(args[6])
T1_nominal <- fread(args[7])
T2_nominal <- fread(args[8])
T3_nominal <- fread(args[9])
T4_nominal <- fread(args[10])

# Get rid of the replicated headings
T0_nominal <- T0_nominal %>% filter(phenotype_id != "phenotype_id")
T1_nominal <- T1_nominal %>% filter(phenotype_id != "phenotype_id")
T2_nominal <- T2_nominal %>% filter(phenotype_id != "phenotype_id")
print(head(T2_nominal))
print(head(T3_nominal))
T3_nominal <- T3_nominal %>% filter(phenotype_id != "phenotype_id")
T4_nominal <- T4_nominal %>% filter(phenotype_id != "phenotype_id")


# Plot gene wise threshold
threshold_plot_T0 <- data_T0 %>% select(2,20)
threshold_plot_T1 <- data_T1 %>% select(2,20)
threshold_plot_T2 <- data_T2 %>% select(2,20)
threshold_plot_T3 <- data_T3 %>% select(2,20)
threshold_plot_T4 <- data_T4 %>% select(2,20)

threshold_plot_T0$Group <- "T0"
threshold_plot_T1$Group <- "T1"
threshold_plot_T2$Group <- "T2"
threshold_plot_T3$Group <- "T3"
threshold_plot_T4$Group <- "T4"

print(head(threshold_plot_T0))

my_palette = c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026")
threshold_plot <- rbind(threshold_plot_T0, threshold_plot_T1, threshold_plot_T2, threshold_plot_T3, threshold_plot_T4)



listInput <- list(T0 = as.vector(threshold_plot_T0$phenotype_id), T1 = as.vector(threshold_plot_T1$phenotype_id), T2 = as.vector(threshold_plot_T2$phenotype_id), T3 = as.vector(threshold_plot_T3$phenotype_id), T4 = as.vector(threshold_plot_T4$phenotype_id))
upset(fromList(listInput), order.by = "freq", sets.bar.color = c(my_palette),
sets.x.label = "cis-eGenes tested", point.size = 4, line.size = 2,
mainbar.y.label = "cis-eGenes tested intersections",
text.scale = 2.5, shade.alpha = 0.5)

write.table(fromList(listInput), file = "/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_03A_upset_input.txt", quote = F, row.names = F, sep = "\t")
