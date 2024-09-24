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
ggplot(data = threshold_plot, aes(x = Group, y = -log10(pval_nominal_threshold), fill = Group)) + 
geom_violin(alpha = 0.5) + 
geom_boxplot(width = 0.07) + 
scale_fill_manual(values = my_palette) + 
theme_bw() +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15))
ggsave(args[11], width = 10, height = 12, dpi = 600)

# get pvalues for significance
all_control_pvalue <- wilcox.test(threshold_plot_T0$pval_nominal_threshold, threshold_plot_T1$pval_nominal_threshold)$p.value
control_reactor_pvalue <- wilcox.test(threshold_plot_T1$pval_nominal_threshold, threshold_plot_T2$pval_nominal_threshold)$p.value
all_reactor_pvalue <- wilcox.test(threshold_plot_T0$pval_nominal_threshold, threshold_plot_T2$pval_nominal_threshold)$p.value

# Plot upset plot of genes tested

listInput <- list(T0 = as.vector(threshold_plot_T0$phenotype_id), T1 = as.vector(threshold_plot_T1$phenotype_id), T2 = as.vector(threshold_plot_T2$phenotype_id), T3 = as.vector(threshold_plot_T3$phenotype_id), T4 = as.vector(threshold_plot_T4$phenotype_id))

pdf(file=args[12], width = 12, height = 12, onefile=FALSE)
upset(fromList(listInput), order.by = "freq", sets.bar.color = c(my_palette),
sets.x.label = "cis-eGenes tested", point.size = 4, line.size = 2,
mainbar.y.label = "cis-eGenes tested intersections",
text.scale = 2.5, shade.alpha = 0.5)


T0_nominal <- T0_nominal %>% filter(phenotype_id != "phenotype_id")
T1_nominal <- T1_nominal %>% filter(phenotype_id != "phenotype_id")
T2_nominal <- T2_nominal %>% filter(phenotype_id != "phenotype_id")
T3_nominal <- T3_nominal %>% filter(phenotype_id != "phenotype_id")
T4_nominal <- T4_nominal %>% filter(phenotype_id != "phenotype_id")


data_T0$is_eGene <- as.character(data_T0$is_eGene)
data_T1$is_eGene <- as.character(data_T1$is_eGene)
data_T2$is_eGene <- as.character(data_T2$is_eGene)
data_T3$is_eGene <- as.character(data_T3$is_eGene)
data_T4$is_eGene <- as.character(data_T4$is_eGene)


data_T0 <- data_T0 %>% filter(is_eGene == "TRUE")
data_T1 <- data_T1 %>% filter(is_eGene == "TRUE")
data_T2 <- data_T2 %>% filter(is_eGene == "TRUE")
data_T3 <- data_T3 %>% filter(is_eGene == "TRUE")
data_T4 <- data_T4 %>% filter(is_eGene == "TRUE")



listInput <- list(T0 = as.vector(data_T0$phenotype_id), 
T1 = as.vector(data_T1$phenotype_id), 
T2 = as.vector(data_T2$phenotype_id), 
T3 = as.vector(data_T3$phenotype_id), 
T4 = as.vector(data_T4$phenotype_id))

print(listInput)
pdf(file=args[19], width = 12, height = 12, onefile=FALSE)
upset(fromList(listInput), order.by = "freq", sets.bar.color = c(my_palette),
sets.x.label = "cis-eGenes tested", point.size = 4, line.size = 2,
mainbar.y.label = "cis-eGenes tested intersections",
text.scale = 2.5, shade.alpha = 0.5)
dev.off()


data_T0 <- left_join(data_T0, T0_nominal, by = c("phenotype_id" = "phenotype_id"))
data_T1 <- left_join(data_T1, T1_nominal, by = c("phenotype_id" = "phenotype_id"))
data_T2 <- left_join(data_T2, T2_nominal, by = c("phenotype_id" = "phenotype_id"))
data_T3 <- left_join(data_T3, T3_nominal, by = c("phenotype_id" = "phenotype_id"))
data_T4 <- left_join(data_T4, T4_nominal, by = c("phenotype_id" = "phenotype_id"))



data_T0 <- data_T0 %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))
data_T1 <- data_T1 %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))
data_T2 <- data_T2 %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))
data_T3 <- data_T3 %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))
data_T4 <- data_T4 %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))

rm(T0_nominal, T1_nominal, T2_nominal,T3_nominal, T4_nominal)

number_T0 <- length(data_T0$variant_id.y)
unique_T0 <- length(unique(data_T0$variant_id.y))
number_T1 <- length(data_T1$variant_id.y)
unique_T1 <- length(unique(data_T1$variant_id.y))
number_T2 <- length(data_T2$variant_id.y)
unique_T2 <- length(unique(data_T2$variant_id.y))
number_T3 <- length(data_T3$variant_id.y)
unique_T3 <- length(unique(data_T3$variant_id.y))
number_T4 <- length(data_T4$variant_id.y)
unique_T4 <- length(unique(data_T4$variant_id.y))

head(data_T0)
data_T0 <- fread(args[1]) %>% filter(is_eGene == TRUE) %>% mutate(variant_gene = paste0(variant_id, "-", phenotype_id), Category = "T0")
data_T1 <- fread(args[2]) %>% filter(is_eGene == TRUE) %>% mutate(variant_gene = paste0(variant_id, "-", phenotype_id), Category = "T1")
data_T2 <- fread(args[3]) %>% filter(is_eGene == TRUE) %>% mutate(variant_gene = paste0(variant_id, "-", phenotype_id), Category = "T2")
data_T3 <- fread(args[4]) %>% filter(is_eGene == TRUE) %>% mutate(variant_gene = paste0(variant_id, "-", phenotype_id), Category = "T3")
data_T4 <- fread(args[5]) %>% filter(is_eGene == TRUE) %>% mutate(variant_gene = paste0(variant_id, "-", phenotype_id), Category = "T4")

head(data_T0)
gtex <- fread("/home/workspace/jogrady/heQTL/data/gtex/Whole_Blood.signifpairs.txt")
gtex$variant_id <- gsub("_b38", "", gtex$variant_id)
gtex$variant_id <- gsub("_", ":", gtex$variant_id)
gtex$variant_id <- gsub("chr", "", gtex$variant_id)
gtex <- gtex %>% mutate(variant_gene = paste0(variant_id, "-", gene_id))

gtex <- gtex %>% select(variant_gene, slope)
colnames(gtex)[2] <- "gtex_slope"


data_T0 <- left_join(data_T0, gtex)
data_T0 <- data_T0 %>% filter(!is.na(gtex_slope))
data_T1 <- left_join(data_T1, gtex)
data_T2 <- left_join(data_T2, gtex)
data_T3 <- left_join(data_T3, gtex)
data_T4 <- left_join(data_T4, gtex)


data <- rbind(data_T0, data_T1, data_T2, data_T3, data_T4)


head(gtex)
head(data_T0)
ggplot(data = data, aes(x = slope, y = gtex_slope, col = Category, shape = Category)) + geom_point() + scale_colour_manual(values = my_palette) + geom_smooth(method = "lm")

head(data_T0)
cor(data_T0$slope, data_T0$gtex_slope)

write.table("Number of signicant cis-eQTLs in T0 group", file = args[13], sep = "\t", col.names = F, append = TRUE)
write.table(number_T0, file = args[13], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in T0 group", file = args[13], sep = "\t", col.names = F, append = T)
write.table(unique_T0, file = args[13], sep = "\t", col.names = F, append = T)

write.table("Number of signicant cis-eQTLs in T1 group", file = args[13], sep = "\t", col.names = F, append = T )
write.table(number_T1, file = args[13], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in T1 group", file = args[13], sep = "\t", col.names = F, append = T)
write.table(unique_T1, file = args[13], sep = "\t", col.names = F, append = T)

write.table("Number of signicant cis-eQTLs in T2 group", file = args[13], sep = "\t", col.names = F, append = T )
write.table(number_T2, file = args[13], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in T2 group", file = args[13], sep = "\t", col.names = F, append = T)
write.table(unique_T2, file = args[13], sep = "\t", col.names = F, append = T)

write.table("Number of signicant cis-eQTLs in T3 group", file = args[13], sep = "\t", col.names = F, append = T )
write.table(number_T3, file = args[13], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in T3 group", file = args[13], sep = "\t", col.names = F, append = T)
write.table(unique_T3, file = args[13], sep = "\t", col.names = F, append = T)

write.table("Number of signicant cis-eQTLs in T4 group", file = args[13], sep = "\t", col.names = F, append = T)
write.table(number_T4, file = args[13], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in T4 group", file = args[13], sep = "\t", col.names = F, append = T)
write.table(unique_T4, file = args[13], sep = "\t", col.names = F, append = T)


# write files with cis-eQTL pairs for TWAS which are significant
write.table(data_T0, file = args[14], sep = "\t", col.names = T, quote = F, row.names = F)
write.table(data_T1, file = args[15],sep = "\t", col.names = T, quote = F, row.names = F)
write.table(data_T2, file = args[16],sep = "\t", col.names = T, quote = F, row.names = F)
write.table(data_T3, file = args[17],sep = "\t", col.names = T, quote = F, row.names = F)
write.table(data_T4, file = args[18],sep = "\t", col.names = T, quote = F, row.names = F)

