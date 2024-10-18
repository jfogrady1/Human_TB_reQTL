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
data <- fread(args[1])
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




rownames(data) <- data$Mixture

transformed_data <- t(data)
colnames(transformed_data) <- transformed_data[1,]
transformed_data <- transformed_data[-1,]
row_cell_types <- rownames(transformed_data)



transformed_data <- apply(transformed_data, 2, as.numeric)
prop.list <- convertDataToList(transformed_data,
                               data.type="proportions", 
                               transform="asin")


prop.list$TransformedProps
rownames(prop.list$TransformedProps) <- row_cell_types
rownames(prop.list$Counts) <- row_cell_types
rownames(prop.list$Proportions) <- row_cell_types
group <- factor(c(rep("T0", 48 ), rep("T1", 48), rep ("T2", 48), rep("T3", 48), rep("T4", 48)), labels = c("T0","T1", "T2", "T3", "T4"))
sample <- factor(gsub("_T.", "", colnames(transformed_data)))
data.frame(group,sample)

# Fit sample ID as a random effect for intra-block correlation assessment
des.tech <- model.matrix(~group)
des.tech
dupcor <- duplicateCorrelation(prop.list$TransformedProps, design = des.tech, block = sample)

dupcor


fit1 <- lmFit(prop.list$TransformedProps, design=des.tech, block=sample, 
              correlation=dupcor$consensus)

fit1 <- eBayes(fit1)
summary(decideTests(fit1))
topTable(fit1, coef=2:5)
Limma_results <- topTable(fit1, coef=c(2:5),number = 16, adjust.method = "BH")

Limma_results
write.table(Limma_results, file = args[2], sep = "\t")
Limma_results$cell_types <- rownames(Limma_results)

Limma_results <- Limma_results %>% mutate(label = case_when(adj.P.Val < 0.05 ~ paste("Padj. < 0.05"),
                                                            adj.P.Val > 0.05 ~ paste("Padj. = ns")))

Limma_results$cell_types <- as.factor(Limma_results$cell_types)
#Limma_results <- Limma_results[order(Limma_results$cell_types),]

# Value for plotting
pval_proportion <- data_plot %>% group_by(cell_types) %>% summarize(proportion = max(proportion)) %>% select(cell_types, proportion)
pval_proportion <- data.frame(pval_proportion)
rownames(pval_proportion) <- pval_proportion$cell_types


#data_plot <- data_plot[order(data_plot$cell_types),]
head(data_plot)
proportion <- pval_proportion$proportion
proportion <- as.vector(proportion)
proportion
data_plot$cell_types
#data_plot <- data_plot[sort(data_plot$cell_types),]

Limma_results$Time = "T3"
Limma_results$cell_types <- as.factor(Limma_results$cell_types)
Limma_results <- left_join(Limma_results, pval_proportion)
Limma_results$proportion <- Limma_results$proportion
Limma_results$cell_types <- as.factor(Limma_results$cell_types)
#Limma_results <- Limma_results[order(Limma_results$cell_types),]

all(unique(data_plot$cell_types) == Limma_results$cell_types)


data_plot$cell_types <- as.factor(data_plot$cell_types)
data_plot <- data_plot[order(data_plot$cell_types),]

rownames(Limma_results) <- Limma_results$cell_types

Limma_results[unique(data_plot$cell_types),]
Limma_results <- Limma_results[order(Limma_results$cell_types),]
Limma_results$proportion <- Limma_results$proportion * 1.15

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
  labs(y = "Deconvolved cell-type proportion",
       x = "Time") +
  geom_text(x = "T3", y = Limma_results$proportion, aes(label = label), data = Limma_results, check_overlap = TRUE) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 1)), # Classical monocyte
      scale_y_continuous(limits = c(0, 0.03)), # Dendritic cell
      scale_y_continuous(limits = c(0, .03)), # MAIT cell
      scale_y_continuous(limits = c(0, 0.1)), # Megakaryocyte
      scale_y_continuous(limits = c(0, 0.03)), # Memory B cell
      scale_y_continuous(limits = c(0, .2)), # Memory CD4+ T cell
      scale_y_continuous(limits = c(0, 0.2)), # Memory CD8+ T cell
      scale_y_continuous(limits = c(0, 0.3)), # Naive B cell
      scale_y_continuous(limits = c(0, 0.125)), # Naive CD4+ T cell
      scale_y_continuous(limits = c(0, .1)), # Naive CD8+ T cell
      scale_y_continuous(limits = c(0, 0.4)), # NK cell
      scale_y_continuous(limits = c(0, .35)), # non classical monocyte
      scale_y_continuous(limits = c(0, 0.02)), # Plasma B cell
      scale_y_continuous(limits = c(0, .15)), # Regualtory T cell
      scale_y_continuous(limits = c(0, .65)) # Undefined T cell
    )
  )

ggsave(args[3], width = 12, height = 12, dpi = 600)
# Correlation to MuSic

counts <- fread(args[4])
rows <- counts[,1]
cols <- colnames(counts)[2:241]
counts <- counts[,-1]
counts <- as.matrix(counts)
rownames(counts) <- rows$GeneSymbol
head(counts)
counts <- counts[unique(rownames(counts)),]
counts<-new("ExpressionSet", exprs=as.matrix(counts))
counts <- exprs(counts)
markers = fread(args[5])

TB.sc.data = readRDS(args[6])
TB.sc.data$cell_type <- Idents(TB.sc.data)
TB.sc.data.sce <- as.SingleCellExperiment(TB.sc.data)

Est.prop.human.TB = music_prop(bulk.mtx = counts, sc.sce  = TB.sc.data.sce, markers = markers$gene, clusters = 'ident', samples = 'sample_id', verbose = T)
music_results <- as.data.frame(Est.prop.human.TB$Est.prop.weighted)
music_results$Mixture <- rownames(music_results)



data_ciber <- fread(args[2])



data_long <- data %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")
music_results_long <- music_results %>% pivot_longer(!Mixture, names_to = "cell_types", values_to = "proportion")

data_long <- data_long[order(data_long$cell_types),]
music_results_long <- music_results_long[order(music_results_long$cell_types),]


music_results_long$Type <- "MuSiC"
data_long$Type <- "Cibersort"
times <- rep(c(rep("T0", 48 ), rep("T1", 48), rep ("T2", 48), rep("T3", 48), rep("T4", 48)),15)
times
data_long$times <- times
music_results_long$times <- times

colnames(music_results_long) <- paste0(colnames(music_results_long), "_MUSIC")
colnames(music_results_long)
data_merged <- cbind(data_long, music_results_long[,c("proportion_MUSIC")])


ggplot(data = data_merged, aes(x = proportion, y = proportion_MUSIC, col = cell_types)) + geom_point() + ggh4x::facet_grid2(cell_types~times, scales = "free", independent = "x") + geom_smooth(method = "lm")

ggsave(args[7], width = 12, height = 12, dpi = 600)

data_merged %>% 
  group_by(cell_types, times) %>% 
  summarize(correlation = cor(proportion, proportion_MUSIC, method = "spearman"),
            pval = cor.test(proportion, proportion_MUSIC, method = "spearman", exact = FALSE)$p.value) %>% 
            mutate(padj = p.adjust(pval, method = "BH")) %>% mutate(shape = if_else(padj < 0.05, "Padj < 0.05", "Padj > 0.05")) %>%
  ungroup() %>%  # Ungroup to prevent any unintended side effects in the plotting
  ggplot(aes(x = times, y = as.numeric(correlation))) + 
    geom_boxplot(aes(fill = as.factor(times), alpha = 0.5)) +  # Ensure fill is mapped correctly, e.g., by times
    geom_jitter(width = 0.1, size = 4, alpha = 1, aes(col = cell_types, shape = shape)) +  # Add jitter with some customization
    scale_fill_manual(values = my_palette) +  # Apply your custom palette
    theme_minimal() +  # Optional: Improve the appearance
    labs(x = "Times", y = "Spearman Correlation") + theme_bw() +
    scale_shape_manual(labels = c("Padj < 0.05", "Padj > 0.05"), values = c(16, 17)) +
      theme(legend.position = "right",
        title = element_text(size = 16),
        legend.background = element_rect(color = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size= 14, colour = "black"),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10))
ggsave(args[8], width = 12, height = 12, dpi = 600)


data_table <- data_merged %>% 
  group_by(cell_types, times) %>% 
  summarize(correlation = cor(proportion, proportion_MUSIC, method = "spearman"),
            pval = cor.test(proportion, proportion_MUSIC, method = "spearman", exact = FALSE)$p.value) 
head(data_table)
write.table(data_table, args[9], sep = "\t", quote = FALSE)


data_merged %>% group_by(cell_types, times) %>% summarize(Median = median(proportion)) %>% filter(Median > 0.01) %>% select(cell_types) %>% table()
 #   Classical Monocyte     Memory CD4+ T cell     Memory CD8+ T cell 
  #                   5                      4                      5 
   #       Naïve B cell                NK cell Non-classical Monocyte 
    #                 5                      5                      5 
     #T regulatory cell       Undefined T cell 
      #               2                      5 

# Classical Monocyte, Memory CD8+ T cell, Naive B cell, NK cell, Non-Classical Monocyte,