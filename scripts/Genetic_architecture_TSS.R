# Genomic architecture of eQTLs and DR eQTLs

library(tidyverse)
library(ggplot2)
library(data.table)

chromosomes <- 1:22
groups <- c("T0", "T1", "T2", "T3", "T4")
data_null <- data.frame(matrix(nrow=0, ncol = 11))

colnames(data_null) <- c("phenotype_id", "variant_id", "start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se", "category", "group")   
for(group in groups) {
  for(chr in chromosomes) {
    data_temp <- fread(paste0("~/heQTL/results/eQTL/", group, ".50000.cis_qtl_pairs.chr",chr,".txt.gz"))
    data_temp$group <- "group"
    data_temp$category <- "ALL"
    data_null <- rbind(data_null, data_temp)
  }
}

data_eQTL <- data.frame(matrix(nrow = 0, ncol = 11))
colnames(data_eQTL) <- colnames(data_null)

head(data_eQTL)
for(group in groups) {
  data_temp <- fread(paste0("~/heQTL/results/eQTL/",group,".50000.cis_qtl_fdr0.1.txt"))
  data_temp$category <- "eQTL"
  print(head(data_temp))
  data_temp <- data_temp %>% filter(is_eGene == TRUE)
  data_temp$group <- group
  data_temp <- data_temp %>% select(colnames(data_eQTL))
  data_eQTL <- rbind(data_eQTL, data_temp)
}


dr_eQTL <- fread("/home/workspace/jogrady/heQTL/results/reQTLs/DR_eQTLs_MASH_1.5_lfsr_0.1.txt") %>% select(2)
rownames(dr_eQTL) <- dr_eQTL$x
dr_eQTL <- dr_eQTL %>% separate(.,x, into = c("variant_id", "phenotype_id"), sep = "-")

dr_eQTL <- left_join(dr_eQTL, data_null)

dr_eQTL <- dr_eQTL[!duplicated(dr_eQTL$variant_id),]

dr_eQTL$category <- "Response"
dr_eQTL$group<- "Response"

library(ggpubr)
my_palette <- c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026", "steelblue")
head(all_ggplot2)
all_ggplot2 <- rbind(data_null, data_eQTL, dr_eQTL)
all_ggplot2$group <- factor(all_ggplot2$group, levels = c("T0", "T1", "T2", "T3", "T4", "Response"))
head(all_ggplot2)
gghistogram(
  all_ggplot2 %>% filter(category == "eQTL" | category == "Response"), x = "start_distance", y = "..density..", bins = 50, fill = "group",
  add_density = TRUE, clinetype = 1, linewidth = 8, alpha = 0.4
) + facet_wrap(~group)

phist <- gghistogram(
  all_ggplot2 %>% filter(category == "eQTL" | category == "Response"), x = "start_distance", rug = FALSE, bins = 99,
  fill = "group", alpha = 0.4) + scale_x_continuous(breaks = c(-50000, -25000, 0, 25000, 50000)) + scale_fill_manual(values = my_palette)


# 2. Create the density plot with y-axis on the right
# Remove x axis elements
pdensity <- ggdensity(
  all_ggplot2 %>% filter(category == "eQTL" | category == "Response"), x = "start_distance", 
  color= "group",
  alpha = 0
) + scale_colour_manual(values = my_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")

library(cowplot)
max(all_ggplot2$start_distance)
aligned_plots <- patchwork::align_patches(phist, pdensity)
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggsave("/home/workspace/jogrady/heQTL/results/reQTLs/Genetic_architecture_TSS_eQTLs_response.pdf", width = 12, height = 12, dpi = 600)
