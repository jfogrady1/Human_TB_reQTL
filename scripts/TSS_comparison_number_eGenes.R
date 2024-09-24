# Script to determine how many SNPs are close to the TSS of genes
# and how many eGenes are detected at each window

# libraries
library(data.table)
library(tidyverse)
library(dplyr)
library(vcfR)

tensor_colnames_egene = c("phenotype_id", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df",
                              "variant_id", "start_distance", "end_distance", "ma_samples", "ma_count", "af",
                              "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta", "qval",
                              "pval_nominal_threshold", "pval_adj_BH", "is_eGene", "timepoint", "distance")

results_df_egene <- as.data.frame(matrix(ncol = 24, nrow = 0))

colnames(results_df_egene) <- tensor_colnames_egene


# function to cycle through results and collate into single object
open_eGene = function(timepoints, distances){
  tensor_colnames_egene = c("phenotype_id", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df",
                            "variant_id", "start_distance", "end_distance", "ma_samples", "ma_count", "af",
                            "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta", "qval",
                            "pval_nominal_threshold", "pval_adj_BH", "is_eGene", "timepoint", "distance")
  
  results_df_egene <- as.data.frame(matrix(ncol = 24, nrow = 0))
  
  colnames(results_df_egene) <- tensor_colnames_egene
  for (t in timepoints) {
    for (d in distances) {
      data_temp = read.table(paste0("/home/workspace/jogrady/heQTL/results/eQTL/", t, ".",d,".cis_qtl_fdr0.1.txt"))
      data_temp$timepoint <- t
      data_temp$distance <- d
      results_df_egene <- rbind(results_df_egene, data_temp)
    }
  }
 return(results_df_egene)
}

egene_results_df <- open_eGene(c("T0", "T1", "T2", "T3","T4"), c("10000","20000","50000","100000","200000", "500000", "1000000"))
dim(egene_results_df)
View(egene_results_df)
eGene_number_total <- egene_results_df %>% filter(is_eGene == TRUE) %>% group_by(distance) %>% summarise(distance = unique(distance),number = length(unique(phenotype_id)))


eGene_number <- egene_results_df %>% filter(is_eGene == TRUE) %>% group_by(timepoint, distance) %>% summarise(number = n())
eGene_number_total <- egene_results_df %>% filter(is_eGene == TRUE) %>% group_by(factor(distance, c("10000","20000","50000","100000","200000", "500000", "1000000"))) %>% summarise(distance = unique(distance), total = length(unique(phenotype_id))) %>% select(distance, total)
eGene_number_total

my_palette = c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026")

eGene_number_total$timepoint = "total"
eGene_number <- ggplot(data = eGene_number, aes(x = factor(distance,  levels =c("10000","20000","50000","100000","200000", "500000", "1000000")) , y = number, group = timepoint)) + 
  geom_line(aes(color=timepoint)) + 
  geom_point(aes(color=timepoint)) + 
  geom_line(as.data.frame(eGene_number_total), mapping = aes(x = factor(distance,  levels =c("10000","20000","50000","100000","200000", "500000", "1000000")), y = total), linetype = 2) +
  geom_point(as.data.frame(eGene_number_total), mapping = aes(x = factor(distance,  levels =c("10000","20000","50000","100000","200000", "500000", "1000000")), y = total)) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Unique cis-eGenes", breaks = c(250, 500, 750, 1000, 1250, 1500, 1750, 2000), labels = c(250, 500, 750, 1000, 1250, 1500, 1750, 2000)),
    breaks = c(250, 500, 750, 1000, 1250, 1500, 1750, 2000)
  ) +
  scale_x_discrete(breaks =c("10000","20000","50000","100000","200000", "500000", "1000000"), labels =c("10","20","50","100","200", "500", "1000")) +
  scale_colour_manual(values = my_palette) +
  labs(y = "cis-eGenes per group",
       x = "Distance +/- TSS (kbps)",
       colour = "Timepoint") + theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size= 12, colour = "black"),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 10))




snps_tested = egene_results_df %>% group_by(factor(distance, c("10000","20000","50000","100000","200000", "500000", "1000000"))) %>% summarise(phenotype_id = phenotype_id, distance = unique(distance),number = num_var, timepoint = timepoint)
snps_tested_total <- egene_results_df %>% group_by(factor(distance, c("10000","20000","50000","100000","200000", "500000", "1000000"))) %>% summarise(distance = unique(distance), total = length(unique(phenotype_id))) %>% select(distance, total)
snps_tested_total$timepoint <- "total"
head(snps_tested)
variant_number <- ggplot(data = snps_tested, aes(x = factor(distance,  levels =c("10000","20000","50000","100000","200000", "500000", "1000000")) , y = number, fill = timepoint)) + 
  geom_boxplot(aes(fill=timepoint), outlier.colour = NA) +
  geom_point(as.data.frame(snps_tested_total), mapping = aes(x = factor(distance,  levels =c("10000","20000","50000","100000","200000", "500000", "1000000")), y = total)) +
  geom_line(as.data.frame(snps_tested_total), mapping = aes(x = factor(distance,  levels =c("10000","20000","50000","100000","200000", "500000", "1000000")), y = total, group = timepoint), linetype = 2) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Genes tested accross groups",
  )) +
  scale_fill_manual(values = c(my_palette, "black")) +
  scale_x_discrete(breaks =c("10000","20000","50000","100000","200000", "500000", "1000000"), labels =c("10","20","50","100","200", "500", "1000")) +
  labs(y = "Variants tested per gene",
       x = "Distance +/- TSS (kbps)",
       fill = "Timepoint") + theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size= 12, colour = "black"),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 10))
library(cowplot)

variant_number
prow <- plot_grid(
  variant_number + theme(legend.position = "none"),
  eGene_number + theme(legend.position = "none"),
  align = 'vh',
  labels = c("A", "B"),
  nrow = 1,
  hjust = -1
)

prow
legend <- get_legend(
  # create some space to the left of the legend
  eGene_number +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 5))) +
    labs(colour = "Timepoint") +
    theme(legend.position = "bottom")
)


plot_grid(prow, legend, ncol = 1, rel_heights = c(1,0.1))
ggsave("/home/workspace/jogrady/heQTL/results/eQTL/plotting/TSS_distance_comparison.pdf", width = 12, height = 10, dpi = 600)
          
eGene_number
ggsave("/home/workspace/jogrady/heQTL/results/eQTL/plotting/TSS_distance_comparison_eGene_only.pdf", width = 12, height = 12, dpi = 600)
