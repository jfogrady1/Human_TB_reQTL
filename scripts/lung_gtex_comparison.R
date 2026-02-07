# Script to compare lung GTEX eQTLs to those from the our study

library(tidyverse)
library(data.table)
library(arrow)
args = CompandArgs(trailingOnly=TRUE)
lung_gtex <- read_parquet(args[1])

lung_gtex$phenotype_id <- gsub("\\..*", "", lung_gtex$phenotype_id)
lung_gtex$variant_id <- gsub("chr", "", lung_gtex$variant_id)
lung_gtex$variant_id <- gsub("_", ":", lung_gtex$variant_id)
lung_gtex$variant_id <- gsub(":b38", "", lung_gtex$variant_id)  

lfsr_resutls <- fread(args[2])
lfsr_resutls$phenotype_id <- gsub("\\..*", "", lfsr_resutls$gene)


head(lfsr_resutls)

lung_gtex$pair <- paste(lung_gtex$variant_id, lung_gtex$phenotype_id, sep = "_")
lfsr_resutls$pair <- paste(lfsr_resutls$variant, lfsr_resutls$phenotype_id, sep = "_")
lfsr_resutls$pair

lfsr_lung_gtex <- lfsr_resutls %>%
  filter(pair %in% lung_gtex$pair)
dim(lfsr_lung_gtex)

head(lfsr_lung_gtex)

# detect columns - lfsr and slope/posterior means
sig_cols <- grep("^T[0-9]+$", names(lfsr_lung_gtex), value = TRUE)
post_cols <- grep("^T[0-9]+_posterior_mean$", names(lfsr_lung_gtex), value = TRUE)

# pivot LFSR (significance) to long
lfsr_lung_gtex_long <- lfsr_lung_gtex %>%
  select(variant_gene, pair, phenotype_id, all_of(sig_cols)) %>%
  pivot_longer(cols = all_of(sig_cols), names_to = "time", values_to = "lfsr")

# pivot posterior means and normalize time names
post_long <- lfsr_lung_gtex %>%
  select(variant_gene, pair, phenotype_id, all_of(post_cols)) %>%
  pivot_longer(cols = all_of(post_cols), names_to = "time", values_to = "posterior_mean") %>%
  mutate(time = sub("_posterior_mean$", "", time))

# join them
lfsr_lung_gtex_long <- lfsr_lung_gtex_long %>%
  left_join(post_long, by = c("variant_gene", "pair", "phenotype_id", "time")) %>%
  mutate(timepoint = as.integer(sub("^T", "", time)))

lfsr_lung_gtex_long <- lfsr_lung_gtex_long %>% left_join(., lung_gtex, by = "pair") 


my_palette = c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026")



  
  ggplot(aes(x = posterior_mean, y = slope)) +
  geom_point() +
  facet_wrap(~time) +
  theme_bw()


lfsr_lung_gtex_long %>%
  filter(lfsr < 0.05) %>% group_by(time) %>% 
  summarise(correlation = cor(posterior_mean, slope, method = 'spearman'),
            p_val = cor.test(posterior_mean, slope, method = 'spearman')$p.value,
            variant_n = n())


#time  correlation p_val variant_n
#T0          0.824     0      1382
#T1          0.826     0      1395
#T2          0.824     0      1370
#T3          0.823     0      1370
#T4          0.825     0      1403

  
lfsr_lung_gtex_long %>%
  filter(lfsr < 0.05) %>%
  ggplot(aes(x = posterior_mean, y = slope, col = time)) +
  geom_point(alpha = 0.7, size =2) +
  scale_colour_manual(values = my_palette) +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, col = "grey60", size = 2, linetype = "dashed") +
  facet_wrap(~time)
ggsave("/home/workspace/jogrady/heQTL/results/Response_1/Figures/lung_gtex_comparison.pdf", width = 12, height = 12, dpi = 600)


# Compare effect directions: study posterior_mean vs GTEx slope
direction_comparison <- lfsr_lung_gtex_long %>%
  mutate(
    study_direction = sign(posterior_mean),
    gtex_direction = sign(slope)
  ) %>%
  mutate(
    same_direction = study_direction == gtex_direction & study_direction != 0
  )

View(direction_comparison)
# summary by time point (lfsr < 0.05)
direction_summary_by_time <- direction_comparison %>%
  filter(lfsr < 0.05, study_direction != 0, gtex_direction != 0) %>%
  group_by(time) %>%
  summarise(
    n_variants = n(),
    n_concordant = sum(same_direction),
    prop_concordant = n_concordant / n_variants,
    .groups = "drop"
  ) %>%
  arrange(time)

direction_summary_by_time
#time  n_variants n_concordant prop_concordant
#1 T0          1382         1303           0.943
#2 T1          1395         1315           0.943
#3 T2          1370         1294           0.945
#4 T3          1370         1291           0.942
#5 T4          1403         1321           0.942