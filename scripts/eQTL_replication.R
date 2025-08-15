# Replication rates human TB eQTL study and power caluclation given dfs and smaple size

library(tidyverse)
library(ggplot2)
library(powerEQTL)
library(ggpubr)
library(data.table)
library(magrittr)
args = commandArgs(trailingOnly = TRUE)

my_palette = c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026")
# read in sig genes
T0_sig = fread(args[1]) %>% filter(is_eGene == TRUE) %>% mutate(Time = "T0")
T1_sig = fread(args[2]) %>% filter(is_eGene == TRUE) %>% mutate(Time = "T1")
T2_sig = fread(args[3]) %>% filter(is_eGene == TRUE) %>% mutate(Time = "T2")
T3_sig = fread(args[4]) %>% filter(is_eGene == TRUE) %>% mutate(Time = "T3")
T4_sig = fread(args[5]) %>% filter(is_eGene == TRUE) %>% mutate(Time = "T4")

median(abs(T0_sig$slope)) # 0.784272
median(abs(T1_sig$slope)) # 0.749377
median(abs(T2_sig$slope)) # 0.768171
median(abs(T3_sig$slope)) # 0.79004
median(abs(T4_sig$slope)) # 0.794482
head(T0_sig)

hist(T0_sig$true_df)
 
all_effect_sizes = rbind(T0_sig, T1_sig, T2_sig, T3_sig, T4_sig)
all_effect_sizes$slope2 <- abs(all_effect_sizes$slope)
my_comparisons = list(c("T0", "T1"),c("T0", "T2"),c("T0", "T3"),c("T0", "T4"),c("T1", "T2"),c("T1", "T3"),c("T1", "T4"),c("T2", "T3"),c("T2", "T4"),c("T3", "T4"))

wilcox_test_internal <- compare_means(slope2 ~ Time, comparisons = my_comparisons, p.adjust.method = "BH", method='wilcox.test', data = all_effect_sizes)
wilcox_test_internal <- wilcox_test_internal %>% mutate(y.position = seq(2.6, 4.85, 0.25)) %>% mutate(Time = group1)



ggplot(all_effect_sizes, aes(x = Time, y = slope2, fill = Time)) + geom_violin(outlier.colour = NA, alpha = 0.8) +  
  geom_boxplot(outlier.colour = NA, width = 0.3, fill = "darkgrey") + 
  scale_fill_manual(values = my_palette) + theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.position = "none",
        axis.text.y = element_text(colour = "black", size = 15)) + stat_pvalue_manual(
          data = wilcox_test_internal,
          label = "p.adj",) + labs(y = "|Effect size (slope)|", x = "Timepoint")
ggsave(args[6], width = 15, height = 12)

############################
# Internal comparison ######
############################


# Read in sum stats
T0 <- fread(args[7])
T1 <- fread(args[8]) 
T2 <- fread(args[9])
T3 <- fread(args[10])  
T4 <- fread(args[11]) 







T0 <- T0 %>% filter(phenotype_id %in% T0_sig$phenotype_id) %>% left_join(T0_sig[,c("phenotype_id", "pval_nominal_threshold")]) %>% filter(as.numeric(pval_nominal) < as.numeric(pval_nominal_threshold)) %>% mutate(Time = "T0")
T1 <- T1 %>% filter(phenotype_id %in% T1_sig$phenotype_id) %>% left_join(T1_sig[,c("phenotype_id", "pval_nominal_threshold")]) %>% filter(as.numeric(pval_nominal) < as.numeric(pval_nominal_threshold)) %>% mutate(Time = "T1")
T2 <- T2 %>% filter(phenotype_id %in% T2_sig$phenotype_id) %>% left_join(T2_sig[,c("phenotype_id", "pval_nominal_threshold")]) %>% filter(as.numeric(pval_nominal) < as.numeric(pval_nominal_threshold)) %>% mutate(Time = "T2")
T3 <- T3 %>% filter(phenotype_id %in% T3_sig$phenotype_id) %>% left_join(T3_sig[,c("phenotype_id", "pval_nominal_threshold")]) %>% filter(as.numeric(pval_nominal) < as.numeric(pval_nominal_threshold)) %>% mutate(Time = "T3")
T4 <- T4 %>% filter(phenotype_id %in% T4_sig$phenotype_id) %>% left_join(T4_sig[,c("phenotype_id", "pval_nominal_threshold")]) %>% filter(as.numeric(pval_nominal) < as.numeric(pval_nominal_threshold)) %>% mutate(Time = "T4")


head(T0)
# Allelic concordance

files <- c(args[12:33])
print(files)
INTERVAL <- lapply(files,function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = TRUE)
})
INTERVAL <- bind_rows(INTERVAL)
head(INTERVAL)
INTERVAL$pair = paste0(INTERVAL$phenotype_id, "-", INTERVAL$chr, ":", INTERVAL$pos_b38)
colnames(INTERVAL)[7:9] <- c("INTERVAL_p_nominal", "INTERVAL_slope", "INTERVAL_slope_se")
colnames(INTERVAL)[12:13] <- c("INTERVAL_effect_allele", "INTERVAL_other_allele")

INTERVAL_signif = read.table(args[34], header = T) %>% filter(qval_sig == TRUE) %>% select(phenotype_id, pval_nominal_threshold)
colnames(INTERVAL_signif)[2] <- "INTERVAL_pval_threshold"


T0$phenotype_id <- gsub("\\..*", "", T0$phenotype_id)
T0$Pair = paste0(T0$phenotype_id, "-", T0$variant_id)

T1$phenotype_id <- gsub("\\..*", "", T1$phenotype_id)
T1$Pair = paste0(T1$phenotype_id, "-", T1$variant_id)



T2$phenotype_id <- gsub("\\..*", "", T2$phenotype_id)
T2$Pair = paste0(T2$phenotype_id, "-", T2$variant_id)


T3$phenotype_id <- gsub("\\..*", "", T3$phenotype_id)
T3$Pair = paste0(T3$phenotype_id, "-", T3$variant_id)


T4$phenotype_id <- gsub("\\..*", "", T4$phenotype_id)
T4$Pair = paste0(T4$phenotype_id, "-", T4$variant_id)



T0_sig$phenotype_id <- gsub("\\..*", "", T0_sig$phenotype_id)
T0_sig$Pair = paste0(T0_sig$phenotype_id, "-", T0_sig$variant_id)


T1_sig$phenotype_id <- gsub("\\..*", "", T1_sig$phenotype_id)
T1_sig$Pair = paste0(T1_sig$phenotype_id, "-", T1_sig$variant_id)



T2_sig$phenotype_id <- gsub("\\..*", "", T2_sig$phenotype_id)
T2_sig$Pair = paste0(T2_sig$phenotype_id, "-", T2_sig$variant_id)


T3_sig$phenotype_id <- gsub("\\..*", "", T3_sig$phenotype_id)
T3_sig$Pair = paste0(T3_sig$phenotype_id, "-", T3_sig$variant_id)


T4_sig$phenotype_id <- gsub("\\..*", "", T4_sig$phenotype_id)
T4_sig$Pair = paste0(T4_sig$phenotype_id, "-", T4_sig$variant_id)



AC_INTERVAL = rbind(T0_sig, T1_sig, T2_sig, T3_sig, T4_sig)
AC_INTERVAL <- AC_INTERVAL %>% separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = ":")
AC_INTERVAL$pair <- paste0(AC_INTERVAL$phenotype_id, "-", AC_INTERVAL$chr, ":", AC_INTERVAL$pos)



AC_INTERVAL <- left_join(AC_INTERVAL, INTERVAL[,c("pair","INTERVAL_p_nominal", "INTERVAL_slope", "INTERVAL_slope_se","INTERVAL_effect_allele", "INTERVAL_other_allele")])
AC_INTERVAL$match = if_else(AC_INTERVAL$ref == AC_INTERVAL$INTERVAL_effect_allele, "TRUE", "FALSE")
AC_INTERVAL <- AC_INTERVAL %>% filter(!is.na(match)) %>% filter(match == "TRUE") %>% mutate("INTERVAL_new_slope" = INTERVAL_slope * -1) # need to change as effect allele is reference
AC_INTERVAL <- left_join(AC_INTERVAL, INTERVAL_signif)
AC_INTERVAL <- AC_INTERVAL %>% filter(INTERVAL_p_nominal < INTERVAL_pval_threshold)




dim(INTERVAL_signif)
table(AC_INTERVAL$Time)
#T0  T1  T2  T3  T4 
#443 436 487 451 464  

AC_INTERVAL <- AC_INTERVAL %>% mutate(AC = case_when(
                                                     slope > 0 & INTERVAL_new_slope > 0 ~ "CONSISTENT",
                                                     slope < 0 & INTERVAL_new_slope < 0 ~ "CONSISTENT",
                                                     slope > 0 & INTERVAL_new_slope < 0 ~ "INCONSISTENT",
                                                     slope < 0 & INTERVAL_new_slope > 0 ~ "INCONSISTENT"))

AC_INTERVAL

FINAL_INTERVAL_AC <- AC_INTERVAL %>%
  group_by(Time) %>%
  summarise(Proportion_Consistent = mean(AC == "CONSISTENT"),
            Count_Consistent = sum(AC == "CONSISTENT"),
            Total_Count = n()) %>% mutate(Study = "INTERVAL")

FINAL_INTERVAL_AC
library(qvalue)

### pi1

AC_INTERVAL = rbind(T0_sig, T1_sig, T2_sig, T3_sig, T4_sig)
AC_INTERVAL <- AC_INTERVAL %>% separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = ":")
AC_INTERVAL$pair <- paste0(AC_INTERVAL$phenotype_id, "-", AC_INTERVAL$chr, ":", AC_INTERVAL$pos)
AC_INTERVAL <- left_join(AC_INTERVAL, INTERVAL[,c("pair","INTERVAL_p_nominal", "INTERVAL_slope", "INTERVAL_slope_se","INTERVAL_effect_allele", "INTERVAL_other_allele")])
AC_INTERVAL$match = if_else(AC_INTERVAL$ref == AC_INTERVAL$INTERVAL_effect_allele, "TRUE", "FALSE")
AC_INTERVAL <- AC_INTERVAL %>% filter(!is.na(match)) %>% filter(match == "TRUE") %>% mutate("INTERVAL_new_slope" = INTERVAL_slope * -1) # need to change as effect allele is reference
AC_INTERVAL <- left_join(AC_INTERVAL, INTERVAL_signif)

head(AC_INTERVAL)


table(AC_INTERVAL$Time)

#T0  T1  T2  T3  T4 
#480 475 516 494 508 

FINAL_INTERVAL_pi1 <- AC_INTERVAL %>% group_by(Time) %>% arrange(desc(pval_nominal)) %>% summarise(pi0 = pi0est(INTERVAL_p_nominal)$pi0,
  pi1 = 1 - pi0est(INTERVAL_p_nominal)$pi0) %>% mutate(Study = "INTERVAL")


FINAL_INTERVAL_pi1

# T0    0.0563 0.944 INTERVAL
# T1    0.0543 0.946 INTERVAL
# T2    0.0359 0.964 INTERVAL
# T3    0.0310 0.969 INTERVAL
# T4    0.0806 0.919 INTERVAL

AC_INTERVAL = rbind(T0_sig, T1_sig, T2_sig, T3_sig, T4_sig)
AC_INTERVAL <- AC_INTERVAL %>% separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = ":")
AC_INTERVAL$pair <- paste0(AC_INTERVAL$phenotype_id, "-", AC_INTERVAL$chr, ":", AC_INTERVAL$pos)
AC_INTERVAL <- left_join(AC_INTERVAL, INTERVAL[,c("pair","INTERVAL_p_nominal", "INTERVAL_slope", "INTERVAL_slope_se","INTERVAL_effect_allele", "INTERVAL_other_allele")])
AC_INTERVAL$match = if_else(AC_INTERVAL$ref == AC_INTERVAL$INTERVAL_effect_allele, "TRUE", "FALSE")
AC_INTERVAL <- AC_INTERVAL %>% filter(!is.na(match)) %>% filter(match == "TRUE") %>% mutate("INTERVAL_new_slope" = INTERVAL_slope * -1) # need to change as effect allele is reference
AC_INTERVAL <- left_join(AC_INTERVAL, INTERVAL_signif)
head(AC_INTERVAL)

dim(AC_INTERVAL)

FINAL_INTERVAL_cor = AC_INTERVAL %>% group_by(Time) %>% summarize(Rho = as.numeric(cor.test(slope, INTERVAL_new_slope, method = "spearman", exact = FALSE)$estimate),
                                             P = as.numeric(cor.test(slope, INTERVAL_new_slope, method = "spearman", exact = FALSE)$p.value)) %>% mutate(Study = "INTERVAL")

# T0    0.798 3.43e-107 INTERVAL
# T1    0.806 5.47e-110 INTERVAL
# T2    0.827 1.73e-130 INTERVAL
# T3    0.812 2.20e-117 INTERVAL
# T4    0.827 9.69e-129 INTERVAL



# Plot the Correlaiton
ggplot(data = AC_INTERVAL, aes(x = slope, y = INTERVAL_new_slope, colour = Time)) + 
  geom_point(alpha = 0.7, size = 2) + scale_color_manual(values = my_palette) +
  theme_bw() +
  labs(colour = "eQTL cohort", x = "Effect size (slope) in discovery cohort", y = "Effect size (slope) in INTERVAL cohort") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "top")
ggsave(args[35], width = 12, height = 12)






##########################
# MASHR comparison
##########################
mashr = read.table(args[36])
mashr$gene <- gsub("\\..*", "", mashr$gene)
mashr <- mashr %>% separate(., variant, into = c("chr", "pos", "ref", "alt"))
mashr$pair <- paste0(mashr$gene, "-", mashr$chr,":",mashr$pos)
mashr <- left_join(mashr, INTERVAL[,c("pair","INTERVAL_p_nominal", "INTERVAL_slope", "INTERVAL_slope_se","INTERVAL_effect_allele", "INTERVAL_other_allele")])
mashr$match = if_else(mashr$ref == mashr$INTERVAL_effect_allele, "TRUE", "FALSE")
mashr <- mashr %>% filter(!is.na(match)) %>% filter(match == "TRUE") %>% mutate("INTERVAL_new_slope" = INTERVAL_slope * -1) # need to change as effect allele is reference

head(mashr)
T0_mashr <- mashr %>% filter(T0 < 0.05) %>% mutate(Time = "T0") %>% select(T0_posterior_mean, INTERVAL_new_slope, gene_name, variant_gene, Time) %>% set_colnames(.,c("Posterior_mean", "External_slope", "gene", "pair", "Time"))
T1_mashr <- mashr %>% filter(T1 < 0.05) %>% mutate(Time = "T1") %>% select(T1_posterior_mean, INTERVAL_new_slope, gene_name, variant_gene, Time) %>% set_colnames(.,c("Posterior_mean", "External_slope", "gene", "pair", "Time"))
T2_mashr <- mashr %>% filter(T2 < 0.05) %>% mutate(Time = "T2") %>% select(T2_posterior_mean, INTERVAL_new_slope, gene_name, variant_gene, Time) %>% set_colnames(.,c("Posterior_mean", "External_slope", "gene", "pair", "Time"))
T3_mashr <- mashr %>% filter(T3 < 0.05) %>% mutate(Time = "T3") %>% select(T3_posterior_mean, INTERVAL_new_slope, gene_name, variant_gene, Time) %>% set_colnames(.,c("Posterior_mean", "External_slope", "gene", "pair", "Time"))
T4_mashr <- mashr %>% filter(T4 < 0.05) %>% mutate(Time = "T4") %>% select(T4_posterior_mean, INTERVAL_new_slope, gene_name, variant_gene, Time) %>% set_colnames(.,c("Posterior_mean", "External_slope", "gene", "pair", "Time"))



mashr_plot <- rbind(T0_mashr,T1_mashr,T2_mashr,T3_mashr,T4_mashr)

head(mashr_plot)

INTERVAL_mashr_cor = mashr_plot %>% group_by(Time) %>% summarize(Rho = as.numeric(cor.test(Posterior_mean, External_slope, method = "spearman", exact = FALSE)$estimate),
                                                                  P = as.numeric(cor.test(Posterior_mean, External_slope, method = "spearman", exact = FALSE)$p.value)) %>% mutate(Study = "INTERVAL")





INTERVAL_mashr_cor
#Time    Rho     P Study   
# T0    0.810     0 INTERVAL
# T1    0.814     0 INTERVAL
# T2    0.817     0 INTERVAL
# T3    0.809     0 INTERVAL
# T4    0.808     0 INTERVAL
mashr_plot$diff = abs(mashr_plot$Posterior_mean - mashr_plot$External_slope)
ggplot(mashr_plot, aes(x = Posterior_mean, y = External_slope, col = Time)) + geom_point(alpha = 0.2, size = 2) + scale_color_manual(values = my_palette) +
  theme_bw() +
  labs(colour = "eQTL cohort", x = "Posterior mean effect size in discovery cohort", y = "Effect size (slope) in INTERVAL cohort") +
  geom_smooth(method = "lm", se = FALSE) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "top") + lims(x = c(-2,2))
ggsave(args[37], width = 12, height = 12)



library(pwr)
pwr.f2.test(u = 1, v = 48 - 12 - 1, sig.level = 0.0001, power = 0.8) # T0
pwr.f2.test(u = 1, v = 48 - 14 - 1, sig.level = 0.0001, power = 0.8) # T1
pwr.f2.test(u = 1, v = 48 - 13 - 1, sig.level = 0.0001, power = 0.8) # T3
pwr.f2.test(u = 1, v = 48 - 13 - 1, sig.level = 0.0001, power = 0.8) # T3
pwr.f2.test(u = 1, v = 48 - 15 - 1, sig.level = 0.0001, power = 0.8) # T4

