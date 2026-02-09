# Multi variate adaptive shrinkage to identify response eQTLs


# MASHR

args <- commandArgs(trailingOnly = TRUE)
#### Load in packages

library(broom)
library(ashr)
library(mashr)
library(data.table)
library(tidyverse)
library(ggVennDiagram) 
library(ggvenn)
library(ggrepel)
library(gprofiler2)
library(viridis)
library(ggpubr)
library(cowplot)
library(broom)
library(vcfR)
library(rstatix)
set.seed(34567)

args[1] <- "T0"
args[2] <- "T1"
args[3] <- "T2"
args[4] <- "T3"
args[5] <- "T4"
args[6] <- "home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"

##### Define functions for and opening necessary files


tensor_colnames_egene = c("phenotype_id", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", # nolint: assignment_linter.
                          "variant_id", "start_distance", "end_distance", "ma_samples", "ma_count", "af",
                          "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta", "qval",
                          "pval_nominal_threshold", "pval_adj_BH", "is_eGene", "timepoint", "distance")

results_df_egene <- as.data.frame(matrix(ncol = 24, nrow = 0))

colnames(results_df_egene) <- tensor_colnames_egene

chromosomes = as.character(c(1:22))


# First read in all the data
# function to cycle through results and collate into single object
open_eGene = function(timepoints, chromosome){
  tensor_colnames_egene = c("phenotype_id", "variant_id", "start_distance", "af", "ma_samples", "ma_count", "af",
                            "pval_nominal", "slope", "slope_se", "timepoint")
  
  results_df_egene <- as.data.frame(matrix(ncol = 11, nrow = 0))
  
  colnames(results_df_egene) <- tensor_colnames_egene
  for (t in timepoints) {
    for (c in chromosomes) {
      data_temp = fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/", t, ".50000.cis_qtl_pairs.chr", c, ".txt.gz"))
      data_temp$timepoint <- t
      results_df_egene <- rbind(results_df_egene, data_temp)
    }
  }
  return(results_df_egene)
}

T0 <- open_eGene(args[1], chromosomes)
T1 <- open_eGene(args[2], chromosomes)
T2 <- open_eGene(args[3], chromosomes)
T3 <- open_eGene(args[4], chromosomes)
T4 <- open_eGene(args[5], chromosomes) # nolint
head(T4)
# Get common genes to all
common_genes = intersect(T0$phenotype_id, intersect(T1$phenotype_id, intersect(T2$phenotype_id, intersect(T3$phenotype_id, T4$phenotype_id))))
length(common_genes)
# 15109 common genes
T0 <- T0 %>% filter(phenotype_id %in% common_genes)
T1 <- T1 %>% filter(phenotype_id %in% common_genes)
T2 <- T2 %>% filter(phenotype_id %in% common_genes)
T3 <- T3 %>% filter(phenotype_id %in% common_genes)
T4 <- T4 %>% filter(phenotype_id %in% common_genes)                         

all(T0$phenotype_id == T1$phenotype_id)
all(T1$phenotype_id == T2$phenotype_id)
all(T2$phenotype_id == T3$phenotype_id)
all(T3$phenotype_id == T4$phenotype_id)


all(T0$variant_id == T1$variant_id)
all(T1$variant_id == T2$variant_id)
all(T2$variant_id == T3$variant_id)
all(T3$variant_id == T4$variant_id)




####### Getting everything set up for MASHR

# Now have all common SNPs tested against all common genes
# set up the mash data frame
simdata = simple_sims(10000,5,1) # simulates data on 40k tests


# get matrix of effect sizes
bhat_df <- cbind(T0$slope, T1$slope, T2$slope, T3$slope, T4$slope)
bhat_df_matrix <- as.matrix(bhat_df)
colnames(bhat_df_matrix) <- c("T0", "T1", "T2", "T3", "T4")
rownames(bhat_df_matrix) <- paste0(T0$variant_id,"-",T0$phenotype_id)



bhat_df_se <- cbind(T0$slope_se, T1$slope_se, T2$slope_se, T3$slope_se, T4$slope_se)
bhat_df_se_matrix <- as.matrix(bhat_df_se)
colnames(bhat_df_se_matrix) <- c("T0", "T1", "T2", "T3", "T4")
rownames(bhat_df_se_matrix) <- paste0(T0$variant_id,"-",T0$phenotype_id)


# function to cycle through results and collate into single object
open_eGene_signif = function(timepoints){
  tensor_colnames_egene = c("phenotype_id", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df",
                            "variant_id", "start_distance", "end_distance", "ma_samples", "ma_count", "af",
                            "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta", "qval",
                            "pval_nominal_threshold", "pval_adj_BH", "is_eGene", "timepoint")
  
  results_df_egene <- as.data.frame(matrix(ncol = 23, nrow = 0))
  
  colnames(results_df_egene) <- tensor_colnames_egene
  for (t in timepoints) {
      data_temp = read.table(paste0("/home/workspace/jogrady/heQTL/results/eQTL/", t, ".50000.cis_qtl_fdr0.1.txt"))
      data_temp$timepoint <- t
      results_df_egene <- rbind(results_df_egene, data_temp)
  }
  return(results_df_egene)
}


# Need two things, strong associations - lead snp per gene tested
# random associations = see below


# Here we will get the most significant true-eQTL SNP per gene accross conditions

signif_eGenes = open_eGene_signif(timepoints = c(args[1], args[2], args[3], args[4], args[5]))

signif_eGenes <- signif_eGenes %>% filter(phenotype_id %in% common_genes)

signif_eGenes #<- signif_eGenes %>% filter(is_eGene == TRUE)


signif_eGenes$variant_gene <- paste0(signif_eGenes$variant_id,"-",signif_eGenes$phenotype_id)

length(unique(signif_eGenes$phenotype_id)) # 1824

signif_variants = signif_eGenes %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(qval == min(qval)) # select most significant variant per gen



# Here are our lead unique SNps per gene accross all condiitons which are eQTLs
strong.subset = as.vector(signif_variants$variant_gene)
names(strong.subset) <- strong.subset


mashr_input <- list(Bhat = bhat_df_matrix,
                    Shat = bhat_df_se_matrix)



random.subset = sample(1:nrow(mashr_input$Bhat),200000)


# learn the canonical covariance matrices
data.temp = mash_set_data(mashr_input$Bhat[random.subset,],mashr_input$Shat[random.subset,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)


data.random = mash_set_data(mashr_input$Bhat[random.subset,],mashr_input$Shat[random.subset,],V=Vhat)
data.strong = mash_set_data(mashr_input$Bhat[strong.subset,],mashr_input$Shat[strong.subset,], V=Vhat)

U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)



U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)


m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

length(get_significant_results(m2, thresh = 0.05, conditions = NULL, sig_fn = get_lfsr))



get_pairwise_sharing(
  m2,
  factor = 2,
  lfsr_thresh = 0.05,
  FUN=abs
)

thresh <- 0.05

lfsr.mash <- m2$result$lfsr
#lfsr.mash
sigmat <- (lfsr.mash < thresh)
nsig <- rowSums(sigmat)
lfsr.mash.sig <- lfsr.mash[nsig > 0,]
dim(lfsr.mash.sig) 
apply(lfsr.mash<thresh, 2, sum)


####################### 
# Downstream stuff after MASHR
#######################

symbols <- fread(args[6]) %>% filter(V3 == "gene") # select(V9)
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)
symbols$length <- abs(symbols$V4 - symbols$V5)
symbols <- symbols %>% dplyr::select(gene_id, gene_name, length)

mas_signif_table <- lfsr.mash.sig

mas_signif_table <- as.data.frame(mas_signif_table)
mas_signif_table$variant_gene <- rownames(mas_signif_table)
mas_signif_table <- mas_signif_table %>% separate(variant_gene, into = c("variant", "gene"), sep = "-")

mas_signif_table <- left_join(mas_signif_table, symbols, by = c("gene" = "gene_id"))

mas_signif_table <- mas_signif_table %>% select(-c("length"))
mas_signif_table <- mas_signif_table %>% select(gene, gene_name, variant, T0, T1, T2, T3, T4)



m2_posterior = mash_compute_posterior_matrices(m,data.strong)

updated_means <- m2_posterior$PosteriorMean
updated_sds <- m2_posterior$PosteriorSD
T0_beta <- T0 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>% filter(variant_gene %in% rownames(m2_posterior$PosteriorMean))  %>% arrange(variant_gene) %>% mutate(Timepoint = "T0", Method = "Original")
T1_beta <- T1 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>% filter(variant_gene %in% rownames(m2_posterior$PosteriorMean))  %>% arrange(variant_gene) %>% mutate(Timepoint = "T1", Method = "Original")
T2_beta <- T2 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>% filter(variant_gene %in% rownames(m2_posterior$PosteriorMean))  %>% arrange(variant_gene) %>% mutate(Timepoint = "T2", Method = "Original")
T3_beta <- T3 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>% filter(variant_gene %in% rownames(m2_posterior$PosteriorMean))  %>% arrange(variant_gene) %>% mutate(Timepoint = "T3", Method = "Original")
T4_beta <- T4 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope)%>% filter(variant_gene %in% rownames(m2_posterior$PosteriorMean))  %>% arrange(variant_gene) %>% mutate(Timepoint = "T4", Method = "Original")


T0_se <- T0 %>% select(1,2,9) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope_se) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T0", Method = "Original")

T1_se <- T1 %>% select(1,2,9) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope_se) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T1", Method = "Original")

T2_se <- T2 %>% select(1,2,9) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope_se) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T2", Method = "Original")

T3_se <- T3 %>% select(1,2,9) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope_se) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T3", Method = "Original")

T4_se <- T4 %>% select(1,2,9) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope_se) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T4", Method = "Original")


new_se_T0 <- updated_sds[,1] %>% as.data.frame() %>% mutate(Timepoint = "T0", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_se_T0)[2] <- "slope_se"
new_se_T1 <- updated_sds[,2] %>% as.data.frame() %>% mutate(Timepoint = "T1", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_se_T1)[2] <- "slope_se"

new_se_T2 <- updated_sds[,3] %>% as.data.frame() %>% mutate(Timepoint = "T2", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_se_T2)[2] <- "slope_se"

new_se_T3 <- updated_sds[,4] %>% as.data.frame() %>% mutate(Timepoint = "T3", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_se_T3)[2] <- "slope_se"

new_se_T4 <- updated_sds[,5] %>% as.data.frame() %>% mutate(Timepoint = "T4", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_se_T4)[2] <- "slope_se"

all(colnames(new_se_T4) == colnames(T4_se))
all(rownames(new_se_T4) == T4_se$variant_gene)



T0_effect <- T0 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T0", Method = "Original")

T1_effect <- T1 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T1", Method = "Original")

T2_effect <- T2 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T2", Method = "Original")

T3_effect <- T3 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T3", Method = "Original")

T4_effect <- T4 %>% select(1,2,8) %>% mutate(variant_gene = paste0(variant_id,"-",phenotype_id)) %>% select(variant_gene, slope) %>%
  filter(variant_gene %in% rownames(m2_posterior$PosteriorMean)) %>% arrange(variant_gene) %>% mutate(Timepoint = "T4", Method = "Original")


new_effect_T0 <- updated_means[,1] %>% as.data.frame() %>% mutate(Timepoint = "T0", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_effect_T0)[2] <- "slope"
new_effect_T1 <- updated_means[,2] %>% as.data.frame() %>% mutate(Timepoint = "T1", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_effect_T1)[2] <- "slope"

new_effect_T2 <- updated_means[,3] %>% as.data.frame() %>% mutate(Timepoint = "T2", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_effect_T2)[2] <- "slope"

new_effect_T3 <- updated_means[,4] %>% as.data.frame() %>% mutate(Timepoint = "T3", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_effect_T3)[2] <- "slope"

new_effect_T4 <- updated_means[,5] %>% as.data.frame() %>% mutate(Timepoint = "T4", Method = "MASHR", variant_gene = rownames(updated_sds)) %>% 
  select(4,1,2,3) %>% arrange(variant_gene)
colnames(new_effect_T4)[2] <- "slope"

all(colnames(new_effect_T4) == colnames(T4_se))
all(rownames(new_effect_T4) == T4_se$variant_gene)

colnames(T0_se)
merged_se <- rbind(T0_se, T1_se,T2_se,T3_se,T4_se, new_se_T0, new_se_T1, new_se_T2, new_se_T3, new_se_T4)
merged_se$Metric = "TensorQTL Standard Error / MASHR Posterior Standard Deviation"
merged_beta <- rbind(T0_beta, T1_beta,T2_beta,T3_beta,T4_beta, new_effect_T0, new_effect_T1, new_effect_T2, new_effect_T3, new_effect_T4)
merged_beta$Metric = "TensorQTL slope / MASHR Posterior Mean"

colnames(merged_beta) <- colnames(merged_se)

mashr_visualisation <- rbind(merged_se, merged_beta)


#######################################
# Check to see variance and median effect of both distributions before and after applying MASHR
#######################################

head(merged_beta)

Mahsh_r_plot <- merged_beta %>% ggplot(aes(x = Timepoint, y = slope_se)) + geom_boxplot(aes(fill = factor(Method, levels = c("Original", "MASHR")), )) + labs(y = "TensorQTL effect size / mashR posterior mean effect size", x = "Timepoint", fill = "Method") + theme_bw() + scale_y_continuous(limits = c(-3,3), breaks = c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5)) +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", face = "bold"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15)) + scale_fill_manual(values = c("#762a83", "#1b7837")) # Stat compare means doesn't do F test, see below for P values
Mahsh_r_plot



Mahsh_r_plot <- merged_se %>% ggplot(aes(x = Timepoint, y = slope_se)) + geom_boxplot(aes(fill = factor(Method, levels = c("Original", "MASHR")), )) + labs(y = "TensorQTL effect size / mashR posterior mean effect size", x = "Timepoint", fill = "Method") + theme_bw() + scale_y_continuous(limits = c(0,1.25), breaks = c(0,0.25,0.5,0.75,1,1.25)) +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", face = "bold"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15)) + scale_fill_manual(values = c("#762a83", "#1b7837")) + stat_compare_means(aes(group = Method), method = "t.test", paired = TRUE)
Mahsh_r_plot


t.test(T0_beta$slope, new_effect_T0$slope, paired = TRUE)


t.test(T0_se$slope_se, new_se_T0$slope_se, paired = TRUE, alternative = "two.sided") 
var.test(T0_beta$slope, new_effect_T0$slope, paired = TRUE, alternative = "two.sided") 


t.test(T1_se$slope_se, new_se_T1$slope_se, paired = TRUE, alternative = "two.sided") 
var.test(T1_beta$slope, new_effect_T1$slope, paired = TRUE, alternative = "two.sided") 


t.test(T2_se$slope_se, new_se_T2$slope_se, paired = TRUE, alternative = "two.sided") 
var.test(T2_beta$slope, new_effect_T2$slope, paired = TRUE, alternative = "two.sided") 


t.test(T3_se$slope_se, new_se_T3$slope_se, paired = TRUE, alternative = "two.sided")
var.test(T3_beta$slope, new_effect_T3$slope, paired = TRUE, alternative = "two.sided") 


t.test(T4_se$slope_se, new_se_T4$slope_se, paired = TRUE, alternative = "two.sided")
var.test(T4_beta$slope, new_effect_T4$slope, paired = TRUE, alternative = "two.sided")




new_effect_T0 <- new_effect_T0 %>% select(1,2)
colnames(new_effect_T0)[2] <- "T0_posterior_mean"
new_effect_T1 <- new_effect_T1 %>% select(1,2)
colnames(new_effect_T1)[2] <- "T1_posterior_mean"
new_effect_T2 <- new_effect_T2 %>% select(1,2)
colnames(new_effect_T2)[2] <- "T2_posterior_mean"
new_effect_T3 <- new_effect_T3 %>% select(1,2)
colnames(new_effect_T3)[2] <- "T3_posterior_mean"
new_effect_T4 <- new_effect_T4 %>% select(1,2)
colnames(new_effect_T4)[2] <- "T4_posterior_mean"
all(rownames(new_effect_T0) == rownames(new_effect_T4))

new_effect_merged <- left_join(new_effect_T0, new_effect_T1)
new_effect_merged <- left_join(new_effect_merged, new_effect_T2)
new_effect_merged <- left_join(new_effect_merged, new_effect_T3)
new_effect_merged <- left_join(new_effect_merged, new_effect_T4)

new_se_T0 <- T0_se %>% select(1,2)
colnames(new_se_T0)[2] <- "T0_posterior_sd"
new_se_T1 <- T1_se %>% select(1,2)
colnames(new_se_T1)[2] <- "T1_posterior_sd"
new_se_T2 <- T2_se %>% select(1,2)
colnames(new_se_T2)[2] <- "T2_posterior_sd"
new_se_T3 <- T3_se %>% select(1,2)
colnames(new_se_T3)[2] <- "T3_posterior_sd"
new_se_T4 <- T4_se %>% select(1,2)
colnames(new_se_T4)[2] <- "T4_posterior_sd"
all(rownames(new_se_T0) == rownames(new_se_T4))

new_se_merged <- left_join(new_se_T0, new_se_T1)
new_se_merged <- left_join(new_se_merged, new_se_T2)
new_se_merged <- left_join(new_se_merged, new_se_T3)
new_se_merged <- left_join(new_se_merged, new_se_T4)



mas_signif_table$variant_gene <- paste0(mas_signif_table$variant, "-", mas_signif_table$gene)
mas_signif_table <- left_join(mas_signif_table, new_effect_merged)
mas_signif_table <- left_join(mas_signif_table, new_se_merged)

t.test(T0$slope_se, new_se_T0$T0_posterior_sd,)

#write.table(mas_signif_table, file = args[8], sep = "\t", quote = FALSE)


 
bhat_df <- cbind(T0$slope, T1$slope, T2$slope, T3$slope, T4$slope)

# Now with correlations
het.norm=function(effectsize){ t(apply(effectsize,1,function(x){x/x[which.max(abs(x))]}))}
pm.mash <- m2_posterior$PosteriorMean # check this
pm.mash.beta <- pm.mash 
pm.mash.beta.sig <- pm.mash.beta[nsig > 0,]
pm.mash.beta.norm <- het.norm(effectsize = pm.mash.beta)
pm.mash.beta.norm <- pm.mash.beta.norm[nsig > 0,]
colnames(pm.mash.beta.norm)





mag_by_tis <- function(t1, t2) {

  if(t2 == "T1") {
    # Filter column for both of interest and merge
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|1", colnames(lfsr.mash.sig))]
    # Count number of significant and not significant
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    nonsig_tis <- rowSums(lfsr.mash_tis > nonsig_thresh)
    # Filter and merge for significant ones
    pm.mash.beta_tis_union <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|1", colnames(pm.mash.beta.sig))]
    dim(pm.mash.beta_tis_union) # significant in at least one tissue

    # Filter for those significant and not significant
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0&nonsig_tis==0, grepl("T0|1", colnames(pm.mash.beta.sig))]
  }
  else if (t2 == "T2") {
        # Filter column for both of interest and merge
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|2", colnames(lfsr.mash.sig))]
    # Count number of significant and not significant
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    nonsig_tis <- rowSums(lfsr.mash_tis > nonsig_thresh)
    # Filter and merge for significant ones
    pm.mash.beta_tis_union <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|2", colnames(pm.mash.beta.sig))]
    
    # Filter for those significant and not significant
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0&nonsig_tis==0, grepl("T0|2", colnames(pm.mash.beta.sig))]
  }
    else if (t2 == "T3") {
        # Filter column for both of interest and merge
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|3", colnames(lfsr.mash.sig))]
    # Count number of significant and not significant
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    nonsig_tis <- rowSums(lfsr.mash_tis > nonsig_thresh)
    # Filter and merge for significant ones
    pm.mash.beta_tis_union <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|3", colnames(pm.mash.beta.sig))]
    
    # Filter for those significant and not significant
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0&nonsig_tis==0, grepl("T0|3", colnames(pm.mash.beta.sig))]
  }
     else if (t2 == "T4") {
        # Filter column for both of interest and merge
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|4", colnames(lfsr.mash.sig))]
    # Count number of significant and not significant
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    nonsig_tis <- rowSums(lfsr.mash_tis > nonsig_thresh)
    # Filter and merge for significant ones
    pm.mash.beta_tis_union <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|4", colnames(pm.mash.beta.sig))]
    # Filter for those significant and not significant
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0&nonsig_tis==0, grepl("T0|4", colnames(pm.mash.beta.sig))]
  }

  # Rename and convert to numeric
  colnames(pm.mash.beta_tis) <- c("T0", "T1")
  # Calculate sharing and if number of eQTLs between both are ahared (within MAG)
  pm.mash.beta_tis <- pm.mash.beta_tis %>% data.frame() %>% mutate(sharing=ifelse(T1/T0 >(1/mag)&T1/T0<mag, TRUE, FALSE))
  dim(pm.mash.beta_tis_union)[1]-sum(pm.mash.beta_tis$sharing) #time-specific eGenes
}

mag=2
thresh=0.05 
nonsig_thresh=1
time_mag_t0_t1 <- mag_by_tis(args[1], args[2])
time_mag_t0_t2 <- mag_by_tis(args[1], args[3])
time_mag_t0_t3 <- mag_by_tis(args[1], args[4])
time_mag_t0_t4 <- mag_by_tis(args[1], args[5])

time_mag_t0_t1
time_mag_t0_t2 

time_mag_t0_t3 
time_mag_t0_t4 



257 + 217 + 236 + 303


#write dynamic eGenes  
#write dynamic eGenes

mag_by_tis_gene <- function(t0, t1) {

  if(t1 == "T1") {
      
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|1", colnames(lfsr.mash.sig))]
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    lfsr.mash_tis <- lfsr.mash_tis[nsig_tis>0,]
    colnames(lfsr.mash_tis) <- c("T0_lfsr", "T1_lfsr")
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|1", colnames(pm.mash.beta.sig))]
    colnames(pm.mash.beta_tis) <- c("T0_beta", "T1_beta")
  }
  else if (t1 == "T2") {
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|2", colnames(lfsr.mash.sig))]
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    lfsr.mash_tis <- lfsr.mash_tis[nsig_tis>0,]
    colnames(lfsr.mash_tis) <- c("T0_lfsr", "T1_lfsr")
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|2", colnames(pm.mash.beta.sig))]
    colnames(pm.mash.beta_tis) <- c("T0_beta", "T1_beta")
  }
  else if (t1 == "T3") {
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|3", colnames(lfsr.mash.sig))]
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    lfsr.mash_tis <- lfsr.mash_tis[nsig_tis>0,]
    colnames(lfsr.mash_tis) <- c("T0_lfsr", "T1_lfsr")
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|3", colnames(pm.mash.beta.sig))]
    colnames(pm.mash.beta_tis) <- c("T0_beta", "T1_beta")
  }
  else if (t1 == "T4") {
    lfsr.mash_tis <- lfsr.mash.sig[,grepl("T0|4", colnames(lfsr.mash.sig))]
    nsig_tis <- rowSums(lfsr.mash_tis < thresh)
    lfsr.mash_tis <- lfsr.mash_tis[nsig_tis>0,]
    colnames(lfsr.mash_tis) <- c("T0_lfsr", "T1_lfsr")
    pm.mash.beta_tis <- pm.mash.beta.sig[nsig_tis>0, grepl("T0|4", colnames(pm.mash.beta.sig))]
    colnames(pm.mash.beta_tis) <- c("T0_beta", "T1_beta")
  }
  pm.mash.beta_tis <- lfsr.mash_tis %>% cbind(pm.mash.beta_tis) %>% as.data.frame() %>% 
  mutate(sharing=ifelse((T0_lfsr<thresh | T1_lfsr<thresh) &
                          T0_lfsr<nonsig_thresh & T1_lfsr < nonsig_thresh & 
                          T1_beta/T0_beta>(1/mag) & T1_beta/T0_beta < mag, TRUE, FALSE),
  effect_abs_diff=abs(T1_beta-T0_beta), mag_diff=abs(T1_beta/T0_beta), gene=as.character(lapply(strsplit(rownames(.), '-'), `[[`, 2)), gene_snp=rownames(.)) %>% 
  filter(!sharing) %>% arrange(desc(effect_abs_diff)) %>% dplyr::select(gene, gene_snp, T0_lfsr, T1_lfsr, T0_beta, T1_beta, effect_abs_diff, mag_diff, sharing)
}
0.00001/0.81

abs(0.81) / (0.0001)

time_mag_t0_t1_gene <- mag_by_tis_gene(args[1], args[2])
time_mag_t0_t2_gene <- mag_by_tis_gene(args[1], args[3]) 
time_mag_t0_t3_gene <- mag_by_tis_gene(args[1], args[4]) 
time_mag_t0_t4_gene <- mag_by_tis_gene(args[1], args[5])
dim(time_mag_t0_t4_gene)


dim(time_mag_t0_t1_gene)




all_genes <- rbind(time_mag_t0_t1_gene, time_mag_t0_t2_gene, time_mag_t0_t3_gene, time_mag_t0_t4_gene)

all_genes <- all_genes %>% arrange(desc(mag_diff))
all_genes
all_genes <- all_genes %>% filter(!duplicated(gene))
dim(all_genes) 
all_genes$gene
colnames(time_mag_t0_t2_gene) <- c("gene", "gene_snp", "T0_lfsr", "T2_lfsr", "T0_beta", "T2_beta", "effect_abs_diff", "mag_diff", "sharing")
colnames(time_mag_t0_t3_gene) <- c("gene", "gene_snp", "T0_lfsr", "T3_lfsr", "T0_beta", "T3_beta", "effect_abs_diff", "mag_diff", "sharing")
colnames(time_mag_t0_t4_gene) <- c("gene", "gene_snp", "T0_lfsr", "T4_lfsr", "T0_beta", "T4_beta", "effect_abs_diff", "mag_diff", "sharing")

x <- list(T0_T1 = time_mag_t0_t1_gene$gene_snp, T0_T2 = time_mag_t0_t2_gene$gene_snp, T0_T3 = time_mag_t0_t3_gene$gene_snp, T0_T4 = time_mag_t0_t4_gene$gene_snp)

DR_eQTLs <- unique(c(time_mag_t0_t1_gene$gene_snp, T0_T2 = time_mag_t0_t2_gene$gene_snp, T0_T3 = time_mag_t0_t3_gene$gene_snp, T0_T4 = time_mag_t0_t4_gene$gene_snp))

overlap_venn <- 
UpSetR::fromList(list(T0_T1 = time_mag_t0_t1_gene$gene_snp, T0_T2 = time_mag_t0_t2_gene$gene_snp, T0_T3 = time_mag_t0_t3_gene$gene_snp, T0_T4 = time_mag_t0_t4_gene$gene_snp))
write.table(overlap_venn, file = '/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_03B_venn.txt', sep = "\t", quote = FALSE)

length(DR_eQTLs)

# Get in annotation
# get all information in one file
symbols <- fread("/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf") %>% filter(V3 == "gene")
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)
symbols$gene_name <- gsub(" ", "", symbols$gene_name)
symbols <- symbols %>% select(gene_id, gene_name)

head(symbols)
time_mag_t0_t1_gene <- left_join(time_mag_t0_t1_gene, symbols, by = c("gene" = "gene_id"))
time_mag_t0_t2_gene <- left_join(time_mag_t0_t2_gene, symbols, by = c("gene" = "gene_id"))
time_mag_t0_t3_gene <- left_join(time_mag_t0_t3_gene, symbols, by = c("gene" = "gene_id"))
time_mag_t0_t4_gene <- left_join(time_mag_t0_t4_gene, symbols, by = c("gene" = "gene_id"))

time_mag_t0_t1_gene$snp <-  as.character(lapply(strsplit(time_mag_t0_t1_gene$gene_snp, '-'), `[[`, 1))
time_mag_t0_t2_gene$snp <-  as.character(lapply(strsplit(time_mag_t0_t2_gene$gene_snp, '-'), `[[`, 1))
time_mag_t0_t3_gene$snp <-  as.character(lapply(strsplit(time_mag_t0_t3_gene$gene_snp, '-'), `[[`, 1))
time_mag_t0_t4_gene$snp <-  as.character(lapply(strsplit(time_mag_t0_t4_gene$gene_snp, '-'), `[[`, 1))
colnames(time_mag_t0_t1_gene)
colnames(time_mag_t0_t4_gene)


time_mag_t0_t1_gene <- time_mag_t0_t1_gene %>% select(2,11,1,10,3,4,5,6)
time_mag_t0_t2_gene <- time_mag_t0_t2_gene %>% select(2,11,1,10,3,4,5,6)
time_mag_t0_t3_gene <- time_mag_t0_t3_gene %>% select(2,11,1,10,3,4,5,6)
time_mag_t0_t4_gene <- time_mag_t0_t4_gene %>% select(2,11,1,10,3,4,5,6)
colnames(time_mag_t0_t4_gene)
dim(time_mag_t0_t1_gene)
#all_genes$gene <- gsub("\\..*", "", all_genes$gene)


#write.table(time_mag_t0_t1_gene, file = args[11], sep = "\t")
#write.table(time_mag_t0_t2_gene, file = args[12], sep = "\t")
#write.table(time_mag_t0_t3_gene, file = args[13], sep = "\t")
#write.table(time_mag_t0_t4_gene, file = args[14], sep = "\t")






intersect_reQTLs <- intersect(time_mag_t0_t1_gene$gene_snp, intersect(time_mag_t0_t2_gene$gene_snp, intersect(time_mag_t0_t3_gene$gene_snp, time_mag_t0_t4_gene$gene_snp)))
length(intersect_reQTLs)
#all_genes <- all_genes[intersect_reQTLs,]
all_genes <- all_genes %>% arrange(desc(mag_diff)) %>% select(gene) %>% as.vector()
t0_t1_common <- time_mag_t0_t1_gene[time_mag_t0_t1_gene$gene_snp %in% intersect_reQTLs,]
t0_t2_common <- time_mag_t0_t2_gene[time_mag_t0_t2_gene$gene_snp %in% intersect_reQTLs,]
t0_t3_common <- time_mag_t0_t3_gene[time_mag_t0_t3_gene$gene_snp %in% intersect_reQTLs,]
t0_t4_common <- time_mag_t0_t4_gene[time_mag_t0_t4_gene$gene_snp %in% intersect_reQTLs,]
colnames(time_mag_t0_t1_gene)

colnames(time_mag_t0_t2_gene)

colnames(time_mag_t0_t3_gene)
colnames(time_mag_t0_t4_gene)


ALL_common <- cbind(t0_t1_common, t0_t2_common[, c(6,8)], t0_t3_common[,c(6,8)], t0_t4_common[,c(6,8)])
ALL_common_betas <- ALL_common %>% select(grep("gene_name|T*_beta", colnames(ALL_common)))
ALL_common_betas <- cbind(ALL_common_betas, as.vector(t0_t1_common$gene_snp))


colnames(ALL_common_betas)[7] <- "Gene_SNP"
ALL_common_betas <- ALL_common_betas %>% mutate(Treated_effect_direction = case_when(
  T1_beta > 0 & T2_beta > 0 & T3_beta > 0 & T4_beta > 0 ~ "Consistently Positive",
  T1_beta < 0 & T2_beta < 0 & T3_beta < 0 & T4_beta < 0 ~ "Consistently Negative",
 .default = "Spurious")) %>% mutate(line_alpha = if_else(Treated_effect_direction == "Consistently Positive" | Treated_effect_direction == "Consistently Negative", 1, 0.2))


ALL_common_betas_long <- pivot_longer(ALL_common_betas, cols = c(colnames(ALL_common_betas)[2:6]), names_to = "Timepoint", values_to = c("Betas"))

ALL_common_betas_long <- ALL_common_betas_long %>% mutate(label = if_else(Timepoint == "T4_beta" & (Treated_effect_direction == "Consistently Positive" | Treated_effect_direction == "Consistently Negative"), gene_name, NA))

my_palette = c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026")

ggplot(ALL_common_betas_long, aes(x = Timepoint, y = Betas, group = Gene_SNP, label = label)) + 
geom_line(data = ALL_common_betas_long, inherit.aes = T, aes(col = Treated_effect_direction, alpha = line_alpha)) +  
geom_point(size = 3, aes(fill = Timepoint), shape = 21) + theme_bw() +
scale_fill_manual(values = my_palette) +
scale_color_manual(values = c("steelblue",  "purple", "black")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size= 12, colour = "black"),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 10)) +
        labs(alpha = "") +
        geom_text_repel(max.overlaps = 40) 

length(unique(ALL_common_betas_long$gene_name))
length(unique(all_genes))




ALL_common_betas_long$Gene_id <- ALL_common_betas_long$Gene_SNP
ALL_common_betas_long$Gene_id <- sub(".*-", "", ALL_common_betas_long$Gene_id)
ALL_common_betas_long$Gene_id <- sub("\\..*", "", ALL_common_betas_long$Gene_id)
write.table(common_genes, args[15], quote = F, sep = "\t", col.names = F, row.names = F)
write.table(unique(ALL_common_betas_long$Gene_id), args[16], quote = F, sep = "\t", col.names = F, row.names = F)
unique(ALL_common_betas_long$Gene_id)



all_genes$gene <- gsub("\\..*", "", all_genes$gene)
common_genes <- gsub("\\..*", "", common_genes)

head(common_genes)

gostres <- gost(query =all_genes$gene,
ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "fdr", 
                custom_bg = common_genes)
df <- as.data.frame(do.call(cbind, gostres$result))
df <- apply(df,2,as.character)
write.table(df, file = args[17], sep = "\t", quote = FALSE)

terms <- c("TNF-alpha/NF-kappa B signaling complex 10",
           "HSP90-CDC37 chaperone complex",
           "Aryl hydrocarbon receptor pathway",
           "Cellular response to chemical stress",
           "regulation of chromatin binding",
           "NOD-like receptor signaling pathway",
           "RIPK1-mediated regulated necrosis",
           "Negative regulation of NOTCH4 signaling",
           "Drug-mediated inhibition of ERBB2 signaling",
           "HIF1A and PPARG regulation of glycolysis",
           "MLL4 complex",
           "PTIP-HMT complex",
           "Interleukin-1 signaling",
           "Proteasome",
           "Hsp90 protein binding"

)



gostres$result <- gostres$result %>%
  mutate(label = ifelse(term_name %in% terms, term_name, NA))
#colnames(result_lables)[11] <- "term_name"
result_df <- gostres$result
    result_df <- result_df %>% mutate(alpha_value = ifelse(term_name %in% terms, 1, 0.5))
    color_palette <- viridis(5)
    pos <- position_jitter(width = 0.15, seed = 3)
    my_plot <- ggplot(result_df, aes(x = source, y = -log10(p_value))) +
    geom_jitter(aes(color = source), alpha = result_df$alpha_value, position = pos, size = 3) +  # Use shape instead of color
    scale_color_brewer(palette = "Dark2" , name = "Source") +  # Apply the color palette
    scale_y_continuous(limits = c(0, 3), breaks = seq(0:3)) +
    labs(y = expression(-log[10](italic(P)[adj])), x = "") +
    theme_bw() +
    labs(size = "Intersection\nSize") +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", face = "bold"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15)) +
    geom_label_repel(aes(x=source, y=-log10(p_value)), label = gostres$result$label, max.overlaps = 20,
                   size=3.0, color='black', fontface = "bold",
                   fill='#FFFFFF33',
                   position = pos,
                   box.padding = 1,
                   point.padding = 0,
                   force = 4
    ) +
    geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
    guides(color = guide_legend(override.aes = list(size = 5)))
my_plot
ggsave(args[18], width = 12, height = 12, dpi = 600)



Mahsh_r_plot
ggsave(args[19], width = 12, height = 12)


# Distance to TSS

tss_mash <- lfsr.mash.sig
tss_mash <- as.data.frame(tss_mash)
tss_mash$variant_gene <- rownames(tss_mash)
tss_mash <- tss_mash %>% pivot_longer(names_to = "Timepoint", values_to = "LFSR", cols = colnames(tss_mash)[1:5])


T0$variant_gene <- paste0(T0$variant_id, "-", T0$phenotype_id)
T0_distance <- T0 %>% select(start_distance, variant_gene)


DR_eQTLs <- as.data.frame(DR_eQTLs)
DR_eQTLs$Timepoint <- "Response-eQTL"
DR_eQTLs$LFSR <- 0

colnames(DR_eQTLs)[1] <- "variant_gene"
tss_mash <- rbind(tss_mash, DR_eQTLs)

tss_mash <- left_join(tss_mash, T0_distance)
table(is.na(tss_mash$start_distance))

tss_mash <- tss_mash %>% filter(LFSR < 0.05)
tss_mash$Timepoint
tss_mash$Timepoint <- factor(tss_mash$Timepoint, levels = c("T0", "T1", "T2", "T3", "T4","Response-eQTL"))

ggplot(data = tss_mash, aes(x = start_distance, fill = Timepoint, col = Timepoint)) +
geom_histogram(alpha= 0.5, position = 'identity', bins = 99) + 
geom_density(data = tss_mash, aes(x=start_distance, fill=Timepoint), alpha=.3) + 
scale_fill_manual(values = c(my_palette, "purple")) + theme_bw()#+ facet_wrap(~ Timepoint)

write.table(tss_mash, file = '/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_03D_TSS_distance.txt', sep = "\t", quote = FALSE, row.names = FALSE)
# 1. Create the histogram plot
phist <- gghistogram(
  tss_mash, x = "start_distance", rug = FALSE,
  fill = "Timepoint", palette = c(my_palette, "purple"), alpha = 0.3, bins = 49
)

# 2. Create the density plot with y-axis on the right
# Remove x axis elements
pdensity <- ggdensity(
tss_mash, x = "start_distance",
  col = "Timepoint", palette = c(my_palette, "purple"), alpha = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")

# 3. Align the two plots and then overlay them.
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
ggsave(args[20], width = 12, height = 12, dpi = 600)


tss_mash %>% group_by(Timepoint) %>% summarise(Mean = mean(start_distance))
tss_mash %>%
  group_by(Timepoint) %>%
  group_modify(~ {
    test_result <- shapiro.test(.x$start_distance)
    tibble(statistic = test_result$statistic, p.value = test_result$p.value)
  }) %>%
  ungroup()4
# A tibble: 6 × 3
#  Timepoint     statistic  p.value
#  <fct>             <dbl>    <dbl>
#1 T0                0.980 1.46e-23
#2 T1                0.980 2.13e-23
#3 T2                0.980 5.07e-23
#4 T3                0.980 4.72e-23
#5 T4                0.980 1.92e-23
#6 Response-eQTL     0.964 7.85e-11
 

test <- tss_mash %>% filter(Timepoint == "Response-eQTL")
shapiro.test(test$start_distance)




# check to see the number of mTB specific genes/ treated genes and dynamic eQTLs

T0_mtb_genes <- time_mag_t0_t1_gene %>% filter(T0_lfsr < 0.05 & T1_lfsr > 0.05) %>% select(gene_name)
T1_mtb_genes <- time_mag_t0_t2_gene %>% filter(T0_lfsr < 0.05 & T2_lfsr > 0.05) %>% select(gene_name)
T2_mtb_genes <- time_mag_t0_t3_gene %>% filter(T0_lfsr < 0.05 & T3_lfsr > 0.05) %>% select(gene_name)
T3_mtb_genes <- time_mag_t0_t4_gene %>% filter(T0_lfsr < 0.05 & T4_lfsr > 0.05) %>% select(gene_name)

mtb_specific <- intersect(T0_mtb_genes$gene_name, intersect(T1_mtb_genes$gene_name, intersect(T2_mtb_genes$gene_name, T3_mtb_genes$gene_name)))
mtb_specific




T0_treat_genes <- time_mag_t0_t1_gene %>% filter(T0_lfsr > 0.05 & T1_lfsr < 0.05) %>% select(gene_name)
T1_treat_genes <- time_mag_t0_t2_gene %>% filter(T0_lfsr > 0.05 & T2_lfsr < 0.05) %>% select(gene_name)
T2_treat_genes <- time_mag_t0_t3_gene %>% filter(T0_lfsr > 0.05 & T3_lfsr < 0.05) %>% select(gene_name)
T3_treat_genes <- time_mag_t0_t4_gene %>% filter(T0_lfsr > 0.05 & T4_lfsr < 0.05) %>% select(gene_name)


treat_specific <- intersect(T0_treat_genes$gene_name, intersect(T1_treat_genes$gene_name, intersect(T2_treat_genes$gene_name, T3_treat_genes$gene_name)))

##############################
# Plotting eQTLs
##############################

# Read in counts
counts0 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T0_residualised_expression.txt") %>% select(-c(1:4, 6)) #%>% as.matrix()
head(counts0)
rownames(counts0) <- counts0$gid

counts1 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T1_residualised_expression.txt") %>% select(-c(1:4, 6)) #%>% as.matrix()
rownames(counts1) <- counts1$gid



counts2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T2_residualised_expression.txt") %>% select(-c(1:4, 6)) #%>% as.matrix()
rownames(counts2) <- counts2$gid




counts3 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T3_residualised_expression.txt") %>% select(-c(1:4,6)) #%>% as.matrix()
rownames(counts3) <- counts3$gid




counts4 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T4_residualised_expression.txt") %>% select(-c(1:4,6)) #%>% as.matrix()
rownames(counts4) <- counts4$gid


# VCF data 
vcf <- vcfR::read.vcfR("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz")

vcf <- cbind(vcf@fix, vcf@gt)
head(vcf)
vcf <- as.data.frame(vcf)

vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("1\\|0:.*", paste0("HET"), x))
head(vcf)
vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("0\\|1:.*", paste0("HET"), x))
vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("0\\|0:.*", paste0("HOM_REF"), x))
vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("1\\|1:.*", paste0("HOM_ALT"), x))
head(vcf)
eQTL_plot <- function(gene_id, SNP_id, gene_name, counts0, counts1, counts2, counts3, counts4, vcf_file, HOM, HET, HOM_ALT) {

  counts0_gene <- counts0 %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 
  counts1_gene <- counts1 %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 
  counts2_gene <- counts2 %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 
  counts3_gene <- counts3 %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 
  counts4_gene <- counts4 %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 
  
  colnames(counts0_gene) <- "Expression"
  colnames(counts1_gene) <- "Expression"
  colnames(counts2_gene) <- "Expression"
  colnames(counts3_gene) <- "Expression"
  colnames(counts4_gene) <- "Expression"
  counts0_gene$Sample <- rownames(counts0_gene)
  counts1_gene$Sample <- rownames(counts1_gene)
  counts2_gene$Sample <- rownames(counts2_gene)
  counts3_gene$Sample <- rownames(counts3_gene)
  counts4_gene$Sample <- rownames(counts4_gene)
  
  counts0_gene$Time <- "T0"
  counts1_gene$Time <- "T1"
  counts2_gene$Time <- "T2"
  counts3_gene$Time <- "T3"
  counts4_gene$Time <- "T4"
  
  vcf_temp <- vcf %>% filter(ID == SNP_id)
  vcf_temp <- pivot_longer(vcf_temp, cols = c(10:57), names_to = "Sample", values_to = "Genotype")
  vcf_temp <- vcf_temp %>% dplyr::select(ID, Sample, Genotype)

  counts0_gene <- left_join(counts0_gene, vcf_temp, by = c("Sample" = "Sample"))
  counts1_gene <- left_join(counts1_gene, vcf_temp, by = c("Sample" = "Sample"))
  counts2_gene <- left_join(counts2_gene, vcf_temp, by = c("Sample" = "Sample"))
  counts3_gene <- left_join(counts3_gene, vcf_temp, by = c("Sample" = "Sample"))
  counts4_gene <- left_join(counts4_gene, vcf_temp, by = c("Sample" = "Sample"))
  
  ALL <- rbind(counts0_gene, counts1_gene, counts2_gene, counts3_gene, counts4_gene)
  ALL[, c(1)] <- sapply(ALL[, c(1)], as.numeric)

  # Normalize genotype strings (remove FORMAT suffix) and map to canonical labels
  ALL$Genotype <- as.character(ALL$Genotype)
  ALL$Genotype <- sub(":.*", "", ALL$Genotype)
  ALL$Genotype <- trimws(ALL$Genotype)
  ALL$Genotype[ALL$Genotype %in% c("0|0", "0/0", "0")] <- "HOM_REF"
  ALL$Genotype[ALL$Genotype %in% c("0|1", "1|0", "0/1", "1/0", "HET")] <- "HET"
  ALL$Genotype[ALL$Genotype %in% c("1|1", "1/1", "1")] <- "HOM_ALT"

  # Set factor and drop unused levels
  ALL$Genotype <- factor(ALL$Genotype, levels = c("HOM_REF", "HET", "HOM_ALT"))
  ALL$Genotype <- droplevels(ALL$Genotype)
  print(ALL$Genotype)
  ALL$Time <- factor(ALL$Time)

  # Build label mapping for observed levels only
  labmap <- c(HOM_REF = HOM, HET = HET, HOM_ALT = HOM_ALT)
  obs <- levels(ALL$Genotype)
  labels_to_use <- labmap[obs]
  print(obs)
  print(labels_to_use)

  p <- ggplot(ALL, aes(y = Expression, x = Genotype)) + 
    geom_boxplot(outlier.colour = NA, aes(fill = Time), trim=FALSE) +
    geom_jitter(shape=16, colour = "darkgrey", position=position_jitter(0.2)) + 
    scale_fill_manual(values = my_palette) +
    facet_wrap(~ Time,nrow = 1) + theme_bw() + xlab(SNP_id) + ylab(paste0("Residualised expression of ", gene_name)) +
    scale_x_discrete(limits = obs, labels = labels_to_use)
  return(list(ALL,p))
}
# TNFRS10A - MTB specific - Activator of apoptosis
TNFRS10A <- eQTL_plot(SNP_id = "8:23224131:A:T", gene_id = "ENSG00000104689.10", gene_name = "TNFRS10A", counts0, counts1, counts2, counts3, counts4, vcf_file = vcf, HOM = "AA", HET = "AT", HOM_ALT = "TT")
TNFRS10A[[2]]
write.table(TNFRS10A[[1]], file = '/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_03E_TNFRS10A.txt', sep = "\t", quote = FALSE, row.names = FALSE)

ggsave(args[27], width = 12, height = 12, dpi = 600)


IFGNR2 <- eQTL_plot(SNP_id = "21:33360260:G:A", gene_id = "ENSG00000159128.16", gene_name = "IFNGR2",counts0, counts1, counts2, counts3, counts4, vcf_file = vcf, HOM = "GG", HET = "GA", HOM_ALT = "AA")

IFGNR2[[2]]
write.table(IFGNR2[[1]], file = '/home/workspace/jogrady/heQTL/results/Response_1/Source_Data/Fig_03F_IFNGR2.txt', sep = "\t", quote = FALSE, row.names = FALSE)

ggsave(args[28], width = 12, height = 12, dpi = 600)
