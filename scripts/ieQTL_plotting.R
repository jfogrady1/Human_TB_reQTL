library(tidyverse)
library(data.table)
library(arrow)
library(ash)
library(ashr)
library(mashr)
library(gtable)
library(cowplot)
library(viridis)


load('/home/workspace/jogrady/heQTL/results/ieQTLs/ieQTL_MASHR.RData')


# VCF data 
library(vcfR)
vcf <- vcfR::read.vcfR("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz")
vcf@gt
vcf <- cbind(vcf@fix, vcf@gt)
head(vcf)
vcf <- as.data.frame(vcf)

# read in expression data
T0_TFR2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T0_residualised_expression.txt") %>% filter(gid == "ENSG00000106327.13") %>% select(5,7:54)
T1_TFR2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T1_residualised_expression.txt") %>% filter(gid == "ENSG00000106327.13") %>% select(5,7:54)
T2_TFR2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T2_residualised_expression.txt") %>% filter(gid == "ENSG00000106327.13") %>% select(5,7:54)
T3_TFR2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T3_residualised_expression.txt") %>% filter(gid == "ENSG00000106327.13") %>% select(5,7:54)
T4_TFR2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T4_residualised_expression.txt") %>% filter(gid == "ENSG00000106327.13") %>% select(5,7:54)

head(T0_TFR2)
T0_TFR2 <- pivot_longer(T0_TFR2, names_to = "Mixture", values_to = "Expression", cols = 2:49) %>% select(2,3)
T1_TFR2 <- pivot_longer(T1_TFR2, names_to = "Mixture", values_to = "Expression", cols = 2:49) %>% select(2,3)
T2_TFR2 <- pivot_longer(T2_TFR2, names_to = "Mixture", values_to = "Expression", cols = 2:49) %>% select(2,3)
T3_TFR2 <- pivot_longer(T3_TFR2, names_to = "Mixture", values_to = "Expression", cols = 2:49) %>% select(2,3)
T4_TFR2 <- pivot_longer(T4_TFR2, names_to = "Mixture", values_to = "Expression", cols = 2:49) %>% select(2,3)
head(T2_TFR2)

# Read in interactions cell type proportion
T0_cell = fread('/home/workspace/jogrady/heQTL/work/ieQTL/T0_Classical_Monocyte_interaction_input_from_cibersort.txt')
T1_cell = fread('/home/workspace/jogrady/heQTL/work/ieQTL/T1_Classical_Monocyte_interaction_input_from_cibersort.txt')
T2_cell = fread('/home/workspace/jogrady/heQTL/work/ieQTL/T2_Classical_Monocyte_interaction_input_from_cibersort.txt')
T3_cell = fread('/home/workspace/jogrady/heQTL/work/ieQTL/T3_Classical_Monocyte_interaction_input_from_cibersort.txt')
T4_cell = fread('/home/workspace/jogrady/heQTL/work/ieQTL/T4_Classical_Monocyte_interaction_input_from_cibersort.txt')

# Read in covariates
T0_cov = fread("/home/workspace/jogrady/heQTL/data/covariate/T0_eqtlcovs.txt") %>% t()
T1_cov = fread("/home/workspace/jogrady/heQTL/data/covariate/T1_eqtlcovs.txt") %>% t()
T2_cov = fread("/home/workspace/jogrady/heQTL/data/covariate/T2_eqtlcovs.txt") %>% t()
T3_cov = fread("/home/workspace/jogrady/heQTL/data/covariate/T3_eqtlcovs.txt") %>% t()
T4_cov = fread("/home/workspace/jogrady/heQTL/data/covariate/T4_eqtlcovs.txt") %>% t()

colnames(T0_cov) <- T0_cov[1,]
colnames(T1_cov) <- T1_cov[1,]
colnames(T2_cov) <- T2_cov[1,]
colnames(T3_cov) <- T3_cov[1,]
colnames(T4_cov) <- T4_cov[1,]

T0_cov <- T0_cov[-1,]
T1_cov <- T1_cov[-1,]
T2_cov <- T2_cov[-1,]
T3_cov <- T3_cov[-1,]
T4_cov <- T4_cov[-1,]

labels <- rownames(T0_cov)

T0_cov <- apply(T0_cov, 2, as.numeric)
rownames(T0_cov) <- labels
T1_cov <- apply(T1_cov, 2, as.numeric)
rownames(T1_cov) <- labels
T2_cov <- apply(T2_cov, 2, as.numeric)
rownames(T2_cov) <- labels
T3_cov <- apply(T3_cov, 2, as.numeric)
rownames(T3_cov) <- labels
T4_cov <- apply(T4_cov, 2, as.numeric)
rownames(T4_cov) <- labels

all(rownames(T0_cov) == T0_cell$Mixture)
# Get everything in order for regression

vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf[,colnames(vcf)[10:57]] <- lapply(vcf[,colnames(vcf)[10:57]], function(x) sub("1\\|1:.*", paste0("2"), x))
vcf_filter <- vcf %>% filter(ID == "7:100615131:C:T")
head(vcf_filter)

vcf_filter <- vcf_filter %>% select(c(3,10:57))
vcf_filter <- vcf_filter %>% pivot_longer(names_to = "Mixture", values_to = "Genotype", cols = 2:49) %>% select(2,3)



ieQTL_plot <- function(gene, variant, gene_name, hom_ref, het, hom_alt, cell_type, groups = c("T0", "T1", "T2", "T3", "T4")) {
  # SERPINB9 gene is interesting
  # Variant = 6:25066467:G:T
  counts0 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T0_residualised_expression.txt") %>% select(5,7:54) %>% as.matrix()
  rownames(counts0) <- counts0[,1]
  counts0 <- as.matrix(counts0[,-1])



  counts1 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T1_residualised_expression.txt") %>% select(5,7:54) %>% as.matrix()
  rownames(counts1) <- counts1[,1]
  counts1 <- as.matrix(counts1[,-1])



  counts2 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T2_residualised_expression.txt") %>% select(5,7:54) %>% as.matrix()
  rownames(counts2) <- counts2[,1]
  counts2 <- as.matrix(counts2[,-1])



  counts3 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T3_residualised_expression.txt") %>% select(5,7:54) %>% as.matrix()
  rownames(counts3) <- counts3[,1]
  counts3 <- as.matrix(counts3[,-1])



  counts4 <- fread("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/T4_residualised_expression.txt") %>% select(5,7:54) %>% as.matrix()
  rownames(counts4) <- counts4[,1]
  counts4 <- as.matrix(counts4[,-1])





  T0 = fread(paste0('/home/workspace/jogrady/heQTL/work/ieQTL/T0_', cell_type, '_interaction_input_from_cibersort.txt'))
  T1 = fread(paste0('/home/workspace/jogrady/heQTL/work/ieQTL/T1_', cell_type, '_interaction_input_from_cibersort.txt'))
  T2 = fread(paste0('/home/workspace/jogrady/heQTL/work/ieQTL/T2_', cell_type, '_interaction_input_from_cibersort.txt'))
  T3 = fread(paste0('/home/workspace/jogrady/heQTL/work/ieQTL/T3_', cell_type, '_interaction_input_from_cibersort.txt'))
  T4 = fread(paste0('/home/workspace/jogrady/heQTL/work/ieQTL/T4_', cell_type, '_interaction_input_from_cibersort.txt'))



  #ENSG00000111913.20 - SERPINB9
  counts0 <- counts0[gene,]
  counts1 <- counts1[gene,]
  counts2 <- counts2[gene,]
  counts3 <- counts3[gene,]
  counts4 <- counts4[gene,]

  counts0 <- as.data.frame(counts0)
  counts1 <- as.data.frame(counts1)
  counts2 <- as.data.frame(counts2)
  counts3 <- as.data.frame(counts3)
  counts4 <- as.data.frame(counts4)
  counts0$Patient <- rownames(counts0)
  counts1$Patient <- rownames(counts1)
  counts2$Patient <- rownames(counts2)
  counts3$Patient <- rownames(counts3)
  counts4$Patient <- rownames(counts4)

  counts0$Time <- "T0"
  counts1$Time <- "T1"
  counts2$Time <- "T2"
  counts3$Time <- "T3"
  counts4$Time <- "T4"


  vcf_test <- vcf %>% filter(ID == variant)

  vcf_test <- t(vcf_test)
  vcf_test <- vcf_test[c(3,10:57), ]
  vcf_test <- as.data.frame(vcf_test)
  
  vcf_test <- vcf_test %>% mutate(genotype = case_when(vcf_test == "0" ~ paste0(hom_ref, ":", hom_ref),
                                          vcf_test == "1" ~ het,
                                          vcf_test == "2" ~ hom_alt))

  if (gene_name == "IFNGR2") {
    vcf_test <- vcf_test %>% mutate(genotype = case_when(vcf_test == "0" ~ paste0(hom_ref, ":", hom_ref),
                                          vcf_test == "1" ~ paste0(het, "/", hom_alt),
                                          vcf_test == "2" ~ paste0(het, "/", hom_alt)))
                                          }
  vcf_test$Patient <- rownames(vcf_test)
  print(head(vcf_test))


  counts0 <- left_join(counts0, T0, by = c("Patient" = "Mixture"))
  counts1 <- left_join(counts1, T1, by = c("Patient" = "Mixture"))
  counts2 <- left_join(counts2, T2, by = c("Patient" = "Mixture"))
  counts3 <- left_join(counts3, T3, by = c("Patient" = "Mixture"))
  counts4 <- left_join(counts4, T4, by = c("Patient" = "Mixture"))



  colnames(counts0)[1] <- gene
  colnames(counts1)[1] <- gene
  colnames(counts2)[1] <- gene
  colnames(counts3)[1] <- gene
  colnames(counts4)[1] <- gene

  colnames(counts1)
  exp <- rbind(counts0, counts1, counts2, counts3, counts4)
  vcf_test[-1,]
  colnames(exp)[1] <- "Gene"
  colnames(exp)[4] <- "Cell_Proportion"
  exp <- left_join(exp, vcf_test)
  print(exp)
  exp <- exp %>% filter(Time %in% groups)
  p <- ggplot(exp, aes(x = Cell_Proportion, y = as.numeric(Gene), col = genotype)) + geom_point(size =3, alpha = 0.8) + theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~Time, nrow = 1) +
  theme_bw() + labs(y = paste0("Residualised ", gene_name, "expression"), x = paste0("Normalised ", cell_type ,"proportion"), colour = variant) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.x = element_text(angle = 0, size = 21, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=5))) +
    scale_colour_viridis(discrete = TRUE)

  if (length(groups) > 1) {
    print(TRUE)
    p + facet_wrap(~Time, nrow = 1)
  }
  return(p)
}

NOD2_CM <- ieQTL_plot(variant = "16:50652646:A:C", gene = "ENSG00000167207.15", gene_name = "NOD2", hom_ref = "A", het = "A:C", hom_alt = "C:C", cell_type = "Classical_Monocyte", groups = c("T0", "T1", "T2", "T3", "T4"))
NOD2_CM + theme_bw() + labs(y = "Normalised NOD2 expression", x = "Normalised Classical Monocyte proportion", colour = "Genotype:\n16:50652646:A:C") +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    axis.title.x = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=5))) +
    scale_colour_viridis(discrete = TRUE)

ggsave("")



IL7R_CM <-  ieQTL_plot(variant = "5:35875125:G:A", gene = "ENSG00000168685.15", gene_name = "IL7R", hom_ref = "G", het = "G:A", hom_alt = "A:A", cell_type = "Classical_Monocyte", groups = c("T0", "T1", "T2", "T3", "T4"))
IL7R_CM + theme_bw() + labs(y = "Normalised IL7R expression", x = "Normalised Classical Monocyte proportion", colour = "Genotype:\n5:35875125:G:A") +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    axis.title.x = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=5))) +
    scale_colour_viridis(discrete = TRUE)




TFR2
TFR2_Classical_Monocyte <- TFR2 + theme_bw() + labs(y = "Normalised NLRC4 expression", x = "Normalised Classical Monocyte proportion", colour = "Genotype:2:32236868:A:AAAAG") +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=5))) +
    scale_colour_viridis(discrete = TRUE)

TFR2_Classical_Monocyte
library(grid)

grid.draw(shift_legend(RIPOR2_NK_cell))
ggsave("/home/workspace/jogrady/heQTL/work/ieQTL/NLRC4_Classical_Monocyte_example_ieQTL.pdf", width = 12, height = 12)
