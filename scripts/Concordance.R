#' concoRdance: A function to determine genotype concordance between 2 files

library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(pheatmap)


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args1 <- args[1]
args2 <- args[2]
out1 <- args[3]
out2 <- args[4]
geno_1 <- read.table(args[1], header = T)
# Filter for common positions present
# Only considered called genotypes in the truth
geno_1[geno_1 == "./."] <- NA
geno_1 <- geno_1 %>% drop_na()
geno_1 <- geno_1 %>% unite("POS", CHROM, POS, sep = ":")

geno_2 <- read.table(args[2], header = T)

geno_2 <- geno_2 %>% unite("POS", CHROM, POS, sep = ":")
geno_2 <- geno_2[geno_2$POS %in% geno_1$POS,]
geno_1 <- geno_1[geno_1$POS %in% geno_2$POS,]
geno_2 <- geno_2[geno_2$POS %in% geno_1$POS,]


# Get sample names from each vcf entry
samples_1 <- colnames(geno_1 %>% select(2:(1+48)))
samples_2 <- colnames(geno_2 %>% select(2:(1+48)))


# Set up final matrix
# dimesions correspond to sample names
concordance_matrix <- matrix(0,nrow = 48, ncol = 48)
colnames(concordance_matrix) <- samples_2
rownames(concordance_matrix) <- samples_1


stats_df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(stats_df) <- c("geno_1_Sample", "geno_2_Sample", "total_right", "total_wrong", "total_na", "prop")


for (geno_1_Sample in samples_1) {
  # Get the reference sample
  geno_1_subset = geno_1 %>% select(all_of(c("POS",geno_1_Sample)))
  for(geno_2_Sample in samples_2) {
    geno_2_subset = geno_2 %>% select(all_of(c("POS",geno_2_Sample)))
    geno_2_subset[geno_2_subset == "./."] <- NA
    geno_2_subset <- geno_2_subset %>% drop_na()
    
    
    geno_1_small <- geno_1_subset[geno_1_subset$POS %in% geno_2_subset$POS,]
    geno_2_subset <- geno_2_subset %>% select(all_of(geno_2_Sample))
    
    geno_1_small <- geno_1_small %>% select(all_of(geno_1_Sample))
    
    total_right = 0
    
    total_wrong = 0
    
    for (row in 1:nrow(geno_2_subset)) {
      current_geno_2 = geno_2_subset[row,]
      print(total_right)
      current_geno_1 = geno_1_small[row,]
      
      # Build the if else ladder - ignore the missing genotypes
      # Remember the phasing
      if (current_geno_1 == "0/0" & current_geno_2 == "0/0") {
        total_right = (total_right + 1)
      } else if (current_geno_1 == "0/1" & current_geno_2 == "0/1" | current_geno_2 == "1/0") {
        total_right = (total_right + 1)
      } else if (current_geno_1 == "1/0" & current_geno_2 == "0/1" | current_geno_2 == "1/0") {
        total_right = (total_right + 1)
      } else if (current_geno_1 == "1/1" & current_geno_2 == "1/1") {
        total_right = (total_right + 1)
        
      } else {
        total_wrong <- (total_wrong + 1)
      }
    }
    prop = round(total_right / (total_right + total_wrong),3)
    concordance_matrix[geno_1_Sample, geno_2_Sample] <- concordance_matrix[geno_1_Sample,geno_2_Sample]+prop
  }
}

png(paste0("/cluster/work/pausch/jogrady/heQTL/work/samplecheck/Concordance_",out1,"_V_",out2,".png"),width=12,height=12,units="in",res=800)
ComplexHeatmap::pheatmap(concordance_matrix, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
            legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
            na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
            name = "Genotype\nConcordance", display_numbers = matrix(ifelse(concordance_matrix > 0.9, "*", ""), nrow(concordance_matrix)), main = paste("Genotype concordance of samples in ", out1, " V ", out2))
dev.off()


final_matrix = concordance_matrix                                   
write.table(x = final_matrix, file = paste0("/cluster/work/pausch/jogrady/heQTL/work/samplecheck/Concordance_Matrix_",out1,"_V_",out2,".txt"), sep = "\t")

