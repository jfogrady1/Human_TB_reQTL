library(data.table)
library(ComplexHeatmap)

T0_V_T1 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T0_V_T1.txt", sep = "\t"))
T0_V_T2 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T0_V_T2.txt", sep = "\t"))
T0_V_T3 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T0_V_T3.txt", sep = "\t"))
T0_V_T4 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T0_V_T4.txt", sep = "\t"))



T1_V_T0 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T1_V_T0.txt", sep = "\t"))
T1_V_T2 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T1_V_T2.txt", sep = "\t"))
T1_V_T3 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T1_V_T3.txt", sep = "\t"))
T1_V_T4 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T1_V_T4.txt", sep = "\t"))

T2_V_T0 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T2_V_T0.txt", sep = "\t"))
T2_V_T1 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T2_V_T1.txt", sep = "\t"))
T2_V_T3 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T2_V_T3.txt", sep = "\t"))
T2_V_T4 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T2_V_T4.txt", sep = "\t"))

T3_V_T0 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T3_V_T0.txt", sep = "\t"))
T3_V_T1 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T3_V_T1.txt", sep = "\t"))
T3_V_T2 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T3_V_T2.txt", sep = "\t"))
T3_V_T4 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T3_V_T4.txt", sep = "\t"))


T4_V_T0 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T4_V_T0.txt", sep = "\t"))
T4_V_T1 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T4_V_T1.txt", sep = "\t"))
T4_V_T2 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T4_V_T2.txt", sep = "\t"))
T4_V_T3 <- as.matrix(read.csv("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_Matrix_T4_V_T3.txt", sep = "\t"))






pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T1.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T0_V_T1, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T0_V_T1 > 0.85, "*", ""), nrow(T0_V_T1)))
dev.off()

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T2.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T0_V_T1, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T0_V_T2 > 0.85, "*", ""), nrow(T0_V_T2)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T3.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T0_V_T1, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T0_V_T3 > 0.85, "*", ""), nrow(T0_V_T3)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T4.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T0_V_T1, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T0_V_T4 > 0.85, "*", ""), nrow(T0_V_T4)))
dev.off()




# T1 

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T0.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T1_V_T0, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T1_V_T0 > 0.85, "*", ""), nrow(T1_V_T0)))
dev.off()

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T2.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T1_V_T2, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T1_V_T2 > 0.85, "*", ""), nrow(T1_V_T2)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T3.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T1_V_T3, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T1_V_T3 > 0.85, "*", ""), nrow(T1_V_T3)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T4.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T1_V_T4, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T1_V_T4 > 0.85, "*", ""), nrow(T1_V_T4)))
dev.off()




# T2

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T0.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T2_V_T0, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T2_V_T0 > 0.85, "*", ""), nrow(T2_V_T0)))
dev.off()

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T1.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T2_V_T1, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T2_V_T1 > 0.85, "*", ""), nrow(T2_V_T1)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T3.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T2_V_T3, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T2_V_T3 > 0.85, "*", ""), nrow(T2_V_T3)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T4.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T2_V_T4, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T2_V_T4 > 0.85, "*", ""), nrow(T2_V_T4)))
dev.off()


# T3

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T0.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T3_V_T0, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T3_V_T0 > 0.85, "*", ""), nrow(T3_V_T0)))
dev.off()

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T1.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T3_V_T1, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T3_V_T1 > 0.85, "*", ""), nrow(T3_V_T1)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T2.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T3_V_T2, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T3_V_T2 > 0.85, "*", ""), nrow(T3_V_T2)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T4.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T3_V_T4, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T3_V_T4 > 0.85, "*", ""), nrow(T3_V_T4)))
dev.off()


# T4

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T0.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T4_V_T0, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T4_V_T0 > 0.85, "*", ""), nrow(T4_V_T0)))
dev.off()

pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T1.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T4_V_T1, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T4_V_T1 > 0.85, "*", ""), nrow(T4_V_T1)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T2.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T4_V_T2, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T4_V_T2 > 0.85, "*", ""), nrow(T4_V_T2)))
dev.off()


pdf(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T3.pdf"),width=12,height=12)
ComplexHeatmap::pheatmap(T4_V_T3, cluster_rows=F, cluster_cols=F,legend_breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
                         legend_labels = c("0","0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1"), border_color = "black", 
                         na_col= "white", fontsize_col = 10, fontsize_row = 10, cellheight = NA, row_names_side = c("left"), 
                         name = "Genotype\nConcordance", display_numbers = matrix(ifelse(T4_V_T3 > 0.85, "*", ""), nrow(T4_V_T3)))
dev.off()

