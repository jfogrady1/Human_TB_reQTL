library(data.table)
library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
mbvFindBestMatch <- function(mbv_df){
  res = dplyr::transmute(mbv_df, mbv_genotype_id = SampleID,
                         het_consistent_frac = n_het_consistent/n_het_covered,
                         hom_consistent_frac = n_hom_consistent/n_hom_covered)
  
  #Identify best het
  best_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  other_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() > 1)
  best_row = dplyr::mutate(best_het, het_min_dist = min(best_het$het_consistent_frac - other_het$het_consistent_frac),
                           hom_min_dist = min(best_het$hom_consistent_frac - other_het$hom_consistent_frac),
                           distance = sqrt(het_min_dist^2 + hom_min_dist^2))
  
  #Compare against best hom
  best_hom = dplyr::arrange(res, -hom_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  if(best_row$mbv_genotype_id != best_hom$mbv_genotype_id){
    best_row = dplyr::mutate(best_row, het_consistent_frac = as.numeric(NA), hom_consistent_frac = as.numeric(NA),
                             het_min_dist = as.numeric(NA), hom_min_dist = as.numeric(NA))
  }
  return(best_row)
}


string_vector = c('SRR12609394','SRR12609405','SRR12609414','SRR12609425','SRR12609434','SRR12609449','SRR12609475','SRR12609489','SRR12609495','SRR12609506',
'SRR12609516','SRR12609526','SRR12609535','SRR12609546','SRR12609574','SRR12609595','SRR12609608','SRR12609626','SRR12609641','SRR12609651','SRR12609676','SRR12609687',
'SRR12609696','SRR12609770','SRR12609780','SRR12609800','SRR12609808','SRR12609817','SRR12609855','SRR12609875','SRR12609885','SRR12609897','SRR12609910',
'SRR12609922','SRR12609934','SRR12609974','SRR12609982','SRR12609994','SRR12610003','SRR12610014','SRR12610032','SRR12610053','SRR12610063','SRR12610071',
'SRR12610092','SRR12610115','SRR12610123','SRR12610148')


string_vector_1 = c('SRR12609772','SRR12609781','SRR12609801','SRR12609809','SRR12609818','SRR12609395','SRR12609406','SRR12609415','SRR12609426','SRR12609435',
'SRR12609856','SRR12609450','SRR12609866','SRR12609476','SRR12609876','SRR12609886','SRR12609898','SRR12609911','SRR12609923','SRR12609935','SRR12609496','SRR12609507',
'SRR12609517','SRR12609527','SRR12609536','SRR12609547','SRR12609575','SRR12609975','SRR12609983','SRR12609995','SRR12610004','SRR12610015','SRR12610054','SRR12609596',
'SRR12609609','SRR12609627','SRR12609642','SRR12609677','SRR12610033','SRR12610064','SRR12610072','SRR12609652','SRR12610093',
'SRR12610116','SRR12610124','SRR12610149','SRR12609688','SRR12609697')

string_vector_2 <- c("SRR12609774", 
"SRR12609783", "SRR12609803", "SRR12609812", "SRR12609820", "SRR12609397", "SRR12609409", "SRR12609418", 
"SRR12609429", "SRR12609438", "SRR12609859", "SRR12609453", "SRR12609869", "SRR12609480", "SRR12609880", "SRR12609890", "SRR12609901", "SRR12609914", "SRR12609926", "SRR12609938", 
"SRR12609499", "SRR12609509", "SRR12609520", "SRR12609530", "SRR12609539", "SRR12609550", "SRR12609578", "SRR12609978", "SRR12609986", "SRR12609998", "SRR12610007", 
"SRR12610018", "SRR12610057", "SRR12609599", "SRR12609612", "SRR12609629", "SRR12609645", "SRR12609680", "SRR12610036", "SRR12610066", "SRR12610075", 
"SRR12609655", "SRR12610096", "SRR12610119", "SRR12610127", "SRR12610152", "SRR12609691","SRR12609700")

string_vector_3 <-c(
'SRR12609775','SRR12609784','SRR12609804','SRR12609813','SRR12609821','SRR12609399','SRR12609411','SRR12609420','SRR12609431','SRR12609439','SRR12609861','SRR12609454',
'SRR12609871','SRR12609482','SRR12609882','SRR12609892','SRR12609903','SRR12609916','SRR12609928','SRR12609940','SRR12609501','SRR12609511','SRR12609522','SRR12609531',
'SRR12609541','SRR12609551','SRR12609580','SRR12609980','SRR12609988','SRR12610000','SRR12610009','SRR12610020','SRR12609587','SRR12609601','SRR12609614',
'SRR12609634','SRR12609647','SRR12610058','SRR12610037','SRR12610068','SRR12610077','SRR12609656','SRR12610098','SRR12610121','SRR12610129','SRR12609682',
'SRR12609693','SRR12609702')

string_vector_4 <-c(
"SRR12609776", 
"SRR12609787", 
"SRR12609805", 
"SRR12609814", 
"SRR12609822", 
"SRR12609401", 
"SRR12609412", 
"SRR12609422", 
"SRR12609433", 
"SRR12609440", 
"SRR12609863", 
"SRR12609457", 
"SRR12609872", 
"SRR12609484", 
"SRR12609883", 
"SRR12609894", 
"SRR12609905", 
"SRR12609918", 
"SRR12609930", 
"SRR12609942", 
"SRR12609503", 
"SRR12609513", 
"SRR12609523", 
"SRR12609533", 
"SRR12609543", 
"SRR12609552", 
"SRR12609581", 
"SRR12609981", 
"SRR12609990", 
"SRR12610001", 
"SRR12610011", 
"SRR12610021", 
"SRR12609589", 
"SRR12609603", 
"SRR12609616", 
"SRR12609636", 
"SRR12609649", 
"SRR12610060", 
"SRR12610038", 
"SRR12610069", 
"SRR12610079", 
"SRR12609658", 
"SRR12610101", 
"SRR12610122", 
"SRR12610131", 
"SRR12609684", 
"SRR12609694", 
"SRR12609704")



ids_timepoint <- c(75, 76, 83, 84, 85, 92, 105, 99, 157, 175, 176, 201, 203, 221, 227, 268, 278, 294, 297, 367, 302, 432, 461, 1, 3, 13, 29, 36, 87, 117, 124, 127, 136, 137, 146, 254, 256, 257, 258, 265, 327, 266, 335, 348, 371, 389, 395, 400)
ids_timepoint_0 <- paste0("P_",ids_timepoint)


ids_timepoint_1 <- paste0("P_",sort(ids_timepoint)) # This is very important, whatever way I ordered this timepoint, the SRR IDs are based on the ordered Patient id

ids_timepoint_2 <- paste0("P_",sort(ids_timepoint))

ids_timepoint_3 <- paste0("P_",sort(ids_timepoint))


ids_timepoint_4 <- paste0("P_", sort(ids_timepoint))


final_df_match <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match <- rbind(final_df_match, data)
}


final_df_match_2 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_2) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_1) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_2 <- rbind(final_df_match_2, data)
}

final_df_match_3 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_3) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_2) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_3 <- rbind(final_df_match_3, data)
}

final_df_match_4 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_4) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_3) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_4 <- rbind(final_df_match_4, data)
}

final_df_match_5 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_5) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_4) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_5 <- rbind(final_df_match_5, data)
}




final_df_match$SRR_PT_ID <- ids_timepoint_0
final_df_match_2$SRR_PT_ID <- ids_timepoint_1
final_df_match_3$SRR_PT_ID <- ids_timepoint_2
final_df_match_4$SRR_PT_ID <- ids_timepoint_3
final_df_match_5$SRR_PT_ID <- ids_timepoint_4

all(final_df_match$mbv_genotype_id == final_df_match$SRR_PT_ID)
all(final_df_match_2$mbv_genotype_id == final_df_match_2$SRR_PT_ID)
all(final_df_match_3$mbv_genotype_id == final_df_match_3$SRR_PT_ID)
all(final_df_match_4$mbv_genotype_id == final_df_match_4$SRR_PT_ID)
all(final_df_match_5$mbv_genotype_id == final_df_match_5$SRR_PT_ID)
head(final_df_match)


data_final <- rbind(final_df_match,
                    final_df_match_2,
                    final_df_match_3,
                    final_df_match_4,
                    final_df_match_5)

data_final <- data_final %>% select(1,2,3,4,8)
data_final <- data_final %>% pivot_longer(, cols = -c("sample_cur", "SRR_PT_ID", "mbv_genotype_id"), names_to ="Metric", values_to = "Proportion")
data_final$Metric = factor(data_final$Metric, levels = c("hom_consistent_frac", "het_consistent_frac"), labels = c("Consistent Homozygous", "Consistent Heterozygous"))
data_final$Label <- if_else(data_final$Proportion < 0.65, data_final$sample_cur, NA)
MBV <- ggplot(data = data_final, aes(x = SRR_PT_ID, y = Proportion, colour = Metric, group = interaction(SRR_PT_ID, Metric))) + 
geom_boxplot(alpha = 0.1, width = 0.75, outlier.color = NA) + 
  geom_point(position = position_dodge(width=0.75)) + 
  geom_text(nudge_x = 0, nudge_y = 0.01,label = data_final$Label, col = "black") +
  theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 90, hjust = 0.5),
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 11),
        axis.text.y = element_text(colour = "black", size = 11)) 

MBV
ggsave("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/MBV_initial.pdf", width = 20, height = 5, dpi = 600)
data_final

final_df_match_2 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_2) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_1) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_2 <- rbind(final_df_match_2, data)
}

sample_cur = as.character("SRR12609781")
data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
head(data)
colnames(data)
data <- data %>% select(1,9,10)
data$mbv_genotype_id <- "P_3"
data$match <- "FALSE"
data$label <- NA
colnames(data)[1] <- "sample_cur"
data$distance <- NA
data$het_min_dist <- NA
data$hom_min_dist <- NA
colnames(data)
data <- data %>% select(1,4,2,3,8,9,7,5,6)
colnames(final_df_match)
colnames(final_df_match)
library(ggrepel)
head(data)
data <- data %>% filter(mbv_genotype_id != sample_cur)


head(final_df_match)
all(final_df_match_2$SRR_PT_ID == final_df_match_2$mbv_genotype_id)
final_df_match_2$match = if_else(final_df_match_2$sample_cur == final_df_match_2$mbv_genotype_id, "True", "False")
final_df_match_2$label = if_else(final_df_match_2$het_consistent_frac < 0.85 | final_df_match_2$hom_consistent_frac < 0.92, final_df_match_2$sample_cur, NA)






colnames(data) <- colnames(final_df_match_2)
final_df_match_2
final_df <- rbind(final_df_match_2, data)
MBV_one_sample <- ggplot(data = final_df, aes(x = het_consistent_frac, hom_consistent_frac, colour = match)) + geom_point(size = 2) + theme_bw() +
  scale_colour_manual(values = c("steelblue", "darkred")) + labs(x = "Concordance at HET", y = "Concordance at HOM", colour = "DNA/RNA\nmatch") +
  ylim(0,1) + xlim(0,1) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 










# It appears that patient 3 has amplification bias - need to investigate this
# Specifcially SRR12609781 from timepoint 1
library("GGally")


# read in filtered genes from DESEQ2
data_genes <- read.table(file = "/home/workspace/jogrady/heQTL/results/DESEQ2/DESEQ2_results_LFC_0.2.txt", sep = "\t")


# Read in all expression data
counts0 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T0.txt", sep = "\t", header = T, row.names = 1))
counts1 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T1.txt", sep = "\t", header = T, row.names = 1))
counts2 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T2.txt", sep = "\t", header = T, row.names = 1))
counts3 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T3.txt", sep = "\t", header = T, row.names = 1))
counts4 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T4.txt", sep = "\t", header = T, row.names = 1))


TPM_bulk = function(counts,lengths) {
  # counts = count matrix of bulk data
  # lenghts = vector of gene lenghts with names corresponding to gene names 
  
  # Intersect for common - cannot TPM normalise genes with no length information
  A = intersect(rownames(counts),names(lengths))
  counts = counts[A,]
  lengths = lengths[A]
  print(all(rownames(counts) == names(lengths)))
  # Perform TPM normalisation
  rate = counts/ lengths
  apply(rate,2,function(x) 1e6*x/sum(x))
  # return
  return(rate)
}



# Necessary files
# Need to convert gene ids to gene symbols
# Get in gene symbol information
symbols <- fread("/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf") %>% filter(V3 == "gene") # select(V9)
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)
symbols$length <- abs(symbols$V4 - symbols$V5)
symbols <- symbols %>% dplyr::select(gene_id, gene_name, length)

genes = unique(data_genes$gene_id)


symbols = symbols %>% filter(gene_id %in% genes)
symbols
# Set up lengths data frame
lengths = symbols$length
names(lengths) = symbols$gene_id


# filter timepoint 1 for expressed genes
intersect(rownames(counts1), names(lengths))

counts1_tpm <- TPM_bulk(counts1, lengths)
dim(counts1_tpm)

data_genes$gene_id <- gsub("\\..*", "", data_genes$gene_id)  


library(GGally)



correlation_plot_res <- cor.mtest(log(counts1_tpm +1), method = "spearman", conf.level = 0.95, exact = FALSE)

correlation_plot_res
correlation_plot <- cor(log(counts1_tpm + 1))
correlation_plot
library(corrplot)
library(grid)
library(gridGraphics)
library(patchwork)
plot_correlation <- wrap_elements(~corrplot(correlation_plot, p.mat = correlation_plot_res$p, method = 'square', order = "original", type = 'upper', insig='blank',
         addCoef.col ='black', number.cex = 0.5, tl.col = 'black', diag=FALSE, is.corr = FALSE, col.lim = c(min(correlation_plot), max(correlation_plot)), col = COL2("RdBu")))


coldata_1 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T1/T1_metadata.txt", sep = "\t", row.names = 1, header = T))         
sample <- rep(colnames(counts1), 1)
colnames(counts1)
Condition = rep("Timepoint_1", 48)
pheno <- cbind(coldata_1,sample, Condition)

dds <- DESeqDataSetFromMatrix(countData = counts1, colData =  pheno, design = ~ sample)
dds_vst <- vst(dds)

pcaData <- plotPCA(dds_vst, intgroup=c("Condition"), returnData=TRUE, ntop = 1500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$col = ifelse(pcaData$name == "P_3", "darkred", "black")
pcaData
pca <- ggplot(pcaData, aes(PC1, PC2, label = name)) +
  geom_point(size=3, colour = pcaData$col) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(nudge_y = 1,label = pcaData$name, col = "black", max.overlaps = 30) + theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 90, hjust = 0.5),
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 11),
        axis.text.y = element_text(colour = "black", size = 11)) 

cor_plot$arg$ty
MBV
pca
library(cowplot)
plot_correlation
grid.echo()
P1 <- grid.grab()


package.version("ggplot2")
corrplot(correlation_plot, p.mat = correlation_plot_res$p, method = 'square', order = "original", type = 'upper', insig='blank',
                              addCoef.col ='black', number.cex = 0.5, tl.col = 'black', diag=FALSE, is.corr = FALSE, col.lim = c(min(correlation_plot), max(correlation_plot)), col = COL2("RdBu"))



library("patchwork")
par(mfrow=c(2,2))

dev.off()
library("ggcorrplot")
plot = ggcorrplot(correlation_plot,
           hc.order = FALSE,
           type = "lower",
           outline.color = "white", lab = FALSE, ggtheme = ggplot2::theme_bw(), lab_size = 2.2) + ggplot2::scale_fill_gradient2(
             low = "lightblue",
             high = "darkred",
             mid = "steelblue",
             midpoint = 0.92,
             limit = c(0.84, 1),
             space = "Lab",
             name = "Correlation") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6.5),
                                           axis.text.y = element_text(size = 6.5))
plot

middle_row <- plot_grid(plot, pca, labels = c('B', 'C'), label_size = 12)

plot_grid(MBV, middle_row, labels = c('A', ''), label_size = 12, ncol = 1)
ggsave("/home/workspace/jogrady/heQTL/results/eQTL/MBV/Initial_MBV_correlation_PCA_analysis.pdf", width = 15, height = 15)



final_df_match <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match <- rbind(final_df_match, data)
}


final_df_match_2 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_2) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_1) {
  if (sample == "SRR12609781") {
    
      sample_cur = as.character(sample)
      data <- fread(paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/WASP/dup_filtered/sorted/",sample_cur,".mbv_output.txt"))
      best_match <- mbvFindBestMatch(data)
      data <- cbind(sample_cur, best_match)
      final_df_match_2 <- rbind(final_df_match_2, data)
  } else {
    sample_cur = as.character(sample)
    data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
    best_match <- mbvFindBestMatch(data)
    data <- cbind(sample_cur, best_match)
    final_df_match_2 <- rbind(final_df_match_2, data)
  }
}

final_df_match_3 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_3) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_2) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_3 <- rbind(final_df_match_3, data)
}

final_df_match_4 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_4) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_3) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_4 <- rbind(final_df_match_4, data)
}

final_df_match_5 <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match_5) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")
for (sample in string_vector_4) {
  sample_cur = as.character(sample)
  data <- fread(paste0("/home/workspace/jogrady/heQTL/results/eQTL/sample_check/",sample_cur,".mbv_output.txt"))
  best_match <- mbvFindBestMatch(data)
  data <- cbind(sample_cur, best_match)
  final_df_match_5 <- rbind(final_df_match_5, data)
}




final_df_match$SRR_PT_ID <- ids_timepoint_0
final_df_match_2$SRR_PT_ID <- ids_timepoint_1
final_df_match_3$SRR_PT_ID <- ids_timepoint_2
final_df_match_4$SRR_PT_ID <- ids_timepoint_3
final_df_match_5$SRR_PT_ID <- ids_timepoint_4

all(final_df_match$mbv_genotype_id == final_df_match$SRR_PT_ID)
all(final_df_match_2$mbv_genotype_id == final_df_match_2$SRR_PT_ID)
all(final_df_match_3$mbv_genotype_id == final_df_match_3$SRR_PT_ID)
all(final_df_match_4$mbv_genotype_id == final_df_match_4$SRR_PT_ID)
all(final_df_match_5$mbv_genotype_id == final_df_match_5$SRR_PT_ID)
head(final_df_match)


data_final <- rbind(final_df_match,
                    final_df_match_2,
                    final_df_match_3,
                    final_df_match_4,
                    final_df_match_5)
MBV
data_final <- data_final %>% select(1,2,3,4,8)
data_final <- data_final %>% pivot_longer(, cols = -c("sample_cur", "SRR_PT_ID", "mbv_genotype_id"), names_to ="Metric", values_to = "Proportion")
data_final$Metric = factor(data_final$Metric, levels = c("hom_consistent_frac", "het_consistent_frac"), labels = c("Consistent Homozygous", "Consistent Heterozygous"))
data_final$Label <- if_else(data_final$sample_cur == "SRR12609781" & data_final$Metric == "Consistent Heterozygous", data_final$sample_cur, NA)
MBV_updated <- ggplot(data = data_final, aes(x = SRR_PT_ID, y = Proportion, colour = Metric, group = interaction(SRR_PT_ID, Metric))) + 
  geom_boxplot(alpha = 0.1, width = 0.75, outlier.color = NA) + 
  geom_point(position = position_dodge(width=0.75)) + 
  geom_text_repel(nudge_x = 0, nudge_y = 0.07,label = data_final$Label, col = "black") +
  ylim(0.4,1) +
  theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 90, hjust = 0.5),
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 11),
        axis.text.y = element_text(colour = "black", size = 11)) 

bottom_row <- plot_grid(MBV_updated, labels = c('D'), label_size = 12)

plot_grid(MBV, middle_row, bottom_row, labels = c('A', '', ''), label_size = 12, ncol = 1)
ggsave("/home/workspace/jogrady/heQTL/results/eQTL/MBV/FINAL_MBV_correlation_PCA_analysis.pdf", width = 15, height = 15)
