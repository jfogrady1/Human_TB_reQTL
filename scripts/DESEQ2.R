# Script for time series DE analysis
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(sva)
library(bladderbatch)
library(pamr)
library(limma)
library(edgeR)
library(splines)
library(data.table)
library(ggrepel)
library(DEGreport)
library("BiocParallel")
library(RColorBrewer)
library("UpSetR")



plotPCA.mystyle <- function(object, intgroup="condition", ntop=500, returnData=FALSE, pcs = c(1,2)) {
  stopifnot(length(pcs) == 2)    ### added this to check number of PCs ####
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
 
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  ########## Here we just use the pcs object passed by the end user ####
  d <- data.frame(PC1=pca$x[,pcs[1]], PC2=pca$x[,pcs[2]], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
}
# Read in the data for each timepoint

# Read in all of the files
counts0 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T0.txt", sep = "\t", header = T, row.names = 1))
counts1 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T1.txt", sep = "\t", header = T, row.names = 1))
counts2 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T2.txt", sep = "\t", header = T, row.names = 1))
counts3 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T3.txt", sep = "\t", header = T, row.names = 1))
counts4 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_T4.txt", sep = "\t", header = T, row.names = 1))

# Read in covariate data
coldata_0 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T0/T0_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_1 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T1/T1_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_2 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T2/T2_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_3 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T3/T3_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_4 <- as.matrix(read.table("/home/workspace/jogrady/heQTL/data/T4/T4_metadata.txt", sep = "\t", row.names = 1, header = T))
coldata_3 <- coldata_3[,-3]
coldata_4 <- coldata_4[,-c(2,4,6)]

# make column for merged dataset
Time0 <- "T0"
Time1 <- "T1"
Time2 <- "T2"
Time3 <- "T3"
Time4 <- "T4"


coldata_0 <- cbind(coldata_0,Time0)
coldata_1 <- cbind(coldata_1,Time1)
coldata_2 <- cbind(coldata_2,Time2)
coldata_3 <- cbind(coldata_3,Time3)
coldata_4 <- cbind(coldata_4,Time4)


rownames(coldata_1) <- paste0(rownames(coldata_1),"_T1")
rownames(coldata_2) <- paste0(rownames(coldata_2),"_T2")
rownames(coldata_3) <- paste0(rownames(coldata_3),"_T3")
rownames(coldata_4) <- paste0(rownames(coldata_4),"_T4")
head(coldata_1)
colnames(counts1) <- paste0(colnames(counts1),"_T1")
colnames(counts2) <- paste0(colnames(counts2),"_T2")
colnames(counts3) <- paste0(colnames(counts3),"_T3")
colnames(counts4) <- paste0(colnames(counts4),"_T4")


# Merge for expression data
edata <- cbind(counts0, counts1, counts2, counts3, counts4)
pheno <- rbind(coldata_0, coldata_1, coldata_2, coldata_3, coldata_4)
sample <- rep(colnames(counts0), 5)
pheno <- cbind(pheno, sample)
head(pheno)
tail(pheno)

pheno <- as.data.frame(pheno)
head(pheno)
pheno$sample <- factor(pheno$sample, levels = rownames(coldata_0))
colnames(pheno)[9] <- "Time"
pheno$Time <- factor(pheno$Time, levels = c("T0", "T1", "T2", "T3", "T4"), labels = c("0","1","2","3","4"))

#-------------------------------------------------------------------------------------------------------------------------------



# Pairwise comparisons

#-------------------------------------------------------------------------------------------------------------------------------

# Set up the dds object
dds <- DESeqDataSetFromMatrix(countData = edata, colData = pheno, design = ~ sample + Time)



# DESEQ analysis for 4 pairewise comparisons all to the reference set before treatment
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 48
dds <- dds[filter,]


dds_pair <- DESeq(dds)

dds_pair
resultsNames(dds_pair)


dds_1_V_0 <- lfcShrink(dds_pair, coef="Time_1_vs_0", type="apeglm")
dds_2_V_0 <- lfcShrink(dds_pair, coef="Time_2_vs_0", type="apeglm")
dds_3_V_0 <- lfcShrink(dds_pair, coef="Time_3_vs_0", type="apeglm")
dds_4_V_0 <- lfcShrink(dds_pair, coef="Time_4_vs_0", type="apeglm")


# When using an alternative hypothesis, we do not shrink the LFCs
dds_1_V_0_lfc = results(dds_pair, lfcThreshold = 0.2, altHypothesis = "greaterAbs", name = "Time_1_vs_0")
dds_2_V_0_lfc = results(dds_pair, lfcThreshold = 0.2, altHypothesis = "greaterAbs", name = "Time_2_vs_0")
dds_3_V_0_lfc = results(dds_pair, lfcThreshold = 0.2, altHypothesis = "greaterAbs", name = "Time_3_vs_0")
dds_4_V_0_lfc = results(dds_pair, lfcThreshold = 0.2, altHypothesis = "greaterAbs", name = "Time_4_vs_0")


# Get in gene symbol information
symbols <- fread("/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf") %>% filter(V3 == "gene") %>% select(V9)
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols <- symbols %>% select(1,3)
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)


# Get geneID into a column in each of the results matrices
dds_1_V_0_lfc <- dds_1_V_0_lfc %>% as.data.frame %>% mutate(gene_id = rownames(dds_1_V_0_lfc))
dds_2_V_0_lfc <- dds_2_V_0_lfc %>% as.data.frame %>% mutate(gene_id = rownames(dds_2_V_0_lfc))
dds_3_V_0_lfc <- dds_3_V_0_lfc %>% as.data.frame %>% mutate(gene_id = rownames(dds_3_V_0_lfc))
dds_4_V_0_lfc <- dds_4_V_0_lfc %>% as.data.frame %>% mutate(gene_id = rownames(dds_4_V_0_lfc))

# Now join in the symbols
dds_1_V_0_lfc <- left_join(dds_1_V_0_lfc, symbols)
dds_2_V_0_lfc <- left_join(dds_2_V_0_lfc, symbols)
dds_3_V_0_lfc <- left_join(dds_3_V_0_lfc, symbols)
dds_4_V_0_lfc <- left_join(dds_4_V_0_lfc, symbols)


# Get deleniating marks for facet wrap
dds_1_V_0_lfc$signif <- if_else(dds_1_V_0_lfc$padj < 0.05, TRUE, FALSE)
dds_2_V_0_lfc$signif <- if_else(dds_2_V_0_lfc$padj < 0.05, TRUE, FALSE)
dds_3_V_0_lfc$signif <- if_else(dds_3_V_0_lfc$padj < 0.05, TRUE, FALSE)
dds_4_V_0_lfc$signif <- if_else(dds_4_V_0_lfc$padj < 0.05, TRUE, FALSE)


# Get colour 
dds_1_V_0_lfc$diffexpressed <- "Not DE"
# if log2Foldchange > 0.6 and pvalue < 0.01, set as "UP" 
dds_1_V_0_lfc$diffexpressed[dds_1_V_0_lfc$log2FoldChange > 0 & dds_1_V_0_lfc$padj < 0.05] <- "DE Up"
# if log2Foldchange < -0.6 and pvalue < 0.01, set as "DOWN"
dds_1_V_0_lfc$diffexpressed[dds_1_V_0_lfc$log2FoldChange < 0 & dds_1_V_0_lfc$padj < 0.05] <- "DE Down"


dds_2_V_0_lfc$diffexpressed <- "Not DE"
# if log2Foldchange > 0.6 and pvalue < 0.01, set as "UP" 
dds_2_V_0_lfc$diffexpressed[dds_2_V_0_lfc$log2FoldChange > 0 & dds_2_V_0_lfc$padj < 0.05] <- "DE Up"
# if log2Foldchange < -0.6 and pvalue < 0.01, set as "DOWN"
dds_2_V_0_lfc$diffexpressed[dds_2_V_0_lfc$log2FoldChange < 0 & dds_2_V_0_lfc$padj < 0.05] <- "DE Down"

dds_3_V_0_lfc$diffexpressed <- "Not DE"
# if log2Foldchange > 0.6 and pvalue < 0.01, set as "UP" 
dds_3_V_0_lfc$diffexpressed[dds_3_V_0_lfc$log2FoldChange > 0 & dds_3_V_0_lfc$padj < 0.05] <- "DE Up"
# if log2Foldchange < -0.6 and pvalue < 0.01, set as "DOWN"
dds_3_V_0_lfc$diffexpressed[dds_3_V_0_lfc$log2FoldChange < 0 & dds_3_V_0_lfc$padj < 0.05] <- "DE Down"


dds_4_V_0_lfc$diffexpressed <- "Not DE"
# if log2Foldchange > 0.6 and pvalue < 0.01, set as "UP" 
dds_4_V_0_lfc$diffexpressed[dds_4_V_0_lfc$log2FoldChange > 0 & dds_4_V_0_lfc$padj < 0.05] <- "DE Up"
# if log2Foldchange < -0.6 and pvalue < 0.01, set as "DOWN"
dds_4_V_0_lfc$diffexpressed[dds_4_V_0_lfc$log2FoldChange < 0 & dds_4_V_0_lfc$padj < 0.05] <- "DE Down"


# Get timepoint for facet
dds_1_V_0_lfc$Contrast <- "T1_V_T0"
dds_2_V_0_lfc$Contrast <- "T2_V_T0"
dds_3_V_0_lfc$Contrast <- "T3_V_T0"
dds_4_V_0_lfc$Contrast <- "T4_V_T0"

dds_df <- rbind(dds_1_V_0_lfc,dds_2_V_0_lfc, dds_3_V_0_lfc, dds_4_V_0_lfc)

dds_df_totals = dds_df %>% group_by(Contrast) %>% filter(signif == T) %>% summarise(total_down = sum(diffexpressed == "DE Down"),
                                                                                    total_up = sum(diffexpressed == "DE Up"))


dds_df_top_10 = dds_df %>% group_by(Contrast) %>% filter(signif == T) %>% group_by(Contrast, diffexpressed) %>% arrange(padj) %>% slice_min(order_by = padj,n = 10)



dds_df_top_10 <- dds_df_top_10 %>% select(gene_id, gene_name, Contrast)
dds_df_top_10 <- dds_df_top_10[,2:4]

colnames(dds_df_top_10)[2] <- "Label"

dds_df <- left_join(dds_df, dds_df_top_10)


write.table(dds_df, file = "/home/workspace/jogrady/heQTL/results/DESEQ2/DESEQ2_results_LFC_0.2.txt", sep = "\t", quote = F)
colnames(dds_df)
# Plotting
ggplot(data=dds_df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = Label)) +
    geom_point(size = 1, alpha = 0.5, shape = 16) + 
    scale_color_manual("Comparison", values=c("#2166ac", "#b2182b", "grey")) +
    labs(x=expression(log[2]("fold change")),
       y=expression(-log[10](italic(P)[adj]))) +
       xlim(c(-5.0,5)) +
    geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
    theme_bw() + # changed to a classic theme for a clean look
    theme(axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black"),
        axis.title.x = element_text(size = 21, color = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        panel.grid.minor = element_blank()) +
        facet_wrap(~Contrast, nrow = 2, scales = "free_y") +
        geom_text_repel(aes(label = Label), max.overlaps = 25) +
         geom_text(data = dds_df_totals, aes(x = 3.5, y = -log10(0.00005), label = paste0("Up = ",total_up), color =NULL,group= NULL)) +
         geom_text(data = dds_df_totals, aes(x = -3.5, y = -log10(0.00005), label = paste0(" Down = ",total_down), color = NULL,group= NULL)) +
    theme(
      strip.text.x = element_text(
        size = 15, color = "black", face = "bold"),
        strip.background = element_rect(
     color="black", size=1.5, linetype="blank"))



ggsave("/home/workspace/jogrady/heQTL/results/DESEQ2/Volcano_plot_top_10_lfc_0.2.pdf", width = 10, height = 10)
min(dds_df$log2FoldChange)

library(tidyverse)

T1 <- dds_1_V_0_lfc %>% filter(signif == TRUE)
T2 <- dds_2_V_0_lfc %>% filter(signif == TRUE)
T3 <- dds_3_V_0_lfc %>% filter(signif == TRUE)
T4 <- dds_4_V_0_lfc %>% filter(signif == TRUE)


library(UpSetR)



ALL_signif <- rbind(T1, T2, T3, T4)



T0_V_T1_up = T1 %>% filter(diffexpressed == "DE Up") %>% select(gene_name) 
T0_V_T2_up = T2 %>% filter(diffexpressed == "DE Up") %>% select(gene_name)
T0_V_T3_up = T3 %>% filter(diffexpressed == "DE Up") %>% select(gene_name) 
T0_V_T4_up = T4 %>% filter(diffexpressed == "DE Up") %>% select(gene_name)
T0_V_T1_down = T1 %>% filter(diffexpressed == "DE Down") %>% select(gene_name)
T0_V_T2_down = T2 %>% filter(diffexpressed == "DE Down") %>% select(gene_name)
T0_V_T3_down = T3 %>% filter(diffexpressed == "DE Down") %>% select(gene_name)
T0_V_T4_down = T4 %>% filter(diffexpressed == "DE Down") %>% select(gene_name)

#save.image(file = "DESEQ2.RData")
#load("DESEQ2.RData")
listInput <- list(T0_V_T1_up = T0_V_T1_up$gene_name,
                  T0_V_T2_up =  T0_V_T2_up$gene_name,
                  T0_V_T3_up = T0_V_T3_up$gene_name,
                  T0_V_T4_up =  T0_V_T4_up$gene_name,
                  T0_V_T1_down =  T0_V_T1_down$gene_name,
                  T0_V_T2_down =  T0_V_T2_down$gene_name,
                  T0_V_T3_down =  T0_V_T3_down$gene_name,
                  T0_V_T4_down = T0_V_T4_down$gene_name)
names(listInput)
library("pals")
library(RColorBrewer)
pdf("/home/workspace/jogrady/heQTL/results/DESEQ2/Upset_DE_Genes_FDR0.05_lfc_0.2.pdf", width = 15, height = 12)
upset(fromList(listInput), order.by = "freq", nintersects = 25, nsets = 10, sets = names(listInput), query.legend = "top",
     point.size = 4, line.size = 2,  text.scale = c(4, 2.5, 1, 1, 2, 2.5), keep.order = TRUE,sets.x.label = "DE genes", mainbar.y.label = "DE gene intersections", sets.bar.color = brewer.rdylbu(n=8),
     queries = list(list(query = intersects, 
     params = list("T0_V_T1_up", "T0_V_T2_up", "T0_V_T3_up", "T0_V_T4_up"), color = "darkred", active = T, 
     query.name = "Consistently up"),
     list(query = intersects, 
     params = list("T0_V_T1_down", "T0_V_T2_down", "T0_V_T3_down", "T0_V_T4_down"), color = "darkblue", active = T, 
     query.name = "Consistently down")))
dev.off()



T0_V_T2_up
T0_V_T2_up_gprofiler <- T2 %>% filter(signif == TRUE) %>% 
  filter(diffexpressed == "DE Up") %>% 
  arrange(desc(padj)) %>% slice_min(., padj, n = 300)

T0_V_T2_down_gprofiler <- T2 %>% filter(signif == TRUE) %>% 
  filter(diffexpressed == "DE Down") %>% 
  arrange(desc(padj)) %>% slice_min(., padj, n = 300)

T0_V_T3_up_gprofiler <- T3 %>% filter(signif == TRUE) %>% 
  filter(diffexpressed == "DE Up") %>% 
  arrange(desc(padj)) %>% slice_min(., padj, n = 300)

T0_V_T3_down_gprofiler <- T3 %>% filter(signif == TRUE) %>% 
  filter(diffexpressed == "DE Down") %>% 
  arrange(desc(padj)) %>% slice_min(., padj, n = 300)

T0_V_T4_up_gprofiler <- T4 %>% filter(signif == TRUE) %>% 
  filter(diffexpressed == "DE Up") %>% 
  arrange(desc(padj)) %>% slice_min(., padj, n = 300)

T0_V_T4_down_gprofiler <- T4 %>% filter(signif == TRUE) %>% 
  filter(diffexpressed == "DE Down") %>% 
  arrange(desc(padj)) %>% slice_min(., padj, n = 300)

dds_1_V_0_lfc$gene_id
dds_2_V_0_lfc$gene_id

# Now onto Gprofiler
library(gprofiler2)
# order the query based on adjusted p value, to do this we will select the most significant association for common genes

background <- dds_1_V_0_lfc$gene_name

gostres_T2_up <- gost(T0_V_T2_up_gprofiler$gene_name,
ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                custom_bg = background,
                sources = c("GO", "REAC", "KEGG", "WP"),
                organism = "hsapiens")

gostres_T2_down <- gost(T0_V_T2_down_gprofiler$gene_name,
                      ordered_query = TRUE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                      user_threshold = 0.05, correction_method = "fdr", 
                      custom_bg = background,
                      sources = c("GO", "REAC", "KEGG", "WP"),
                      organism = "hsapiens")


gostres_T3_up <- gost(T0_V_T3_up_gprofiler$gene_name,
                        ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        custom_bg = background,
                        sources = c("GO", "REAC", "KEGG", "WP"),
                        organism = "hsapiens")


gostres_T3_down <- gost(T0_V_T3_down_gprofiler$gene_name,
                        ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        custom_bg = background,
                        sources = c("GO", "REAC", "KEGG", "WP"),
                        organism = "hsapiens")


gostres_T4_up <- gost(T0_V_T4_up_gprofiler$gene_name,
                        ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        custom_bg = background,
                        sources = c("GO", "REAC", "KEGG", "WP"),
                        organism = "hsapiens")


gostres_T4_down <- gost(T0_V_T4_down_gprofiler$gene_name,
                        ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = TRUE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        custom_bg = background,
                        sources = c("GO", "REAC", "KEGG", "WP"),
                        organism = "hsapiens")



# Real figure - T3
gostres_T3_up_res <- gostres_T3_up$result
gostres_T3_up_res$query <- "T3_V_T0_up"

gostres_T3_down_res <- gostres_T3_down$result
gostres_T3_down_res$query <- "T3_V_T0_down"



terms <- c("defense response",
           "innate immune response",
           "B cell differentiation",
           "B cell proliferation",
           "B cell activation",
           "Extrafollicular and follicular B cell activation by SARS CoV 2",
           "Hematopoietic cell lineage",
           "response to type II interferon",
           "pyroptosis",
           "interleukin-1 production",
           "Immune response to tuberculosis",
           "NOD-like receptor signaling pathway",
           "canonical inflammasome complex",
           "Toll-like receptor 4 binding",
           "high-affinity IgG receptor activity",
           "Th17 cell differentiation",
           "regulation of transformation of host cell by virus",
           "African trypanosomiasis",
           "The AIM2 inflammasome",
           "Type II interferon signaling",
           "IPAF inflammasome complex"
)





figure_2_res <- rbind(gostres_T3_up_res, gostres_T3_down_res)
figure_2_res$query <- factor(figure_2_res$query, levels = c("T3_V_T0_up", "T3_V_T0_down"))
figure_2_res  <- figure_2_res  %>%
  mutate(label = ifelse(term_name %in% terms, term_name, NA))
figure_2_res <- figure_2_res %>% mutate(alpha_value = ifelse(term_name %in% terms, 1, 0.5))
color_palette <- viridis(5)
pos <- position_jitter(width = 0.3, seed = 3)
ggplot(figure_2_res, aes(x = source, y = -log10(p_value))) +
  geom_jitter(aes(color = source), alpha = figure_2_res$alpha_value, position = pos, size = 3) +  # Use shape instead of color
  scale_color_brewer(palette = "Dark2" , name = "Source") +  # Apply the color palette
  labs(y = expression(-log[10](italic(P)[adj])), x = "") +
  theme_bw() + facet_wrap(~ query, scales = "free_y") +
  labs(size = "Intersection\nSize") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", face = "bold"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15)) +
  geom_text_repel(aes(x=source, y=-log10(p_value)), label = figure_2_res$label, max.overlaps = 20,
                   size=3.0, color='black', fontface = "bold",
                   fill='#FFFFFF33',
                   position = pos,
                  segment.linewidth = 0.2,
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 4
  ) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("/home/workspace/jogrady/heQTL/results/DESEQ2/Gprofiler_results_T3_up_and_down.pdf", width = 15, height = 10)


# Focus on genes which are up - may be markers of tuberculosis treatment (successful if related back to a phenotype)

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



head(symbols)
symbols$length
symbols$gene_id <- gsub("\\..*", "", symbols$gene_id)
symbols_length <- as.numeric(symbols$length)
symbols_length
names(symbols_length) <- symbols$gene_id
names(symbols_length)




B = intersect(T1$gene_name, intersect(T2$gene_name, intersect(T3$gene_name, T4$gene_name)))
B
T1_intersect = T1 %>% filter(gene_name %in% B) %>% mutate(T0 = 0) %>% select(gene_name, T0, log2FoldChange)
T2_intersect = T2 %>% filter(gene_name %in% B) %>% select(gene_name, log2FoldChange)
T3_intersect = T3 %>% filter(gene_name %in% B) %>% select(gene_name, log2FoldChange)
T4_intersect = T4 %>% filter(gene_name %in% B) %>% select(gene_name, log2FoldChange)

T1_intersect
ALL_LFC <- cbind(T1_intersect, T2_intersect$log2FoldChange, T3_intersect$log2FoldChange, T4_intersect$log2FoldChange)
colnames(ALL_LFC) <- c("gene_name", "T0", "T1VT0", "T2VT0", "T3VT0", "T4VT0")

ALL_LFC <- ALL_LFC %>% pivot_longer(cols = colnames(ALL_LFC)[2:6], names_to = "Time", values_to = "LFC")





library(dplyr)
library(lubridate)
library(ggplot2)
library("ggpp")
library(grid)
my_palette = c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026")
ALL_LFC$Label = if_else(ALL_LFC$Time == "T4VT0", ALL_LFC$gene_name, NA)

ALL_LFC <- ALL_LFC %>% arrange(c(Time)) %>% arrange(LFC)
head(ALL_LFC)

y_points = seq(-3, 2, by = 0.04)
y_points
x_points = rep(c(5.1, 5.4), 52)
x_points <- x_points[1:103]
y_points <- sort(y_points, decreasing = TRUE)
y_points <- y_points[1:103]
x_points

ALL_LFC <- ALL_LFC %>%
  arrange(desc(Time), desc(LFC))

ALL_LFC_points = ALL_LFC %>% filter(Time == "T4VT0") %>% arrange(desc(LFC)) %>% mutate(x_coord = x_points, 
                                                                      y_coord = y_points) %>% select(Time, Label, x_coord, y_coord)
y_points

ggplot(ALL_LFC, aes(x = Time, y = LFC, fill = Time, group = gene_name, label = Label)) + geom_line(linetype = "twodash", linewidth = 0.25) + geom_point(shape = 21)+ theme_bw() + # changed to a classic theme for a clean look
  theme(axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black"),
        axis.title.x = element_text(size = 21, color = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        panel.grid.minor = element_blank()) +
  geom_text_s(aes(label = Label),
              size = 2.5,
              position = position_nudge_to(x = c(5.3, 5.1)),
              box.padding = 0.5, segment.linewidth = 0.1) +
  scale_fill_manual(values = my_palette)

ggsave("/home/workspace/jogrady/heQTL/results/DESEQ2/LFC_tracking_contrasts.pdf", width = 12, height = 10)


genes_of_interest = c("BAK1",
                      "SCARF1",
                      "ICAM1",
                      "LIMK1",
                      "SEPT4",
                      "FAM20A",
                      "CD274",
                      "NRN1",
                      "IGF2BP3",
                      "CASP5",
                      "TIFA",
                      "FCGR1A",
                      "GBP5",
                      "C1QC",
                      "SDC3",
                      "CDCP1",
                      "C1QB",
                      "P2RY14",
                      "GRAMD1C",
                      "GBP6",
                      "SOCS3",
                      "PRTN3",
                      "PDCD1LG2",
                      "FCGR1B",
                      "GRIN3A",
                      "GBP1P1",
                      "FCGR1CP")



B # my list
genes_of_interest  # genes_of interest
background # background




A_fisher <- sum(B %in% genes_of_interest)                                # Genes in both identified and reference
B_fisher <- sum((background %in% genes_of_interest) & !(background %in% B))  # In reference, not in identified
C_fisher <- sum((background %in% B) & !(background %in% genes_of_interest))  # In identified, not in reference
D_fisher <- length(background) - A_fisher - B_fisher - C_fisher  


fisher_result <- fisher.test(matrix(c(A_fisher, B_fisher, C_fisher, D_fisher), nrow = 2), alternative = "greater")
fisher_result







# PCA now
vsd <- vst(dds_pair, blind=FALSE)
pcaData <- plotPCA.mystyle(vsd, intgroup=c("Time", "gender", "subgroup", "smear_results"), returnData=TRUE, pcs = c(1,2))
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Time, shape=gender)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[2],"% variance")) + 
  coord_fixed()




common <- intersect(T4$gene_name, intersect(T3$gene_name, intersect(T1$gene_name, T2$gene_name)))
View(dds_df)
head(dds_df)
intersection <- dds_df %>% group_by(gene_name) %>% filter(gene_name %in% common)
intersection <- intersection %>% group_by(gene_name) %>% reframe(lFC = log2FoldChange)
head(intersection)

intersection$Time <- rep(c("T1", "T2", "T3", "T4"), 961)

intersection$Time <- factor(intersection$Time)
ggplot(data = intersection, aes(x = Time, y = lFC, group = gene_name)) + geom_line()
duplicated(intersection$gene_name)
table(intersection$gene_name) == 8
View(intersection)
length(unique(dds_df$gene_name))
dim(intersection)
16013 * 4
listInput <- list(T1 = T1$gene_name, T2 = T2$gene_name, T3 = T3$gene_name, T4 = T4$gene_name)
upset(fromList(listInput), order.by = "freq", sets.bar.color = c("#542788","#2166ac", "#b2182b", "yellow"),
sets.x.label = "cis-eGenes", point.size = 4, line.size = 2,
mainbar.y.label = "cis-eGene intersections",
text.scale = 2.5, shade.alpha = 0.5)





ggsave("Results.pdf", width = 12, height = 12, dpi = 600)
