# MASHR for ieQTL results

# First get the top ieQTL for each Gene across timepoints

library(tidyverse)
library(data.table)
library(arrow)
library(ash)
library(ashr)
library(mashr)
library(gtable)
library(cowplot)
set.seed(34567)

args = commandArgs(trailingOnly=TRUE)
shift_legend <- function(p){

  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }

  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }

  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")

  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")

  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")

  return(gp)
}
set.seed(50)
MASH = function(celltype) {
  chromosomes = as.character(c(1:22))
  open_eGene = function(timepoints, chromosomes, cell_type = NULL){
    tensor_colnames_egene = c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_g","b_g","b_g_se",        
                              "pval_i","b_i","b_i_se","pval_gi","b_gi","b_gi_se")
    
    results_df_egene <- as.data.frame(matrix(ncol = 15, nrow = 0))
    
    colnames(results_df_egene) <- tensor_colnames_egene
    for (t in timepoints) {
      for (c in chromosomes) {
        data_temp = read_parquet(paste0("/home/workspace/jogrady/heQTL/work/ieQTL/", t, "_",cell_type,".cis_qtl_pairs.chr", c, ".parquet"))
        data_temp$timepoint <- t
        results_df_egene <- rbind(results_df_egene, data_temp)
      }
    }
    return(results_df_egene)
  }

  T0 = open_eGene(timepoints = "T0", chromosomes = c(1:22), cell_type = celltype)


  T0_nk <- open_eGene("T0", chromosomes, cell_type = celltype)
  T1_nk <- open_eGene("T1", chromosomes, cell_type = celltype)
  T2_nk <- open_eGene("T2", chromosomes, cell_type = celltype)
  T3_nk <- open_eGene("T3", chromosomes, cell_type = celltype)
  T4_nk <- open_eGene("T4", chromosomes, cell_type = celltype)

  problematic_variants <- c()
  T0_nk_variants <- T0_nk %>% filter(is.na(b_gi_se)) %>% select(variant_id) %>% as.vector()
  T1_nk_variants <- T1_nk %>% filter(is.na(b_gi_se)) %>% select(variant_id) %>% as.vector()
  T2_nk_variants <- T2_nk %>% filter(is.na(b_gi_se)) %>% select(variant_id) %>% as.vector()
  T3_nk_variants <- T3_nk %>% filter(is.na(b_gi_se)) %>% select(variant_id) %>% as.vector()
  T4_nk_variants <- T4_nk %>% filter(is.na(b_gi_se)) %>% select(variant_id) %>% as.vector()

  problematic_variants <- c(T0_nk_variants$variant_id, T1_nk_variants$variant_id, T2_nk_variants$variant_id, T3_nk_variants$variant_id, T4_nk_variants$variant_id)
  problematic_variants
  # Get common genes to all
  T0_nk$variant_gene <- paste0(T0_nk$variant_id, "-", T0_nk$phenotype_id)
  T1_nk$variant_gene <- paste0(T1_nk$variant_id, "-", T1_nk$phenotype_id)
  T2_nk$variant_gene <- paste0(T2_nk$variant_id, "-", T2_nk$phenotype_id)
  T3_nk$variant_gene <- paste0(T3_nk$variant_id, "-", T3_nk$phenotype_id)
  T4_nk$variant_gene <- paste0(T4_nk$variant_id, "-", T4_nk$phenotype_id)


  T0_nk <- T0_nk[!(T0_nk$variant_id %in% problematic_variants),]
  T1_nk <- T1_nk[!(T1_nk$variant_id %in% problematic_variants),]
  T2_nk <- T2_nk[!(T2_nk$variant_id %in% problematic_variants),]
  T3_nk <- T3_nk[!(T3_nk$variant_id %in% problematic_variants),]
  T4_nk <- T4_nk[!(T4_nk$variant_id %in% problematic_variants),]


  common_variant_gene = intersect(T0_nk$variant_gene, intersect(T1_nk$variant_gene, intersect(T2_nk$variant_gene, intersect(T3_nk$variant_gene, T4_nk$variant_gene))))
  length(common_variant_gene)



  # 15109 common genes
  T0_nk <- T0_nk %>% filter(variant_gene %in% common_variant_gene)
  T1_nk <- T1_nk %>% filter(variant_gene %in% common_variant_gene)
  T2_nk <- T2_nk %>% filter(variant_gene %in% common_variant_gene)
  T3_nk <- T3_nk %>% filter(variant_gene %in% common_variant_gene)
  T4_nk <- T4_nk %>% filter(variant_gene %in% common_variant_gene)                         
  length(common_variant_gene)


  all(T0_nk$phenotype_id == T1_nk$phenotype_id)
  all(T1_nk$phenotype_id == T2_nk$phenotype_id)
  all(T2_nk$phenotype_id == T3_nk$phenotype_id)
  all(T3_nk$phenotype_id == T4_nk$phenotype_id)


  all(T0_nk$variant_id == T1_nk$variant_id)
  all(T1_nk$variant_id == T2_nk$variant_id)
  all(T2_nk$variant_id == T3_nk$variant_id)
  all(T3_nk$variant_id == T4_nk$variant_id)


  # need to drop NAs if the standard error is 0


  # get matrix of effect sizes
  bhat_df <- cbind(T0_nk$b_gi, T1_nk$b_gi, T2_nk$b_gi, T3_nk$b_gi, T4_nk$b_gi)
  bhat_df_matrix <- as.matrix(bhat_df)
  colnames(bhat_df_matrix) <- c("T0", "T1", "T2", "T3", "T4")
  rownames(bhat_df_matrix) <- paste0(T0_nk$variant_id,"-",T0_nk$phenotype_id)


  bhat_df_se <- cbind(T0_nk$b_gi_se, T1_nk$b_gi_se, T2_nk$b_gi_se, T3_nk$b_gi_se, T4_nk$b_gi_se)
  bhat_df_se_matrix <- as.matrix(bhat_df_se)
  colnames(bhat_df_se_matrix) <- c("T0", "T1", "T2", "T3", "T4")
  rownames(bhat_df_se_matrix) <- paste0(T0_nk$variant_id,"-",T0_nk$phenotype_id)

  head(bhat_df_matrix)

  # function to cycle through results and collate into single object
  open_eGene_signif = function(timepoints, cell_type = NULL){
    tensor_colnames_egene = c("phenotype_id",   "variant_id",     "start_distance", "af",             "ma_samples",     "ma_count",       
                              "pval_g",         "b_g",            "b_g_se",        "pval_i",         "b_i",            "b_i_se",         
                              "pval_gi",        "b_gi",           "b_gi_se",        "tests_emt",      "pval_emt",       "pval_adj_bh" )
    
    results_df_egene <- as.data.frame(matrix(ncol = 18, nrow = 0))
    
    colnames(results_df_egene) <- tensor_colnames_egene
    for (t in timepoints) {
      data_temp = read.table(paste0("/home/workspace/jogrady/heQTL/work/ieQTL/", t, "_", cell_type,".cis_qtl_top_assoc.txt.gz"), skip = 1)
      colnames(data_temp) <- tensor_colnames_egene
      data_temp$timepoint <- t
      results_df_egene <- rbind(results_df_egene, data_temp)
    }
    return(results_df_egene)
  }

  signif_eGenes = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = celltype)
  signif_eGenes <- signif_eGenes[!(signif_eGenes$variant_id %in% problematic_variants),]
  signif_eGenes$variant_gene <- paste0(signif_eGenes$variant_id, "-", signif_eGenes$phenotype_id)
  signif_eGenes <- signif_eGenes %>% filter(variant_gene %in% common_variant_gene)

  signif_eGenes #<- signif_eGenes %>% filter(is_eGene == TRUE)
  signif_eGenes$variant_gene <- paste0(signif_eGenes$variant_id,"-",signif_eGenes$phenotype_id)

  length(unique(signif_eGenes$phenotype_id))

  signif_variants = signif_eGenes %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.

  # Here are our lead unique SNps per gene accross all condiitons which are eQTLs
  strong.subset = as.vector(signif_variants$variant_gene)
  names(strong.subset) <- strong.subset

  length(strong.subset)
  T4_nk %>% filter(variant_id == "1:45684814:C:A")
  mashr_input <- list(Bhat = bhat_df_matrix,
                      Shat = bhat_df_se_matrix)
  random.subset = sample(1:nrow(mashr_input$Bhat), 200000)

  dim(bhat_df_matrix)
  head(bhat_df_matrix)
  dim(bhat_df_se_matrix)
  head(bhat_df_se_matrix)
  print(celltype)
  # learn the canonical covariance matrices
  data.temp = mash_set_data(mashr_input$Bhat[random.subset,],mashr_input$Shat[random.subset,])
  Vhat = estimate_null_correlation_simple(data.temp)

  rm(data.temp)

  data.random = mash_set_data(mashr_input$Bhat[random.subset,],mashr_input$Shat[random.subset,],V=Vhat)
  data.strong = mash_set_data(mashr_input$Bhat[strong.subset,],mashr_input$Shat[strong.subset,], V=Vhat)
  U.pca = cov_pca(data.strong,5)
  U.ed = cov_ed(data.strong, U.pca)

  print("HERE")

  U.c = cov_canonical(data.random)
  m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)


  print("COMPLETED")

  m


  m2_nk = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

  length(get_significant_results(m2_nk, thresh = 0.05, conditions = NULL, sig_fn = get_lfsr))

  thresh = 0.01
  lfsr.mash_nk <- m2_nk$result$lfsr
  sigmat <- (lfsr.mash_nk < thresh)
  nsig <- rowSums(sigmat)
  lfsr.mash.sig <- lfsr.mash_nk[nsig > 0,]
  dim(lfsr.mash.sig) 
  apply(lfsr.mash_nk<thresh, 2, sum)


  m2_posterior_nk = mash_compute_posterior_matrices(m,data.strong)

  updated_means <- m2_posterior_nk$PosteriorMean
  updated_sds <- m2_posterior_nk$PosteriorSD
  head(updated_means)
  updated_means %>% as.data.frame() %>% filter(T0 == min(T0))
  return(list(MASH = m2_nk, Posterior = m2_posterior_nk, Posterior_means = updated_means, Posterior_sds = updated_sds))
}


nk_cell_mash <- MASH(args[4])
Classical_Monocyte <- MASH(args[1])
Naïve_B_cell_mash <- MASH(celltype = args[3])
memory_CD8_T_cell_mash <- MASH(args[2])
Non_classical_Monocyte <- MASH(args[5])



length(get_significant_results(nk_cell_mash$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))
length(get_significant_results(Classical_Monocyte$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))
length(get_significant_results(Naïve_B_cell_mash$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))
length(get_significant_results(Non_classical_Monocyte$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))
length(get_significant_results(memory_CD8_T_cell_mash$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))






library(tidyverse)
library(viridis)

# NK cell
nk_signif <- get_significant_results(nk_cell_mash$MASH, thresh = 0.01, conditions = NULL)
nk_signif <- names(nk_signif)
nk_signif <- nk_cell_mash$MASH$result$lfsr[nk_signif,]
nk_signif_df <- as.data.frame(nk_signif)
nk_signif_df$Term <- rownames(nk_signif)
nk_signif_df <- pivot_longer(nk_signif_df, cols = colnames(nk_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
nk_signif_df$signif <- if_else(nk_signif_df$LFSR < 0.01, TRUE, FALSE)
nk_signif_df$cell <- "NK cell"



# Monocyte 
monocyte_signif <- get_significant_results(Classical_Monocyte$MASH, thresh = 0.01, conditions = NULL)
monocyte_signif <- names(monocyte_signif)
monocyte_signif <- Classical_Monocyte$MASH$result$lfsr[monocyte_signif,]
monocyte_signif_df <- as.data.frame(monocyte_signif)
monocyte_signif_df$Term <- rownames(monocyte_signif)
monocyte_signif_df <- pivot_longer(monocyte_signif_df, cols = colnames(monocyte_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
monocyte_signif_df$signif <- if_else(monocyte_signif_df$LFSR < 0.01, TRUE, FALSE)
monocyte_signif_df$cell <- "Classical Monocyte"


# B cell
b_signif <- get_significant_results(Naïve_B_cell_mash$MASH, thresh = 0.01, conditions = NULL)
b_signif <- names(b_signif)
b_signif <- Naïve_B_cell_mash$MASH$result$lfsr[b_signif,]
b_signif_df <- as.data.frame(b_signif)
b_signif_df$Term <- rownames(b_signif)
b_signif_df <- pivot_longer(b_signif_df, cols = colnames(b_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
b_signif_df$signif <- if_else(b_signif_df$LFSR < 0.01, TRUE, FALSE)
b_signif_df$cell <- "Naive B Cell"


#Non_classical_Monocyte
Non_classical_Monocytes <- get_significant_results(Non_classical_Monocyte$MASH, thresh = 0.01, conditions = NULL)
Non_classical_Monocytes <- names(Non_classical_Monocytes)
Non_classical_Monocytes <- Non_classical_Monocyte$MASH$result$lfsr[Non_classical_Monocytes,]

Non_classical_Monocytes_df <- as.data.frame(Non_classical_Monocytes)
Non_classical_Monocytes_df$Term <- rownames(Non_classical_Monocytes)



Non_classical_Monocytes_df <- pivot_longer(Non_classical_Monocytes_df, cols = colnames(Non_classical_Monocytes)[1:5], names_to = "Timepoint", values_to = "LFSR")
Non_classical_Monocytes_df$signif <- if_else(Non_classical_Monocytes_df$LFSR < 0.01, TRUE, FALSE)
Non_classical_Monocytes_df$cell <- "Non_classical_Monocyte"

#CD8+ tcells
cd8_signif <- get_significant_results(memory_CD8_T_cell_mash$MASH, thresh = 0.01, conditions = NULL)
cd8_signif <- names(cd8_signif)
cd8_signif <- memory_CD8_T_cell_mash$MASH$result$lfsr[cd8_signif,]
cd8_signif_df <- as.data.frame(cd8_signif)
cd8_signif_df$Term <- rownames(cd8_signif)
cd8_signif_df <- pivot_longer(cd8_signif_df, cols = colnames(cd8_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
cd8_signif_df$signif <- if_else(cd8_signif_df$LFSR < 0.01, TRUE, FALSE)
cd8_signif_df$cell <- "Memory CD8+ T cell"




final_df <- rbind(nk_signif_df, monocyte_signif_df, Non_classical_Monocytes_df,b_signif_df,cd8_signif_df)


final_df <- final_df %>% filter(LFSR < 0.01) %>% group_by(Timepoint, cell) %>% summarize(value = n()) %>% as.data.frame()

head(final_df)
final_df$cell <- as.factor(final_df$cell)
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(final_df$cell), ncol(final_df)) )
colnames(to_add) <- colnames(final_df)
to_add$cell <- rep(levels(final_df$cell), each=empty_bar)
final_df <- rbind(final_df, to_add)
final_df <- final_df %>% arrange(cell)
final_df$id <- seq(1, nrow(final_df))

final_df
# Get the name and the y position of each label
label_data <- final_df
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- final_df %>% 
  group_by(cell) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
# Make the plot
grid_data

my_palette = rep(c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026", "white", "white"), 5)
ggplot(final_df, aes(x=as.factor(id), y=value, fill=cell)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=cell), stat="identity", alpha=0.5, width = 0.9) +  theme_minimal() + ylim(-800, 800) +
   scale_fill_manual(values = c("#40004b", "#762a83", "#9970ab",  "#5aae61", "#1b7837", "#00441b")) +
   geom_text(data=label_data, aes(x=id, y=value + 10, label=value, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  theme(axis.text = element_blank(),
  axis.title = element_blank(),
  panel.grid = element_blank(),
  legend.position = "",
  plot.margin = unit(rep(-15,10), "cm")) +
  coord_polar(start = 0) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(final_df$id),4), y = c(25,50, 75, 100), label = c("25", "50", "75", "100") , color="grey", size=5 , angle=0, fontface="bold", hjust=1) +
  geom_text(data=label_data, aes(x=id, y=value+100, label=Timepoint, hjust=hjust), color=my_palette, fontface="bold",alpha=1, size=7, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start -0.5, y = -5.5, xend = end +0.5, yend = -5.5), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=cell), hjust=c(1,1,0.5,0, 0.5), vjust = c(4, 0, -3, 0, 1), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)
ggsave(args[6], width = 12, height = 12, dpi = 600) 



symbols <- fread(args[7]) %>% filter(V3 == "gene")
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)
symbols$gene_name <- gsub(" ", "", symbols$gene_name)
symbols <- symbols %>% select(gene_id, gene_name)


#T0_nk$variant_gene <- paste0(T0_nk$variant_id, "-", T0_nk$phenotype_id)
final_df <- rbind(nk_signif_df, cd8_signif_df, monocyte_signif_df, b_signif_df, Non_classical_Monocytes_df)
head(final_df)
final_df <- separate(final_df, col = Term, into = c("variant", "gene_id"), sep = "-")
final_df
final_df$Term <- paste0(final_df$variant, "-", final_df$gene_id)

final_df <- left_join(final_df, symbols)


open_eGene_signif = function(timepoints, cell_type = NULL){
  tensor_colnames_egene = c("phenotype_id",   "variant_id",     "start_distance", "af",             "ma_samples",     "ma_count",       
                            "pval_g",         "b_g",            "b_g_se",        "pval_i",         "b_i",            "b_i_se",         
                            "pval_gi",        "b_gi",           "b_gi_se",        "tests_emt",      "pval_emt",       "pval_adj_bh" )
  
  results_df_egene <- as.data.frame(matrix(ncol = 18, nrow = 0))
  
  colnames(results_df_egene) <- tensor_colnames_egene
  for (t in timepoints) {
    data_temp = read.table(paste0("/home/workspace/jogrady/heQTL/work/ieQTL/", t, "_", cell_type,".cis_qtl_top_assoc.txt.gz"), skip = 1)
    colnames(data_temp) <- tensor_colnames_egene
    data_temp$timepoint <- t
    results_df_egene <- rbind(results_df_egene, data_temp)
  }
  return(results_df_egene)
}

signif_eGenes_nk = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "NK_cell")
dim(signif_eGenes_nk)
signif_eGenes_nk$variant_gene <- paste0(signif_eGenes_nk$variant_id, "-", signif_eGenes_nk$phenotype_id)
nk_genes <- get_significant_results(nk_cell_mash$MASH, thresh = 1)
names(nk_genes)
signif_eGenes_nk <- signif_eGenes_nk %>% filter(variant_gene %in% names(nk_genes))



signif_eGenes_nk #<- signif_eGenes %>% filter(is_eGene == TRUE)


signif_variants_nk = signif_eGenes_nk %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.


head(nk_signif_df)
head(signif_variants_nk)
nk_signif_df_merged <- left_join(nk_signif_df, signif_variants_nk, by = c("Term" = "variant_gene"))





signif_eGenes_cd8 = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Memory_CD8+_T_cell")
dim(signif_eGenes_cd8)
signif_eGenes_cd8$variant_gene <- paste0(signif_eGenes_cd8$variant_id, "-", signif_eGenes_cd8$phenotype_id)
cd8_genes <- get_significant_results(memory_CD8_T_cell_mash$MASH, thresh = 1)
names(cd8_genes)
signif_eGenes_cd8 <- signif_eGenes_cd8 %>% filter(variant_gene %in% names(cd8_genes))

signif_variants_cd8 = signif_eGenes_cd8 %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
cd8_signif_df_merged <- left_join(cd8_signif_df, signif_variants_cd8, by = c("Term" = "variant_gene"))




signif_eGenes_cm = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Classical_Monocyte")
dim(signif_eGenes_cm)
signif_eGenes_cm$variant_gene <- paste0(signif_eGenes_cm$variant_id, "-", signif_eGenes_cm$phenotype_id)
cm_genes <- get_significant_results(Classical_Monocyte$MASH, thresh = 1)
names(cm_genes)
signif_eGenes_cm <- signif_eGenes_cm %>% filter(variant_gene %in% names(cm_genes))

signif_variants_cm = signif_eGenes_cm %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
cm_signif_df_merged <- left_join(monocyte_signif_df, signif_variants_cm, by = c("Term" = "variant_gene"))





signif_eGenes_ncm = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Non-classical_Monocyte")
dim(signif_eGenes_ncm)
signif_eGenes_ncm$variant_gene <- paste0(signif_eGenes_ncm$variant_id, "-", signif_eGenes_ncm$phenotype_id)
ncm_genes <- get_significant_results(Non_classical_Monocyte$MASH, thresh = 1)
names(ncm_genes)
signif_eGenes_ncm <- signif_eGenes_ncm %>% filter(variant_gene %in% names(ncm_genes))

signif_variants_ncm = signif_eGenes_ncm %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
ncm_signif_df_merged <- left_join(Non_classical_Monocytes_df, signif_variants_ncm, by = c("Term" = "variant_gene"))



signif_eGenes_nb = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Naïve_B_cell")
dim(signif_eGenes_nb)
signif_eGenes_nb$variant_gene <- paste0(signif_eGenes_nb$variant_id, "-", signif_eGenes_nb$phenotype_id)
nb_genes <- get_significant_results(Naïve_B_cell_mash$MASH, thresh = 1)
signif_eGenes_nb <- signif_eGenes_nb %>% filter(variant_gene %in% names(nb_genes))
signif_variants_nb = signif_eGenes_nb %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
nb_signif_df_merged <- left_join(b_signif_df, signif_variants_nb, by = c("Term" = "variant_gene"))



nk_signif_df_merged <- left_join(nk_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
cd8_signif_df_merged <- left_join(cd8_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
cm_signif_df_merged <- left_join(cm_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
ncm_signif_df_merged <- left_join(ncm_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
nb_signif_df_merged <- left_join(nb_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))



final_merged_df <- rbind(nk_signif_df_merged, cm_signif_df_merged, ncm_signif_df_merged,nb_signif_df_merged, cd8_signif_df_merged)
dim(final_merged_df)
head(final_merged_df)
final_merged_df_write <- final_merged_df %>% group_by(Term, cell) %>% pivot_wider(., names_from = Timepoint, values_from = c(LFSR), id_cols = c("Term", "variant_id", "phenotype_id", "cell", "af", "ma_samples", "ma_count", "gene_name")) #%>% fill(everything(), .direction = "down")

write.table(final_merged_df_write, file = args[8], sep = "\t", row.names = FALSE, quote = FALSE)


colnames(final_merged_df)[24] <- "Timepoint_most_significant"
final_merged_df <- final_merged_df %>% select(1,2,3,4,5,6,25,7:24)

final_merged_df_mtb_private <- final_merged_df %>%
  group_by(variant_id, gene_name, cell, phenotype_id, pval_gi) %>%  # Grouping by unique identifiers
  summarise(
    T0 = any(Timepoint == "T0" & LFSR < 0.01),
    T1 = any(Timepoint == "T1" & LFSR < 0.01),
    T2 = any(Timepoint == "T2" & LFSR < 0.01),
    T3 = any(Timepoint == "T3" & LFSR < 0.01),
    T4 = any(Timepoint == "T4" & LFSR < 0.01),
    Term = first(Term), 
    LFSR = first(LFSR)  # Keeping one LFSR value per variant
  ) %>%
  filter(T0 == TRUE & T1 == FALSE & T2 == FALSE & T3 == FALSE & T4 == FALSE) %>%  # Keep only T0-only variants
  ungroup()

# View result
print(final_merged_df_mtb_private)

write.table(final_merged_df_mtb_private, file = args[9], sep = "\t", row.names = FALSE, quote = FALSE)



final_merged_df_treat = final_merged_df %>%
  group_by(variant_id, gene_name, cell, phenotype_id, pval_gi) %>%
  summarise(
    T0 = any(Timepoint == "T0" & LFSR < 0.01),
    T1 = any(Timepoint == "T1" & LFSR < 0.01),
    T2 = any(Timepoint == "T2" & LFSR < 0.01),
    T3 = any(Timepoint == "T3" & LFSR < 0.01),
    T4 = any(Timepoint == "T4" & LFSR < 0.01),
    Term = first(Term),
    LFSR = first(LFSR)
  ) %>%
  filter(T0 == FALSE & (T1 == TRUE & T2 == TRUE & T3 == TRUE & T4 == TRUE)) %>%
  ungroup()

write.table(final_merged_df_treat, file = args[10], sep = "\t", row.names = FALSE, quote = FALSE)





genes_signif_in_treat <- final_merged_df %>% filter(signif == TRUE & Timepoint != "T0") %>% filter(Timepoint != "T0") %>% group_by(cell, gene_name) %>% summarise(Number = n()) %>% filter(Number >= 3) %>% select(gene_name) %>% as.vector()

head(genes_signif_in_treat)
final_merged_df %>% filter(final_merged_df$gene_name %in% genes_signif_in_treat$gene_name) 
library(ggpubr)
library(cowplot)

final_merged_df$Timepoint <- factor(final_merged_df$Timepoint)

final_merged_df_signif <- final_merged_df %>% filter(signif == TRUE)
dim(final_merged_df_signif)
# 1. Create the histogram plot
phist <- gghistogram(
  final_merged_df_signif, x = "start_distance", rug = FALSE,
  fill = "Timepoint", palette = c(my_palette), alpha = 0.3, bins = 19
) + facet_wrap(~cell)
phist
# 2. Create the density plot with y-axis on the right
# Remove x axis elements
pdensity <- ggdensity(
final_merged_df_signif, x = "start_distance",
  col = "Timepoint", palette = c(my_palette), alpha = 0
) + facet_wrap(~cell) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")
pdensity
# 3. Align the two plots and then overlay them.
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
ggsave(args[11], width = 12, height = 8, dpi = 600)

save.image(file=args[12])