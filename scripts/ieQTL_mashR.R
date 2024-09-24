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
  dim(signif_eGenes)
  signif_eGenes$variant_gene <- paste0(signif_eGenes$variant_id, "-", signif_eGenes$phenotype_id)
  signif_eGenes <- signif_eGenes %>% filter(variant_gene %in% common_variant_gene)

  signif_eGenes #<- signif_eGenes %>% filter(is_eGene == TRUE)
  dim(signif_eGenes)

  signif_eGenes$variant_gene <- paste0(signif_eGenes$variant_id,"-",signif_eGenes$phenotype_id)

  length(unique(signif_eGenes$phenotype_id)) # 1824

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

#signif_eGenes = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Classical_Monocyte")

#signif_eGenes %>% filter(phenotype_id == "ENSG00000167468.20") %>% View()
#min(signif_eGenes$pval_gi)
nk_cell_mash <- MASH("NK_cell")
memory_CD4_Tcell_mash <- MASH("Memory_CD4+_T_cell")
memory_Classical_Monocyte <- MASH("Classical_Monocyte")
memory_Naive_B_cell_mash <- MASH(celltype = "Naive_B_cell")
memory_CD8_T_cell_mash <- MASH("Memory_CD8+_T_cell")
Megakaryocyte <- MASH("Megakaryocyte")



length(get_significant_results(nk_cell_mash$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))

length(get_significant_results(memory_Classical_Monocyte$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))


length(get_significant_results(memory_Naive_B_cell_mash$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))


length(get_significant_results(Megakaryocyte$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))

length(get_significant_results(memory_CD4_Tcell_mash$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))
length(get_significant_results(memory_CD8_T_cell_mash$MASH, thresh = 0.01, conditions = NULL, sig_fn = get_lfsr))






library(tidyverse)
library(viridis)
nk_signif <- get_significant_results(nk_cell_mash$MASH, thresh = 0.01, conditions = NULL)
nk_signif <- names(nk_signif)
nk_signif <- nk_cell_mash$MASH$result$lfsr[nk_signif,]
nk_signif_df <- as.data.frame(nk_signif)
nk_signif_df$Term <- rownames(nk_signif)



nk_signif_df <- pivot_longer(nk_signif_df, cols = colnames(nk_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
nk_signif_df$signif <- if_else(nk_signif_df$LFSR < 0.01, TRUE, FALSE)
nk_signif_df$cell <- "NK cell"


cd4_signif <- get_significant_results(memory_CD4_Tcell_mash$MASH, thresh = 0.01, conditions = NULL)
cd4_signif <- names(cd4_signif)
cd4_signif <- memory_CD4_Tcell_mash$MASH$result$lfsr[cd4_signif,]

cd4_signif_df <- as.data.frame(cd4_signif)
cd4_signif_df$Term <- rownames(cd4_signif)



cd4_signif_df <- pivot_longer(cd4_signif_df, cols = colnames(cd4_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
cd4_signif_df$signif <- if_else(cd4_signif_df$LFSR < 0.01, TRUE, FALSE)
cd4_signif_df$cell <- "Memory CD4+ T cell"



# Monocyte 
monocyte_signif <- get_significant_results(memory_Classical_Monocyte$MASH, thresh = 0.01, conditions = NULL)
monocyte_signif <- names(monocyte_signif)
monocyte_signif <- memory_Classical_Monocyte$MASH$result$lfsr[monocyte_signif,]

monocyte_signif_df <- as.data.frame(monocyte_signif)
monocyte_signif_df$Term <- rownames(monocyte_signif)



monocyte_signif_df <- pivot_longer(monocyte_signif_df, cols = colnames(monocyte_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
monocyte_signif_df$signif <- if_else(monocyte_signif_df$LFSR < 0.01, TRUE, FALSE)
monocyte_signif_df$cell <- "Classical Monocyte"


# B cell
b_signif <- get_significant_results(memory_Naive_B_cell_mash$MASH, thresh = 0.01, conditions = NULL)
b_signif <- names(b_signif)
b_signif <- memory_Naive_B_cell_mash$MASH$result$lfsr[b_signif,]

b_signif_df <- as.data.frame(b_signif)
b_signif_df$Term <- rownames(b_signif)



b_signif_df <- pivot_longer(b_signif_df, cols = colnames(b_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
b_signif_df$signif <- if_else(b_signif_df$LFSR < 0.01, TRUE, FALSE)
b_signif_df$cell <- "Naive B Cell"


#Megakaryocyte
megakaryocytes <- get_significant_results(Megakaryocyte$MASH, thresh = 0.01, conditions = NULL)
megakaryocytes <- names(megakaryocytes)
megakaryocytes <- Megakaryocyte$MASH$result$lfsr[megakaryocytes,]

megakaryocytes_df <- as.data.frame(megakaryocytes)
megakaryocytes_df$Term <- rownames(megakaryocytes)



megakaryocytes_df <- pivot_longer(megakaryocytes_df, cols = colnames(megakaryocytes)[1:5], names_to = "Timepoint", values_to = "LFSR")
megakaryocytes_df$signif <- if_else(megakaryocytes_df$LFSR < 0.01, TRUE, FALSE)
megakaryocytes_df$cell <- "Megakaryocyte"

#CD8+ tcells


cd8_signif <- get_significant_results(memory_CD8_T_cell_mash$MASH, thresh = 0.01, conditions = NULL)
cd8_signif <- names(cd8_signif)
cd8_signif <- memory_CD8_T_cell_mash$MASH$result$lfsr[cd8_signif,]

cd8_signif_df <- as.data.frame(cd8_signif)
cd8_signif_df$Term <- rownames(cd8_signif)



cd8_signif_df <- pivot_longer(cd8_signif_df, cols = colnames(cd8_signif)[1:5], names_to = "Timepoint", values_to = "LFSR")
cd8_signif_df$signif <- if_else(cd8_signif_df$LFSR < 0.01, TRUE, FALSE)
cd8_signif_df$cell <- "Memory CD8+ T cell"




final_df <- rbind(nk_signif_df, monocyte_signif_df, megakaryocytes_df,b_signif_df, cd4_signif_df, cd8_signif_df)
dim(final_df)
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

my_palette = rep(c("#ffeda0", "#feb24c", "#fc4e2a", "#bd0026", "#800026", "white", "white"), 6)
ggplot(final_df, aes(x=as.factor(id), y=value, fill=cell)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=cell), stat="identity", alpha=0.5, width = 0.9) +  theme_minimal() + ylim(-800, 800) +
   scale_fill_manual(values = c("#40004b", "#762a83", "#9970ab",  "#5aae61", "#1b7837", "#00441b")) +
   geom_text(data=label_data, aes(x=id, y=value-50, label=value, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
  
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  theme(axis.text = element_blank(),
  axis.title = element_blank(),
  panel.grid = element_blank(),
  legend.position = "",
  plot.margin = unit(rep(-15,10), "cm")) +
  coord_polar(start = 0) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(final_df$id),3), y = c(50,100,150), label = c("50","100", "150") , color="grey", size=5 , angle=0, fontface="bold", hjust=1) +
  geom_text(data=label_data, aes(x=id, y=value+10, label=Timepoint, hjust=hjust), color=my_palette, fontface="bold",alpha=1, size=7, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start -0.5, y = -5.5, xend = end +0.5, yend = -5.5), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=cell), hjust=c(1,1,0.5,0, 0.5, 0.5), vjust = c(4, 0, -3, 0, 1, 1), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)
ggsave("/home/workspace/jogrady/heQTL/work/ieQTL/MASHR_LFSR_0.01_6_Cell_types.pdf", width = 12, height = 12, dpi = 600) 



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


#T0_nk$variant_gene <- paste0(T0_nk$variant_id, "-", T0_nk$phenotype_id)
final_df <- rbind(nk_signif_df, cd4_signif_df, cd8_signif_df, monocyte_signif_df, b_signif_df, megakaryocytes_df)
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


nk_signif_df_merged <- left_join(nk_signif_df, signif_variants_nk, by = c("Term" = "variant_gene"))
View(nk_signif_df_merged)




signif_eGenes_cd8 = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Memory_CD8+_T_cell")
dim(signif_eGenes_cd8)
signif_eGenes_cd8$variant_gene <- paste0(signif_eGenes_cd8$variant_id, "-", signif_eGenes_cd8$phenotype_id)
cd8_genes <- get_significant_results(memory_CD8_T_cell_mash$MASH, thresh = 1)
names(cd8_genes)
signif_eGenes_cd8 <- signif_eGenes_cd8 %>% filter(variant_gene %in% names(cd8_genes))

signif_variants_cd8 = signif_eGenes_cd8 %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
cd8_signif_df_merged <- left_join(cd8_signif_df, signif_variants_cd8, by = c("Term" = "variant_gene"))


signif_eGenes_cd4 = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Memory_CD4+_T_cell")
dim(signif_eGenes_cd4)
signif_eGenes_cd4$variant_gene <- paste0(signif_eGenes_cd4$variant_id, "-", signif_eGenes_cd4$phenotype_id)
cd4_genes <- get_significant_results(memory_CD4_Tcell_mash$MASH, thresh = 1)
names(cd4_genes)
signif_eGenes_cd4 <- signif_eGenes_cd4 %>% filter(variant_gene %in% names(cd4_genes))

signif_variants_cd4 = signif_eGenes_cd4 %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
cd4_signif_df_merged <- left_join(cd4_signif_df, signif_variants_cd4, by = c("Term" = "variant_gene"))
View(cd4_signif_df_merged)


signif_eGenes_cm = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Classical_Monocyte")
dim(signif_eGenes_cm)
signif_eGenes_cm$variant_gene <- paste0(signif_eGenes_cm$variant_id, "-", signif_eGenes_cm$phenotype_id)
cm_genes <- get_significant_results(memory_Classical_Monocyte$MASH, thresh = 1)
names(cm_genes)
signif_eGenes_cm <- signif_eGenes_cm %>% filter(variant_gene %in% names(cm_genes))

signif_variants_cm = signif_eGenes_cm %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
cm_signif_df_merged <- left_join(monocyte_signif_df, signif_variants_cm, by = c("Term" = "variant_gene"))





signif_eGenes_mk = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Megakaryocyte")
dim(signif_eGenes_mk)
signif_eGenes_mk$variant_gene <- paste0(signif_eGenes_mk$variant_id, "-", signif_eGenes_mk$phenotype_id)
mk_genes <- get_significant_results(Megakaryocyte$MASH, thresh = 1)
names(mk_genes)
signif_eGenes_mk <- signif_eGenes_mk %>% filter(variant_gene %in% names(mk_genes))

signif_variants_mk = signif_eGenes_mk %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
mk_signif_df_merged <- left_join(megakaryocytes_df, signif_variants_mk, by = c("Term" = "variant_gene"))



signif_eGenes_nb = open_eGene_signif(timepoints = c("T0", "T1", "T2", "T3", "T4"), cell_type = "Naive_B_cell")
dim(signif_eGenes_nb)
signif_eGenes_nb$variant_gene <- paste0(signif_eGenes_nb$variant_id, "-", signif_eGenes_nb$phenotype_id)
nb_genes <- get_significant_results(memory_Naive_B_cell_mash$MASH, thresh = 1)
names(nb_genes)
signif_eGenes_nb <- signif_eGenes_nb %>% filter(variant_gene %in% names(nb_genes))

signif_variants_nb = signif_eGenes_nb %>% arrange(phenotype_id) %>% group_by(phenotype_id) %>% filter(pval_emt == min(pval_emt)) %>% group_by(phenotype_id) %>% filter(pval_gi == min(pval_gi)) # select most significant variant per gene and if there is a tie, select variant with lowest nominal association.
nb_signif_df_merged <- left_join(b_signif_df, signif_variants_nb, by = c("Term" = "variant_gene"))

head(nb_signif_df_merged)


nk_signif_df_merged <- left_join(nk_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
cd4_signif_df_merged <- left_join(cd4_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
cd8_signif_df_merged <- left_join(cd8_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
cm_signif_df_merged <- left_join(cm_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
mk_signif_df_merged <- left_join(mk_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))
nb_signif_df_merged <- left_join(nb_signif_df_merged, symbols, by = c("phenotype_id" = "gene_id"))



final_merged_df <- rbind(nk_signif_df_merged, cm_signif_df_merged, mk_signif_df_merged,nb_signif_df_merged, cd4_signif_df_merged, cd8_signif_df_merged)
tail(final_merged_df)
colnames(final_merged_df)[24] <- "Timepoint_most_significant"
final_merged_df <- final_merged_df %>% select(1,2,3,4,5,6,25,7:24)

colnames(final_merged_df)

final_merged_df_mtb_private <- final_merged_df %>% group_by(gene_name, cell, pval_gi, LFSR) %>% select(c(1:10)) %>% summarise(
  T0 = if_else(Timepoint == "T0" & LFSR < 0.01, TRUE, FALSE),
  T1 = if_else(Timepoint == "T1" & LFSR < 0.01, TRUE, FALSE),
  T2 = if_else(Timepoint == "T2" & LFSR < 0.01, TRUE, FALSE),
  T3 = if_else(Timepoint == "T3" & LFSR < 0.01, TRUE, FALSE),
  T4 = if_else(Timepoint == "T4" & LFSR < 0.01, TRUE, FALSE),
  Term = Term, variant_id = variant_id, phenotype_id = phenotype_id, pval_gi = pval_gi, LFSR = LFSR) %>% filter(T0 == TRUE & T1 == FALSE & T2 == FALSE & T3 == FALSE & T4 == FALSE)


write.table(final_merged_df_mtb_private, file = "/home/workspace/jogrady/heQTL/results/ieQTLs/mTB_private_ieQTLs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
View(final_merged_df_mtb_private)
final_merged_df %>% filter(gene_name == "ARID3A") %>% View()

final_merged_df_treat_private <- final_merged_df %>% group_by(gene_name, cell) %>% summarise(
  T0 = if_else(Timepoint == "T0" & LFSR < 0.01, TRUE, FALSE),
  T1 = if_else(Timepoint == "T1" & LFSR < 0.01, TRUE, FALSE),
  T2 = if_else(Timepoint == "T2" & LFSR < 0.01, TRUE, FALSE),
  T3 = if_else(Timepoint == "T3" & LFSR < 0.01, TRUE, FALSE),
  T4 = if_else(Timepoint == "T4" & LFSR < 0.01, TRUE, FALSE)) %>% filter(T0 == FALSE)
View(final_merged_df_treat_private)

View(final_merged_df)

final_merged_df %>% filter(LFSR < 0.01, Timepoint != "T0" %>% )


View(final_merged_df)
genes_signif_in_treat <- final_merged_df %>% filter(signif == TRUE & Timepoint != "T0") %>% filter(Timepoint == "T0") %>% group_by(cell, gene_name) %>% summarise(Number = n()) %>% filter(Number >= 3) %>% select(gene_name) %>% as.vector()
final_merged_df %>% View()

head(genes_signif_in_treat)
final_merged_df %>% filter(final_merged_df$gene_name %in% genes_signif_in_treat$gene_name) %>% View()
library(ggpubr)
library(cowplot)
View(final_merged_df)
head(final_merged_df)

final_merged_df$Timepoint <- factor(final_merged_df$Timepoint)

final_merged_df_signif <- final_merged_df %>% filter(signif == TRUE)
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
ggsave("/home/workspace/jogrady/heQTL/results/ieQTLs/Genomic_architecture_ieQTLs.pdf", width = 12, height = 8, dpi = 600)


# Lets look at TFR2 in classical monocytes
# variant = 7:100615131:C:T
# gene = ENSG00000106327.13




# VCF data 
library(vcfR)
vcf <- vcfR::read.vcfR("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz")
vcf@gt
vcf <- cbind(vcf@fix, vcf@gt)
head(vcf)
vcf <- as.data.frame(vcf)

head(T0_TFR2)
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

NLRC4_CM <- ieQTL_plot(variant = "2:32236868:A:AAAAG", gene = "ENSG00000091106.19", gene_name = "NLRC4", hom_ref = "A", het = "A:AAAAG", hom_alt = "AAAAGA:AAAAG", cell_type = "Classical_Monocyte", groups = c("T0", "T1", "T2", "T3", "T4"))
NLRC4_CM



GPX4_CM <- ieQTL_plot(variant = "19:1096925:G:GA", gene = "ENSG00000167468.20", gene_name = "GPX4", hom_ref = "G", het = "G:GA", hom_alt = "GA:GA", cell_type = "Classical_Monocyte", groups = c("T0", "T1", "T2", "T3", "T4"))
GPX4_CM


ggsave(GPX4_CM, file = "/home/workspace/jogrady/heQTL/results/ieQTLs/GPX4_Classical_Monocyte_ieQTL.pdf", width = 12, height = 8, dpi = 600)



IFNGR2_CM <- ieQTL_plot(variant = "21:33353924:T:C", gene = "ENSG00000159128.16", gene_name = "IFNGR2", hom_ref = "T", het = "T:C", hom_alt = "C:C", cell_type = "Classical_Monocyte", groups = c("T0", "T1", "T2", "T3", "T4"))
IFNGR2_CM

ggsave("/home/workspace/jogrady/heQTL/results/ieQTLs/IFNGR2_Classical_Monocyte_ieQTL.pdf", width = 12, height = 8, dpi = 600)

CD53_NK <- ieQTL_plot(variant = "1:110899907:A:AAC", gene = "ENSG00000167468.20", gene_name = "CD53", hom_ref = "A", het = "A:AAC", hom_alt = "AAC:AAC", cell_type = "NK_cell", groups = c("T0", "T1", "T2", "T3", "T4"))
CD53_NK






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
