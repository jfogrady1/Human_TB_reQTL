suppressPackageStartupMessages({
  library(regioneR)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(readr)
  library(data.table)
  library(vcfR)
})


cis_eqtl_file      <- fread("/home/workspace/jogrady/heQTL/data/epigenomic/cis_eQTL_input.txt")
response_eqtl_file     <- fread("/home/workspace/jogrady/heQTL/data/epigenomic/reQTL_input.txt")
ieqtl_file         <- fread("/home/workspace/jogrady/heQTL/data/epigenomic/ieQTL_input.txt")
background_file    <- vcfR::read.vcfR("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz")
dcar_file          <- fread("/home/workspace/jogrady/heQTL/data/epigenomic/Correa-macado-2020-ATAC-seq.txt") %>% select(1:7) %>% filter(Q.adj.P.Val_HC.Mtb < 0.05)
h3k27ac_file       <- fread("/home/workspace/jogrady/heQTL/data/epigenomic/Correa-macado-2020-H3K27ac.txt") %>% select(1:4) %>% filter(Q_adj.P.Val_HC.Mtb < 0.05)

head(dcar_file)
head(background_file@fix)


background_GR <- as.data.frame(background_file@fix) %>%
  mutate(
    CHROM = gsub("chr", "", CHROM),
    START = as.numeric(POS) - 1,
    END = as.numeric(POS)
  ) %>%
  makeGRangesFromDataFrame(
    seqnames.field = "CHROM",
    start.field = "START",
    end.field = "END",
    keep.extra.columns = FALSE
  )


length(unique(cis_eqtl_file$variant)) # 5,287 variants unqiue for cis-eQTL

cis_eqtl_file <- cis_eqtl_file %>%
  separate(variant, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
  distinct(variant, .keep_all = TRUE) %>%
  mutate(
    chr = as.character(chr),
    start = as.numeric(pos) - 1,  # Convert to 0-based (BED format)
    end = as.numeric(pos)
  ) %>%
  select(chr, start, end, variant, gene_name) %>%
  arrange(as.numeric(chr), start)
head(cis_eqtl_file)

response_eqtl_file <- response_eqtl_file %>%
  separate(snp, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
  distinct(snp, .keep_all = TRUE) %>%
  mutate(
    chr = as.character(chr),
    start = as.numeric(pos) - 1,  # Convert to 0-based (BED format)
    end = as.numeric(pos)
  ) %>%
  select(chr, start, end, snp, gene_name) %>%
  arrange(as.numeric(chr), start)

dim(response_eqtl_file) # 790

ieqtl_file <- ieqtl_file %>%
  separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
  distinct(variant_id, .keep_all = TRUE) %>%
  mutate(
    chr = as.character(chr),
    start = as.numeric(pos) - 1,  # Convert to 0-based (BED format)
    end = as.numeric(pos)
  ) %>%
  select(chr, start, end, variant_id, gene_name) %>%
  arrange(as.numeric(chr), start)

dim(ieqtl_file) # unique variants.

eqtl_GR <- makeGRangesFromDataFrame(
  cis_eqtl_file,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

reqtl_GR <- makeGRangesFromDataFrame(
  response_eqtl_file,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)
ieQTL_GR <- makeGRangesFromDataFrame(
  ieqtl_file,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)



ATAC_GR <- dcar_file %>%
  separate(PeakID, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(
    chr = gsub("chr", "", chr),
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE)

head(h3k27ac_file)
H3K27AC_GR <- h3k27ac_file %>%
  separate(PeakID, into = c("chr", "start", "end"), sep = "_") %>%
  mutate(
    chr = gsub("chr", "", chr),
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE)
    


# Expand features by +/- 10Kb
expand_region <- function(gr, flank_bp = 5000) {
  return(GenomicRanges::resize(gr, width = width(gr) + 2*flank_bp, fix = "center"))
}

ATAC_expanded <- expand_region(ATAC_GR, flank_bp = 5000)
H3K27AC_expanded <- expand_region(H3K27AC_GR, flank_bp = 5000)
ATAC_expanded


seqlevels(eqtl_GR)
seqlevels(background_GR)
seqlevels(ATAC_expanded)
# Check for overlaps manually (should be > 0)
sum(countOverlaps(eqtl_GR, ATAC_expanded, ignore.strand=TRUE) > 0)


  # Manual enrichment test - more robust
set.seed(42)
run_perm_test <- function(query_GR, feature_GR, background_GR, ntimes = 1000, pseudocount = 0) {

  observed <- sum(countOverlaps(query_GR, feature_GR, ignore.strand = TRUE) > 0)
  print(observed)
  permuted <- replicate(ntimes, {
    shuffled <- background_GR[sample(length(background_GR), length(query_GR))]
    sum(countOverlaps(shuffled, feature_GR, ignore.strand = TRUE) > 0)
  })
  
  perm_mean <- mean(permuted)
  perm_sd <- sd(permuted)
  
  p_value <- (sum(permuted >= observed) + 1) / (ntimes + 1)
  z_score <- (observed - perm_mean) / perm_sd
  
  # Fold enrichment stats (stable)
  # Observed FE = observed / mean(permuted), with pseudocount to avoid divide-by-zero
  mean_fold_enrichment <- (observed + pseudocount) / (perm_mean + pseudocount)

  # Uncertainty around observed FE from permutation denominator distribution
  fold_enrichments <- (observed + pseudocount) / (permuted + pseudocount)
  print(fold_enrichments)
  se_fold_enrichment <- sd(fold_enrichments, na.rm = TRUE)
  ci_lower <- as.numeric(quantile(fold_enrichments, 0.025, na.rm = TRUE))
  ci_upper <- as.numeric(quantile(fold_enrichments, 0.975, na.rm = TRUE))
  
  return(list(
    observed = observed,
    perm_mean = perm_mean,
    perm_sd = perm_sd,
    p_value = p_value,
    z_score = z_score,
    fold_enrichment = mean_fold_enrichment,
    se_fold_enrichment = se_fold_enrichment,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    permuted_dist = permuted,
    fold_enrichment_dist = fold_enrichments
  ))
}


# Run enrichment tests for ATAC and H3K27ac only
n_perm <- 10000
qtl_sets <- list(
  "cis-eQTL" = eqtl_GR,
  "reQTL" = reqtl_GR,
  "ieQTL" = ieQTL_GR
)

feature_sets <- list(
  "ATAC" = ATAC_expanded,
  "H3K27ac" = H3K27AC_expanded
)

results_list <- list()

for (qtl_name in names(qtl_sets)) {
  for (feature_name in names(feature_sets)) {
    test_result <- run_perm_test(
      query_GR = qtl_sets[[qtl_name]],
      feature_GR = feature_sets[[feature_name]],
      background_GR = background_GR,
      ntimes = n_perm
    )

    key <- paste(qtl_name, feature_name, sep = "__")
    results_list[[key]] <- test_result
  }
}

names(results_list)
results_list[["cis-eQTL__ATAC"]][1] # 825
results_list[["cis-eQTL__H3K27ac"]][1] # 199
results_list[["reQTL__ATAC"]][1] # 122
results_list[["reQTL__H3K27ac"]][1] # 30
results_list[["ieQTL__ATAC"]][1] # 142
results_list[["ieQTL__H3K27ac"]][1] # 36


# Summary table for plot
enrichment_df <- bind_rows(lapply(names(results_list), function(key) {
  split_key <- strsplit(key, "__")[[1]]
  res <- results_list[[key]]

  data.frame(
    QTLType = split_key[1],
    Feature = split_key[2],
    FoldEnrichment = res$fold_enrichment,
    CILower = res$ci_lower,
    CIUpper = res$ci_upper,
    PValue = res$p_value,
    stringsAsFactors = FALSE
  )
}))

enrichment_df$Feature <- factor(enrichment_df$Feature, levels = c("ATAC", "H3K27ac"))
enrichment_df$QTLType <- factor(enrichment_df$QTLType, levels = c("cis-eQTL", "reQTL", "ieQTL"))

enrichment_df <- enrichment_df %>%
  mutate(
    Category = factor(paste(QTLType, Feature, sep = " | "),
                      levels = rev(paste(rep(c("cis-eQTL", "reQTL", "ieQTL"), each = 2),
                                         rep(c("ATAC", "H3K27ac"), times = 3), sep = " | "))),
    PLabel = paste0("p=", signif(PValue, 2))
  )

print(enrichment_df)

# Source data for the figure: summary values used in the plot
plot_source_df <- enrichment_df %>%
  select(QTLType, Feature, Category, FoldEnrichment, CILower, CIUpper, PValue, PLabel)

# Optional raw permutation data for each category
plot_permutation_df <- bind_rows(lapply(names(results_list), function(key) {
  split_key <- strsplit(key, "__")[[1]]
  res <- results_list[[key]]

  data.frame(
    QTLType = split_key[1],
    Feature = split_key[2],
    Category = paste(split_key[1], split_key[2], sep = " | "),
    Permutation = seq_along(res$permuted_dist),
    PermutedOverlap = as.integer(res$permuted_dist),
    PermutedFoldEnrichment = (res$observed) / (res$permuted_dist),
    stringsAsFactors = FALSE
  )
}))

enrichment_plot <- ggplot(enrichment_df, aes(y = Category, x = FoldEnrichment, fill = Feature)) +
  geom_col(width = 0.68, alpha = 0.95) +
  geom_errorbar(aes(xmin = CILower, xmax = CIUpper), width = 0.2, linewidth = 0.5) +
  geom_point(aes(x = FoldEnrichment), size = 2, color = "black") +
  geom_text(aes(x = pmax(FoldEnrichment, CIUpper) + 0.05, label = PLabel), size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey35") +
  scale_fill_manual(values = c("ATAC" = "#4DBBD5FF", "H3K27ac" = "#E64B35FF")) +
  labs(x = "Fold Enrichment", y = NULL,
       subtitle = "Bars: observed enrichment; error bars: 95% permutation interval") +
  theme_classic(base_size = 11) +
  theme(legend.title = element_blank(), legend.position = "right")

print(enrichment_plot)

write.table(plot_source_df,
            file = "/home/workspace/jogrady/heQTL/results/epigenomic_enrichment_plot_source_data.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(plot_permutation_df,
            file = "/home/workspace/jogrady/heQTL/results/epigenomic_enrichment_permutation_raw_values.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(enrichment_df,
            file = "/home/workspace/jogrady/heQTL/results/epigenomic_enrichment_summary.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

ggsave(
  filename = "/home/workspace/jogrady/heQTL/results/epigenomic_enrichment_horizontal_barplot.pdf",
  plot = enrichment_plot,
  width = 15,
  height = 12,
  dpi = 600
)
