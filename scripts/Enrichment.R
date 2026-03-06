# Script to see if eQTL/reQTL/ieQTL are enriched in epigenomic peaks

library(tidyverse)
library(ggplot2)
library(data.table)
library(regioneR)
library(vcfR)
library(BiocParallel)
options("mc.cores" = 50)
detectCores()


# Load in the variants we need to test for enrichment
# Including the background set of variants (all variants tested in the eQTL analysis)
load('/home/workspace/jogrady/heQTL/results/ieQTLs/ieQTL_MASHR.RData')
theloadedobjects <- load('/home/workspace/jogrady/heQTL/results/ieQTLs/ieQTL_MASHR.RData')
all_variants <- vcfR::read.vcfR('/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz')
eQTL <- fread('/home/workspace/jogrady/heQTL/results/reQTLs/MashR_new_LFSR_values_annotated.txt')
reQTL <- fread('/home/workspace/jogrady/heQTL/results/reQTLs/DR_eQTLs_MASH_2_lfsr_0.05.txt')
reQTL <- reQTL %>% separate(., col = x, into = c("variant", "gene"), sep = "-")



# All variants variant information store in @fix$ID column in the format chr:pos:ref:alt
# eQTL = variant information (chr:pos:ref:alt) stored in variant column
# reQTL = variant information (chr:pos:ref:alt) stored in variant column
# ieQTL = variant information stored in variant_id column

ieQTL <- final_merged_df %>% filter(signif == TRUE)

rm(theloadedobjects)

# Load in the regions from epigenomic data
ATAC_peaks <- fread('/home/workspace/jogrady/heQTL/data/epigenomic/Correa-macado-2020-ATAC-seq.txt') %>% select(1:7) %>% 
filter(Q.adj.P.Val_HC.Mtb < 0.05)


H3K27ac_peaks <- fread('/home/workspace/jogrady/heQTL/data/epigenomic/Correa-macado-2020-H3K27ac.txt') %>% select(c(1,11,12,2,3,4)) %>%
filter(Q_adj.P.Val_HC.Mtb < 0.05)
H3K27ac_peaks

table(sign(ATAC_peaks$log2FC_HC.Mtb))
table(sign(H3K27ac_peaks$logFC_HC.Mtb))

ATAC_peaks$PeakID <- gsub("chr", "", ATAC_peaks$PeakID)
ATAC_peaks <- ATAC_peaks %>% separate(., col = PeakID, into = c("chr", 'range'), sep = ":") 
ATAC_peaks <- ATAC_peaks %>% separate(., col = range, into = c("start", "end"), sep = "-")

H3K27ac_peaks$PeakID <- gsub("chr", "", H3K27ac_peaks$PeakID)
H3K27ac_peaks <- H3K27ac_peaks %>% separate(., col = PeakID, into = c("chr", 'start','end'), sep = "_")

# Now we have chr, start, end for each peak.

# ============================================================================
# Function to convert variant data to GRanges
# ============================================================================
extract_genomic_ranges <- function(variant_info, format_type = "chr:pos:ref:alt") {
  # Parse variant information and create GRanges objects
  if (format_type == "chr:pos:ref:alt") {
    variant_split <- str_split_fixed(variant_info, ":", 4)
    chr <- variant_split[, 1]
    pos <- as.numeric(variant_split[, 2])
    ref <- variant_split[, 3]
    alt <- variant_split[, 4]
  }
  
  # Create GRanges for variant position
  ranges <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = pos - 10000, end = pos + 10000)
  )
  return(ranges)
}

# Convert peaks to GRanges
ATAC_gr <- GRanges(
  seqnames = ATAC_peaks$chr,
  ranges = IRanges(start = as.numeric(ATAC_peaks$start), 
                   end = as.numeric(ATAC_peaks$end))
)
ATAC_gr
H3K27ac_gr <- GRanges(
  seqnames = H3K27ac_peaks$chr,
  ranges = IRanges(start = as.numeric(H3K27ac_peaks$start), 
                   end = as.numeric(H3K27ac_peaks$end))
)

# ============================================================================
# Function to perform enrichment analysis
# ============================================================================
perform_enrichment <- function(qtl_variants, background_variants, peaks, qtl_type, peak_type) {
  # Get the test set (QTL variants)
  test_ranges <- try(extract_genomic_ranges(qtl_variants))
  # Get the background set
  bg_ranges <- try(extract_genomic_ranges(background_variants))
  
    # Perform permutation test using regioneR
    pt <- permTest(
    A = test_ranges,
    B = peaks,
    genome = bg_ranges,
    randomize.function = randomizeRegions,
    evaluate.function = numOverlaps,
    ntimes = 1000,verbose = TRUE, force.parallel = TRUE, mc.cores = 50)
  
  return(list(
    qtl_type = qtl_type,
    peak_type = peak_type,
    permtest = pt
  ))
}

# ============================================================================
# Prepare variant data
# ============================================================================

cat("Preparing variant datasets...\n")

# Extract background variants (all tested variants)
all_var_ids <- all_variants@fix[, "ID"]

# eQTL variants
eqtl_variants <- eQTL$variant

# reQTL variants
reqtl_variants <- reQTL$variant

# ieQTL variants
ieqtl_variants <- ieQTL$variant_id

cat(sprintf("Total variants tested: %d\n", length(all_var_ids)))
cat(sprintf("eQTL variants: %d\n", length(eqtl_variants)))
cat(sprintf("reQTL variants: %d\n", length(reqtl_variants)))
cat(sprintf("ieQTL variants: %d\n", length(ieqtl_variants)))

# ============================================================================
# Perform enrichment tests
# ============================================================================

cat("\nPerforming enrichment analysis...\n")

# Test eQTL enrichment
cat("Testing eQTL enrichment in ATAC peaks...\n")
eqtl_atac <- perform_enrichment(eqtl_variants, all_var_ids, ATAC_gr, "eQTL", "ATAC")

cat("Testing eQTL enrichment in H3K27ac peaks...\n")
eqtl_h3k27ac <- perform_enrichment(eqtl_variants, all_var_ids, H3K27ac_gr, "eQTL", "H3K27ac")

# Test reQTL enrichment
cat("Testing reQTL enrichment in ATAC peaks...\n")
reqtl_atac <- perform_enrichment(reqtl_variants, all_var_ids, ATAC_gr, "reQTL", "ATAC")

cat("Testing reQTL enrichment in H3K27ac peaks...\n")
reqtl_h3k27ac <- perform_enrichment(reqtl_variants, all_var_ids, H3K27ac_gr, "reQTL", "H3K27ac")

# Test ieQTL enrichment
cat("Testing ieQTL enrichment in ATAC peaks...\n")
ieqtl_atac <- perform_enrichment(ieqtl_variants, all_var_ids, ATAC_gr, "ieQTL", "ATAC")

cat("Testing ieQTL enrichment in H3K27ac peaks...\n")
ieqtl_h3k27ac <- perform_enrichment(ieqtl_variants, all_var_ids, H3K27ac_gr, "ieQTL", "H3K27ac")

# ============================================================================
# Compile and summarize results
# ============================================================================

enrichment_results <- list(
  eqtl_atac = eqtl_atac,
  eqtl_h3k27ac = eqtl_h3k27ac,
  reqtl_atac = reqtl_atac,
  reqtl_h3k27ac = reqtl_h3k27ac,
  ieqtl_atac = ieqtl_atac,
  ieqtl_h3k27ac = ieqtl_h3k27ac
)
saveRDS(enrichment_results, file = '/home/workspace/jogrady/heQTL/results/Response_1/outputs/Enrichment_Results.rds')
enrichment_results$eqtl_h3k27ac$permtest
# Create summary table
enrichment_summary <- bind_rows(
  tibble(
    QTL_Type = "eQTL",
    Peak_Type = "ATAC",
    N_QTL = eqtl_atac$n_qtl,
    N_Peaks = eqtl_atac$n_peaks,
    Observed_Overlaps = eqtl_atac$observed,
    Expected_Overlaps = mean(eqtl_atac$permtest$numOverlaps),
    P_value = eqtl_atac$permtest$pval,
    Fold_Change = eqtl_atac$observed / mean(eqtl_atac$permtest$numOverlaps)
  ),
  tibble(
    QTL_Type = "eQTL",
    Peak_Type = "H3K27ac",
    N_QTL = eqtl_h3k27ac$n_qtl,
    N_Peaks = eqtl_h3k27ac$n_peaks,
    Observed_Overlaps = eqtl_h3k27ac$observed,
    Expected_Overlaps = mean(eqtl_h3k27ac$permtest$numOverlaps),
    P_value = eqtl_h3k27ac$permtest$pval,
    Fold_Change = eqtl_h3k27ac$observed / mean(eqtl_h3k27ac$permtest$numOverlaps)
  ),
  tibble(
    QTL_Type = "reQTL",
    Peak_Type = "ATAC",
    N_QTL = reqtl_atac$n_qtl,
    N_Peaks = reqtl_atac$n_peaks,
    Observed_Overlaps = reqtl_atac$observed,
    Expected_Overlaps = mean(reqtl_atac$permtest$numOverlaps),
    P_value = reqtl_atac$permtest$pval,
    Fold_Change = reqtl_atac$observed / mean(reqtl_atac$permtest$numOverlaps)
  ),
  tibble(
    QTL_Type = "reQTL",
    Peak_Type = "H3K27ac",
    N_QTL = reqtl_h3k27ac$n_qtl,
    N_Peaks = reqtl_h3k27ac$n_peaks,
    Observed_Overlaps = reqtl_h3k27ac$observed,
    Expected_Overlaps = mean(reqtl_h3k27ac$permtest$numOverlaps),
    P_value = reqtl_h3k27ac$permtest$pval,
    Fold_Change = reqtl_h3k27ac$observed / mean(reqtl_h3k27ac$permtest$numOverlaps)
  ),
  tibble(
    QTL_Type = "ieQTL",
    Peak_Type = "ATAC",
    N_QTL = ieqtl_atac$n_qtl,
    N_Peaks = ieqtl_atac$n_peaks,
    Observed_Overlaps = ieqtl_atac$observed,
    Expected_Overlaps = mean(ieqtl_atac$permtest$numOverlaps),
    P_value = ieqtl_atac$permtest$pval,
    Fold_Change = ieqtl_atac$observed / mean(ieqtl_atac$permtest$numOverlaps)
  ),
  tibble(
    QTL_Type = "ieQTL",
    Peak_Type = "H3K27ac",
    N_QTL = ieqtl_h3k27ac$n_qtl,
    N_Peaks = ieqtl_h3k27ac$n_peaks,
    Observed_Overlaps = ieqtl_h3k27ac$observed,
    Expected_Overlaps = mean(ieqtl_h3k27ac$permtest$numOverlaps),
    P_value = ieqtl_h3k27ac$permtest$pval,
    Fold_Change = ieqtl_h3k27ac$observed / mean(ieqtl_h3k27ac$permtest$numOverlaps)
  )
)

ieqtl_h3k27ac$permtest$numOverlaps$Z-score

# Add adjusted p-values
enrichment_summary <- enrichment_summary %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  mutate(Significant = ifelse(FDR < 0.05, "Yes", "No"))

cat("\n=== ENRICHMENT ANALYSIS RESULTS ===\n")
print(enrichment_summary)

# Save results
write.csv(enrichment_summary, 
          file = '/home/workspace/jogrady/heQTL/results/enrichment/Enrichment_Summary.csv',
          row.names = FALSE)

# ============================================================================
# Visualization
# ============================================================================

# Plot 1: Fold change comparison
enrichment_plot <- enrichment_summary %>%
  ggplot(aes(x = QTL_Type, y = Fold_Change, fill = Peak_Type)) +
  geom_col(position = "dodge") +
  labs(title = "Enrichment of QTLs in Epigenomic Peaks",
       x = "QTL Type",
       y = "Fold Change (Observed/Expected)",
       fill = "Peak Type") +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5)

ggsave(filename = '/home/workspace/jogrady/heQTL/results/enrichment/Enrichment_FoldChange.pdf',
       plot = enrichment_plot, width = 8, height = 6)

# Plot 2: Significance heatmap
enrichment_heatmap <- enrichment_summary %>%
  pivot_wider(names_from = Peak_Type, values_from = P_value) %>%
  column_to_rownames("QTL_Type") %>%
  select(-starts_with("N_"), -starts_with("Observed"), -starts_with("Expected"), -starts_with("Fold"), -starts_with("FDR"), -Significant) %>%
  as.matrix()

pdf(file = '/home/workspace/jogrady/heQTL/results/enrichment/Enrichment_PValues.pdf')
heatmap(enrichment_heatmap, main = "P-values for QTL Enrichment in Epigenomic Peaks")
dev.off()

cat("\n=== Plots saved ===\n")
cat("Fold change: /home/workspace/jogrady/heQTL/results/enrichment/Enrichment_FoldChange.pdf\n")
cat("P-values: /home/workspace/jogrady/heQTL/results/enrichment/Enrichment_PValues.pdf\n")


# Get the test set (QTL variants)
test_ranges <- try(extract_genomic_ranges(reqtl_variants))
if (class(test_ranges) == "try-error") {
cat(sprintf("Error processing %s variants\n", qtl_type))
return(NULL)
}

# Get the background set
bg_ranges <- try(extract_genomic_ranges(all_var_ids))
if (class(bg_ranges) == "try-error") {
cat(sprintf("Error processing background variants\n"))
return(NULL)
}
bg_ranges
# Perform permutation test using regioneR
pt <- permTest(
A = test_ranges,
B = ATAC_gr,
genome = bg_ranges,
randomize.function = randomizeRegions,
evaluate.function = numOverlaps,
ntimes = 1000,verbose = TRUE, force.parallel = TRUE, mc.cores = 50)

plot(pt)

return(list(
qtl_type = qtl_type,
peak_type = peak_type,
permtest = pt,
n_qtl = length(test_ranges),
n_peaks = length(peaks),
observed = length(subsetByOverlaps(test_ranges, peaks))
))

ieqtl_atac
