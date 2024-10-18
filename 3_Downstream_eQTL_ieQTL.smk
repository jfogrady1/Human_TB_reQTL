# Main Snakefile
configfile: "/home/workspace/jogrady/heQTL/config/Alignment.yaml"  # Primary config file


autosome = [str(i) for i in range(1, 23)]

include: "snakepit/downstream/quantification.smk"
include: "snakepit/downstream/eqtl.smk"
include: "snakepit/downstream/response_eQTL.smk"
include: "snakepit/downstream/scRNA_deconvolution.smk"
include: "snakepit/downstream/ieQTL.smk"



rule all:
    input:
        #quantification with featureCounts
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/align/T0/{sample0}.bam", sample0 = config["fastqfiles_0"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/{sample1}.bam", sample1 = config["fastqfiles_1"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/align/T2/{sample2}.bam", sample2 = config["fastqfiles_2"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/align/T3/{sample3}.bam", sample3 = config["fastqfiles_3"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/align/T4/{sample4}.bam", sample4 = config["fastqfiles_4"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_{time}.txt", time = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt", timepoint = config["TIME"]),
        "/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_ALL.txt",
        # eQTL analysis
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt", timepoint = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_TPM_{timepoint}.txt", timepoint = config["TIME"]),
        "/home/workspace/jogrady/heQTL/data/covariate/Genotypes.PCA_eigenvect.txt",
        "/home/workspace/jogrady/heQTL/data/ref_genome/TSS_annotation.bed.gz",
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz", timepoint = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/data/covariate/Scree_plot_PCA_{timepoint}.pdf", timepoint = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs.txt", timepoint = config["TIME"]),
        multiext("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/eQTL_genotypes_R0.6", ".bed", ".bim", ".fam"),          
        expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl.txt.gz", timepoint = config["TIME"], window = config["cis_window"]),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl_fdr0.1.txt", window = config["cis_window"], timepoint = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl_pairs.chr{chr}.parquet", timepoint = config["TIME"], window = config["cis_window"], chr = autosome),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.{window}.cis_qtl_pairs.chr{chr}.txt.gz", timepoint = config["TIME"], window = config["cis_window"], chr = autosome),
        expand("/home/workspace/jogrady/heQTL/results/eQTL/{timepoint}.50000.cis_qtl_pairs.ALL.txt", timepoint = config["TIME"]),
        "/home/workspace/jogrady/heQTL/results/eQTL/numbers_stats_cis.txt",
        # response eQTL analysis
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Normalised_counts_{timepoint}_qtltools.bed.gz", timepoint = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/data/covariate/{timepoint}_eqtlcovs_qtltools.txt", timepoint = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/{timepoint}_residualised_expression.txt", timepoint = config["TIME"]),
        "/home/workspace/jogrady/heQTL/results/reQTLs/DR_eQTLs_MASH_2_lfsr_0.05.txt",
        # scRNA_seq work
        expand("/home/workspace/jogrady/heQTL/work/scRNA_seq/{singlecell_ids}", singlecell_ids=config["sc_ids"]),
        "/home/workspace/jogrady/heQTL/work/scRNA_seq/TB.combined.metadata.txt",
        "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/ALL_benchmark_results.txt",
        "/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/Raw_count_matrix_all_samples.txt",
        "/home/workspace/jogrady/heQTL/work/scRNA_seq/Benchmark/CIBERSORTx_Music_correlation_res.txt",
        # ieQTLs
        expand("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_{timepoint}.bed.gz", timepoint = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_NK_cell_interaction_input_from_cibersort.txt", t = config["TIME"], celltype = config["cell"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_NK_cell.cis_qtl_top_assoc.txt.gz", t = config["TIME"], celltype = config["cell"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Memory_CD8+_T_cell.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Naïve_B_cell.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Classical_Monocyte.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_Non-classical_Monocyte.cis_qtl_top_assoc.txt.gz", t = config["TIME"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_{celltype}_interaction_input_from_cibersort.txt", t = config["TIME"], celltype = config["cell"]),
        expand("/home/workspace/jogrady/heQTL/work/ieQTL/{t}_cell_interaction_input_from_cibersort.txt", t = config["TIME"]),
        "/home/workspace/jogrady/heQTL/results/ieQTLs/ieQTL_MASHR.RData",