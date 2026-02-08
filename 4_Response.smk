rule all:
    input:
        "/home/workspace/jogrady/heQTL/results/Response_1/Figures/lung_gtex_comparison.pdf",
        '/home/workspace/jogrady/heQTL/results/Response_1/Figures/individual_correlation_heatmap.pdf'

rule lung_GTEX_comparison:
    input:
        LFSR_results = "/home/workspace/jogrady/heQTL/results/reQTLs/MashR_new_LFSR_values_annotated.txt"
        GTEX_results = "/home/workspace/jogrady/heQTL/data/gtex/analysis/gtex_summary_statistics/GTEx_Analysis_v11_eQTL/Lung.v11.eQTLs.signif_pairs.parquet"
        script = "/home/workspace/jogrady/heQTL/scripts/lung_gtex_comparison.R"
    output:
        pdf =  "/home/workspace/jogrady/heQTL/results/Response_1/Figures/lung_gtex_comparison.pdf"
    shell:
        """
         Rscript {input.script} {input.GTEX_results} {input.LFSR_results}{output.pdf}
        """

correlation_between_covariates:
    input:
        script = "/home/workspace/jogrady/heQTL/scripts/covariate_correlation.R",
        genotype_pcs = '/home/workspace/jogrady/heQTL/data/covariate/Genotypes.PCA_eigenvect.txt'
    output:
        pdf = '/home/workspace/jogrady/heQTL/results/Response_1/Figures/individual_correlation_heatmap.pdf',

    params:
        counts = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_'
        sample = ['T0', 'T1', 'T2', 'T3', 'T4'],
        cov_path = '/home/workspace/jogrady/heQTL/data/covariate/'
    shell:
        """
        Rscript {input.script} {params.counts} {params.cov_path} {input.genotype_pcs} {params.sample} {output.pdf}
        """