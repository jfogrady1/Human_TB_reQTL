rule all:
    input:
        "/home/workspace/jogrady/heQTL/results/Response_1/Figures/lung_gtex_comparison.pdf"

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